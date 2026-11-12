"""Compile a Boolean spec for a given delivery format into a circuit design."""

from __future__ import annotations

from dataclasses import dataclass, field

from . import leakage, report as report_mod
from .architecture import Design
from .delivery import DeliveryVerdict, evaluate_all, minimum_vectors_needed, vector_catalogue
from .logic import BooleanFunction, Cover, parse_expression, parse_table
from .parts import Parts, Payload
from .synthesis import (
    Candidate,
    CircuitSpec,
    PartChoice,
    SynthesisError,
    best_effort,
    best_fitting,
    build_design,
    pareto_front,
    search,
)


@dataclass
class CompileResult:
    """The chosen design, its delivery verdicts and its leakage estimates."""

    spec: CircuitSpec
    design: Design
    covers: dict[str, Cover]
    choice: PartChoice
    target_vector: str
    target_verdict: DeliveryVerdict
    all_verdicts: list[DeliveryVerdict]
    leakage: dict[str, leakage.LeakageReport]
    candidates: list[Candidate] = field(default_factory=list)
    parts: Parts = field(default_factory=Parts)
    searched: bool = False

    @property
    def feasible(self) -> bool:
        return self.target_verdict.fits

    @property
    def circuit_bp(self) -> int:
        return self.design.total_bp

    @property
    def pareto(self) -> list[Candidate]:
        return pareto_front(self.candidates) if self.candidates else []

    def report(self, width: int = 80) -> str:
        return report_mod.render(self, width=width)

    def to_genbank(self, name: str = "cargologic_design") -> str:
        return self.design.to_genbank(name)

    def __str__(self) -> str:  # pragma: no cover
        return self.report()


def _coerce_functions(
    logic: str | dict[str, str],
    inputs: list[str] | None,
    as_table: bool,
) -> list[BooleanFunction]:
    parse = parse_table if as_table else parse_expression
    if isinstance(logic, str):
        return [parse(logic, name="out", inputs=inputs)]

    names = list(inputs) if inputs else None
    if names is None:
        # All outputs must share one input list.
        names = []
        for expr in logic.values():
            for v in parse(expr, name="tmp").inputs:
                if v not in names:
                    names.append(v)
    return [parse(expr, name=nm, inputs=names) for nm, expr in logic.items()]


def _coerce_payloads(
    payload: str | int | dict[str, str | int], functions: list[BooleanFunction], parts: Parts
) -> dict[str, Payload]:
    if isinstance(payload, dict):
        missing = {f.name for f in functions} - set(payload)
        if missing:
            raise SynthesisError(f"no payload given for output(s) {sorted(missing)}")
        return {k: parts.resolve_payload(v) for k, v in payload.items()}
    resolved = parts.resolve_payload(payload)
    return {f.name: resolved for f in functions}


def compile_circuit(
    logic: str | dict[str, str],
    payload: str | int | dict[str, str | int] = "EGFP",
    vector: str = "aav_single",
    inputs: list[str] | None = None,
    recombinases: list[str] | None = None,
    include_recombinases: bool = False,
    optimise: bool = True,
    min_dynamic_range: float = 10.0,
    parts: Parts | None = None,
    as_table: bool = False,
    promoter: str | None = None,
    blocker: str | None = None,
    terminator: str | None = None,
    n_samples: int = 2000,
) -> CompileResult:
    """Compile a Boolean specification into a recombinase circuit design."""
    parts = parts or Parts()
    cat = vector_catalogue()
    if vector not in cat:
        raise KeyError(f"unknown vector {vector!r}; choose from {sorted(cat)}")
    vec = cat[vector]

    functions = _coerce_functions(logic, inputs, as_table)
    spec = CircuitSpec(
        inputs=functions[0].inputs,
        functions=functions,
        payloads=_coerce_payloads(payload, functions, parts),
        recombinase_names=recombinases,
    )
    spec.validate()

    covers = {f.name: f.minimise() for f in functions}

    candidates: list[Candidate] = []
    searched = False

    if optimise and not (promoter or blocker or terminator):
        candidates = search(
            spec,
            parts,
            vec,
            include_recombinases=include_recombinases,
            min_dynamic_range=min_dynamic_range,
            n_samples=max(200, n_samples // 5),
        )
        searched = True
        pick = best_fitting(candidates, min_dynamic_range) or best_effort(candidates)
        choice, design = pick.choice, pick.design
    else:
        choice = PartChoice()
        if promoter:
            choice.promoter = promoter
        if blocker:
            choice.blocker = blocker
        if terminator:
            choice.terminator = terminator
        design = build_design(spec, parts, choice, include_recombinases)

    leak = {
        f.name: leakage.predict(design, f.name, n_samples=n_samples, seed=7)
        for f in functions
    }

    return CompileResult(
        spec=spec,
        design=design,
        covers=covers,
        choice=choice,
        target_vector=vector,
        target_verdict=vec.evaluate(design.total_bp),
        all_verdicts=evaluate_all(design.total_bp),
        leakage=leak,
        candidates=candidates,
        parts=parts,
        searched=searched,
    )


def minimum_aav_count(result: CompileResult) -> int:
    return minimum_vectors_needed(result.circuit_bp)
