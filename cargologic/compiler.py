"""Top-level entry point: a Boolean spec and a delivery format in, a verdict out."""

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
    """Everything the compiler worked out, plus a rendered report."""

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

    # Multi-output: parse each, then reconcile the input list across outputs.
    names = list(inputs) if inputs else None
    if names is None:
        seen: list[str] = []
        for expr in logic.values():
            f = parse(expr, name="tmp")
            for v in f.inputs:
                if v not in seen:
                    seen.append(v)
        names = seen
    return [parse(expr, name=nm, inputs=names) for nm, expr in logic.items()]


def _coerce_payloads(
    payload: str | int | dict[str, str | int], functions: list[BooleanFunction], parts: Parts
) -> dict[str, Payload]:
    if isinstance(payload, dict):
        missing = {f.name for f in functions} - set(payload)
        if missing:
            raise SynthesisError(
                f"no payload given for output(s) {sorted(missing)}"
            )
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
    """Compile a Boolean specification into a recombinase circuit design.

    Parameters
    ----------
    logic
        A Boolean expression (``"A and not B"``), a named gate (``"AND"``, ``"XOR"``,
        ``"MAJORITY3"``), or a dict mapping output names to expressions for a
        multi-output circuit. With ``as_table=True`` these are parsed as truth tables
        instead (``"00:0,01:1,10:1,11:0"``, ``"0110"`` or a minterm list ``"1,2"``).
    payload
        A library name (``"CAR_2g"``), a raw bp count, or a dict per output.
    vector
        Target delivery format. See ``cargologic.delivery.vector_catalogue()``.
    include_recombinases
        Count the recombinase coding sequences against the packaging budget. Usually
        false: published devices supply them separately.
    optimise
        Search part choices for a design that fits. Set false to use the defaults as-is.

    Returns
    -------
    CompileResult
        Call ``.report()`` for the human-readable verdict.

    Examples
    --------
    >>> r = compile_circuit("A and not B", payload="CAR_2g", vector="aav_single")
    >>> r.feasible
    True
    """
    parts = parts or Parts()
    cat = vector_catalogue()
    if vector not in cat:
        raise KeyError(f"unknown vector {vector!r}; choose from {sorted(cat)}")
    vec = cat[vector]

    functions = _coerce_functions(logic, inputs, as_table)
    payloads = _coerce_payloads(payload, functions, parts)
    spec = CircuitSpec(
        inputs=functions[0].inputs,
        functions=functions,
        payloads=payloads,
        recombinase_names=recombinases,
    )
    spec.validate()

    covers = {f.name: f.minimise() for f in functions}

    explicit = PartChoice()
    if promoter:
        explicit.promoter = promoter
    if blocker:
        explicit.blocker = blocker
    if terminator:
        explicit.terminator = terminator
    user_fixed = any([promoter, blocker, terminator])

    candidates: list[Candidate] = []
    searched = False

    if optimise and not user_fixed:
        candidates = search(
            spec,
            parts,
            vec,
            include_recombinases=include_recombinases,
            min_dynamic_range=min_dynamic_range,
            n_samples=max(200, n_samples // 5),
        )
        searched = True
        pick = best_fitting(candidates, min_dynamic_range)
        if pick is None:
            # Nothing both fits and clears the dynamic-range floor. Take the best
            # dynamic range available, breaking ties-within-tolerance by size.
            pick = best_effort(candidates)
        choice = pick.choice
        design = pick.design
    else:
        choice = explicit
        design = build_design(spec, parts, choice, include_recombinases)

    target_verdict = vec.evaluate(design.total_bp)
    all_verdicts = evaluate_all(design.total_bp)

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
        target_verdict=target_verdict,
        all_verdicts=all_verdicts,
        leakage=leak,
        candidates=candidates,
        parts=parts,
        searched=searched,
    )


def minimum_aav_count(result: CompileResult) -> int:
    """How many AAVs this design would need."""
    return minimum_vectors_needed(result.circuit_bp)
