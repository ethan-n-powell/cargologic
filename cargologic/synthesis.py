"""Boolean spec to DNA design, plus a search over part choices for one that fits."""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

from . import leakage
from .architecture import (
    build_blade,
    Design,
    build_cassette,
    build_recombinase_elements,
    build_shared_payload_tail,
    can_share_payload,
)
from .delivery import Vector, vector_catalogue
from .logic import BooleanFunction, Cover
from .parts import (
    DEFAULT_BLOCKER,
    DEFAULT_LINKER,
    DEFAULT_PROMOTER,
    DEFAULT_TERMINATOR,
    PROMOTER_DOWNSIZE_ORDER,
    Parts,
    Payload,
)


class SynthesisError(ValueError):
    pass


@dataclass
class CircuitSpec:
    inputs: list[str]
    functions: list[BooleanFunction]
    payloads: dict[str, Payload]
    recombinase_names: list[str] | None = None

    @property
    def n(self) -> int:
        return len(self.inputs)

    def validate(self) -> None:
        if not self.functions:
            raise SynthesisError("no outputs specified")
        for f in self.functions:
            if f.inputs != self.inputs:
                raise SynthesisError(
                    f"output {f.name!r} is over inputs {f.inputs} but the circuit "
                    f"declares {self.inputs}"
                )
            if f.name not in self.payloads:
                raise SynthesisError(f"no payload given for output {f.name!r}")


@dataclass
class PartChoice:
    """The part selections the search explores."""

    promoter: str = DEFAULT_PROMOTER
    blocker: str = DEFAULT_BLOCKER
    terminator: str = DEFAULT_TERMINATOR
    linker: str = DEFAULT_LINKER

    def key(self) -> tuple[str, str, str]:
        return (self.promoter, self.blocker, self.terminator)


def _group_outputs(
    spec: CircuitSpec, covers: dict[str, Cover]
) -> list[tuple[list[str], Cover]]:
    """Outputs with identical minimised covers share one cassette via a 2A linker."""
    groups: dict[tuple, list[str]] = {}
    for f in spec.functions:
        sig = tuple(sorted(t.pattern() for t in covers[f.name].terms))
        groups.setdefault(sig, []).append(f.name)
    return [(names, covers[names[0]]) for names in groups.values()]



def _try_blade(
    spec: CircuitSpec,
    group,
    promoter,
    terminator,
    recs,
    linker,
    include_recombinases: bool,
    parts: Parts,
    choice: PartChoice,
) -> Design | None:
    """Build the BLADE template for a single output group, or None if inapplicable."""
    names, cover = group
    if not cover.terms:
        return None
    on_minterms = {m for t in cover.terms for m in t.minterms()}
    if not on_minterms:
        return None

    payloads = [spec.payloads[nm] for nm in names]
    cassettes, ordered = build_blade(
        n_inputs=spec.n,
        on_minterms=on_minterms,
        outputs=names,
        payloads=payloads,
        promoter=promoter,
        terminator=terminator,
        recombinases=recs,
        linker=linker if len(payloads) > 1 else None,
    )

    rec_elements = []
    notes = [
        f"BLADE template: one promoter upstream of {1 << spec.n} addresses with "
        f"{(1 << spec.n) - 1} recombination site pairs. Chosen because the cover "
        f"does not factor well enough for product terms to pay off."
    ]
    if include_recombinases:
        rec_elements = build_recombinase_elements(recs, promoter, terminator, linker)
        notes.append(
            "Recombinase coding sequences are included in the budget. Published devices "
            "almost always supply them separately, and carrying them is what usually "
            "makes a design exceed AAV."
        )

    return Design(
        inputs=list(spec.inputs),
        recombinases=recs,
        cassettes=cassettes,
        strategy="blade",
        recombinase_elements=rec_elements,
        include_recombinases=include_recombinases,
        notes=notes,
        ordered_elements=ordered,
    )


def _blade_is_smaller(
    blade: Design,
    spec: CircuitSpec,
    parts: Parts,
    choice: PartChoice,
    include_recombinases: bool,
) -> bool:
    """True when the BLADE template beats the product-term layout on length."""
    product = build_design(
        spec, parts, choice, include_recombinases, force_strategy="__product__"
    )
    return blade.total_bp < product.total_bp


def build_design(
    spec: CircuitSpec,
    parts: Parts,
    choice: PartChoice,
    include_recombinases: bool = False,
    force_strategy: str | None = None,
) -> Design:
    """Lay out a complete design for the given part choices."""
    spec.validate()

    rec_names = spec.recombinase_names or parts.orthogonal_set(spec.n)
    if len(rec_names) != spec.n:
        raise SynthesisError(f"{len(rec_names)} recombinases given for {spec.n} inputs")
    recs = [parts.recombinases[r] for r in rec_names]

    promoter = parts.promoters[choice.promoter]
    blocker = parts.terminators[choice.blocker]
    terminator = parts.terminators[choice.terminator]
    linker = parts.linkers[choice.linker]

    covers = {f.name: f.minimise() for f in spec.functions}

    groups = _group_outputs(spec, covers)

    # The BLADE template is a whole-circuit layout, not a per-group one, so it only
    # applies to a single output group. Cheap to price, so try it and keep it if smaller.
    if force_strategy == "__product__":
        force_strategy = None
    elif force_strategy == "blade" or (force_strategy is None and len(groups) == 1):
        blade = _try_blade(
            spec, groups[0], promoter, terminator, recs, linker,
            include_recombinases, parts, choice,
        )
        if blade is not None and (force_strategy == "blade" or _blade_is_smaller(
            blade, spec, parts, choice, include_recombinases
        )):
            return blade
        if force_strategy == "blade":
            raise SynthesisError(
                "blade requires a single output group with at least one ON state"
            )

    cassettes = []
    shared_tail = []
    notes: list[str] = []
    strategies_used = set()

    for names, cover in groups:
        payloads = [spec.payloads[nm] for nm in names]
        group_linker = linker if len(payloads) > 1 else None
        if len(names) > 1:
            notes.append(
                f"Outputs {', '.join(names)} fire in identical states, so they share "
                f"one cassette joined by {linker.name}. That saves "
                f"{len(names) - 1} promoter/polyA pair(s)."
            )

        if not cover.terms:
            notes.append(
                f"Output(s) {', '.join(names)} are never ON. No cassette is emitted; "
                f"the design simply omits that payload."
            )
            continue

        strategy = force_strategy or (
            "tandem_promoters"
            if can_share_payload(cover) and cover.n_terms > 1
            else "independent"
        )
        tandem = strategy == "tandem_promoters"
        if tandem and not can_share_payload(cover):
            raise SynthesisError(
                "tandem_promoters requires every product term to share the same "
                "negative-literal set; this cover does not"
            )
        strategies_used.add(strategy)

        for term in cover.terms:
            cassettes.append(
                build_cassette(
                    implicant=term,
                    outputs=names,
                    payloads=payloads,
                    promoter=promoter,
                    blocker=blocker,
                    terminator=terminator,
                    recombinases=recs,
                    linker=group_linker,
                    carries_payload=not tandem,
                )
            )

        if tandem:
            shared_tail.extend(
                build_shared_payload_tail(
                    neg_lits=cover.terms[0].negative_literals(),
                    outputs=names,
                    payloads=payloads,
                    terminator=terminator,
                    recombinases=recs,
                    linker=group_linker,
                )
            )
            notes.append(
                f"Output(s) {', '.join(names)}: {cover.n_terms} gated promoters in "
                f"tandem share one payload copy, saving "
                f"{(cover.n_terms - 1) * (sum(p.bp for p in payloads) + terminator.bp)} bp "
                f"against one unit per term."
            )

    rec_elements = []
    if include_recombinases:
        rec_elements = build_recombinase_elements(recs, promoter, terminator, linker)
        notes.append(
            "Recombinase coding sequences are included in the budget. Published devices "
            "almost always supply them separately, and carrying them is what usually "
            "makes a design exceed AAV."
        )

    strategy = (
        "mixed"
        if len(strategies_used) > 1
        else (strategies_used.pop() if strategies_used else "empty")
    )

    return Design(
        inputs=list(spec.inputs),
        recombinases=recs,
        cassettes=cassettes,
        strategy=strategy,
        shared_payload_elements=shared_tail,
        recombinase_elements=rec_elements,
        include_recombinases=include_recombinases,
        notes=notes,
    )


@dataclass
class Candidate:
    design: Design
    choice: PartChoice
    circuit_bp: int
    dynamic_range: float
    fits: bool
    headroom: int | None

    @property
    def score(self) -> tuple:
        """Sort key: fitting first, then larger dynamic range, then smaller size."""
        return (0 if self.fits else 1, -self.dynamic_range, self.circuit_bp)


# Smallest first. Smaller polyA signals leak more, trading size for dynamic range.
_TERMINATOR_SEARCH = ["SPA", "SV40pA", "bGHpA", "hGHpA", "tandem"]


def search(
    spec: CircuitSpec,
    parts: Parts,
    vector: Vector,
    include_recombinases: bool = False,
    min_dynamic_range: float = 10.0,
    n_samples: int = 400,
    promoters: list[str] | None = None,
) -> list[Candidate]:
    """Enumerate part choices and score each resulting design."""
    spec.validate()
    prom_opts = promoters or [p for p in PROMOTER_DOWNSIZE_ORDER if p in parts.promoters]
    cands: list[Candidate] = []

    for prom, blk, term in itertools.product(
        prom_opts, _TERMINATOR_SEARCH, _TERMINATOR_SEARCH
    ):
        if blk not in parts.terminators or term not in parts.terminators:
            continue
        choice = PartChoice(promoter=prom, blocker=blk, terminator=term)
        try:
            design = build_design(spec, parts, choice, include_recombinases)
        except SynthesisError:  # pragma: no cover
            continue
        bp = design.total_bp
        verdict = vector.evaluate(bp)

        primary = spec.functions[0].name
        rep = leakage.predict(design, primary, n_samples=n_samples, seed=1)
        dr = rep.dynamic_range_median

        cands.append(
            Candidate(
                design=design,
                choice=choice,
                circuit_bp=bp,
                dynamic_range=1e9 if dr == float("inf") else dr,
                fits=verdict.fits,
                headroom=verdict.headroom_bp,
            )
        )

    cands.sort(key=lambda c: c.score)
    return cands


# Relative gap below which two dynamic ranges count as equal. The leakage model's
# parameters are order-of-magnitude estimates, so finer differences are not real.
DYNAMIC_RANGE_TOLERANCE = 0.10


def _dr_equivalent(a: float, b: float, tol: float = DYNAMIC_RANGE_TOLERANCE) -> bool:
    if a == b:
        return True
    hi = max(a, b)
    return hi > 0 and abs(a - b) / hi <= tol


def best_fitting(
    cands: list[Candidate], min_dynamic_range: float = 10.0
) -> Candidate | None:
    """The smallest design that fits and clears the dynamic-range floor."""
    ok = [c for c in cands if c.fits and c.dynamic_range >= min_dynamic_range]
    if not ok:
        return None
    return min(ok, key=lambda c: (c.circuit_bp, -c.dynamic_range))


def best_effort(cands: list[Candidate]) -> Candidate:
    """Pick when nothing clears the floor: maximise range, break near-ties by size."""
    fitting = [c for c in cands if c.fits] or cands
    best_dr = max(c.dynamic_range for c in fitting)
    tied = [c for c in fitting if _dr_equivalent(c.dynamic_range, best_dr)]
    return min(tied, key=lambda c: c.circuit_bp)


def pareto_front(cands: list[Candidate]) -> list[Candidate]:
    """Designs not dominated on both size and dynamic range."""
    front: list[Candidate] = []
    for c in sorted(cands, key=lambda x: (x.circuit_bp, -x.dynamic_range)):
        dominated = any(
            o.circuit_bp <= c.circuit_bp
            and (
                o.dynamic_range >= c.dynamic_range
                or _dr_equivalent(o.dynamic_range, c.dynamic_range)
            )
            for o in front
        )
        if not dominated:
            front.append(c)
    return front
