"""Synthesis: Boolean spec -> DNA design, and the search for one that fits.

Two jobs live here.

1. BUILD. Minimise each output's function, assign orthogonal recombinases to inputs,
   pick a sharing strategy, and lay out cassettes.

2. SEARCH. When the design does not fit the target vector, look for one that does. The
   levers are a smaller promoter, a smaller blocker/terminator, and payload sharing --
   and they are not free: a smaller polyA leaks more, which lowers dynamic range. So the
   search is genuinely multi-objective and `pareto_front` returns the trade-off curve
   rather than pretending there is a single right answer.

The thing worth noticing is that the single biggest size lever is upstream of all of
this: the number of product terms after logic minimisation. Shaving a promoter saves a
few hundred bp; eliminating a product term saves a kilobase and change.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

from . import leakage
from .architecture import (
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
    """What the user asked for."""

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
    """The tunable part selections the search explores."""

    promoter: str = DEFAULT_PROMOTER
    blocker: str = DEFAULT_BLOCKER
    terminator: str = DEFAULT_TERMINATOR
    linker: str = DEFAULT_LINKER

    def key(self) -> tuple[str, str, str]:
        return (self.promoter, self.blocker, self.terminator)


def _group_outputs(
    spec: CircuitSpec, covers: dict[str, Cover]
) -> list[tuple[list[str], Cover]]:
    """Outputs whose minimised covers are identical can share cassettes via a 2A linker.

    Two payloads expressed in exactly the same input states do not need separate
    promoters, blockers and polyA signals -- they need one cassette and a 2A peptide.
    This is a pure win and the compiler takes it automatically.
    """
    groups: dict[tuple, list[str]] = {}
    for f in spec.functions:
        cov = covers[f.name]
        sig = tuple(sorted(t.pattern() for t in cov.terms))
        groups.setdefault(sig, []).append(f.name)
    out = []
    for sig, names in groups.items():
        out.append((names, covers[names[0]]))
    return out


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
        raise SynthesisError(
            f"{len(rec_names)} recombinases given for {spec.n} inputs"
        )
    recs = [parts.recombinases[r] for r in rec_names]

    promoter = parts.promoters[choice.promoter]
    blocker = parts.terminators[choice.blocker]
    terminator = parts.terminators[choice.terminator]
    linker = parts.linkers[choice.linker]

    covers = {f.name: f.minimise() for f in spec.functions}
    groups = _group_outputs(spec, covers)

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

        shareable = can_share_payload(cover) and cover.n_terms > 1
        strategy = force_strategy or ("tandem_promoters" if shareable else "independent")
        if strategy == "tandem_promoters" and not can_share_payload(cover):
            raise SynthesisError(
                "tandem_promoters requires every product term to share the same "
                "negative-literal set; this cover does not"
            )
        strategies_used.add(strategy)

        if strategy == "tandem_promoters":
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
                        carries_payload=False,
                    )
                )
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
        else:
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
                        carries_payload=True,
                    )
                )

    rec_elements = []
    if include_recombinases:
        rec_elements = build_recombinase_elements(
            recs, promoter, terminator, parts.linkers[choice.linker]
        )
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


# --------------------------------------------------------------------------------------
# Search
# --------------------------------------------------------------------------------------


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
        """Prefer fitting, then larger dynamic range, then smaller."""
        return (0 if self.fits else 1, -self.dynamic_range, self.circuit_bp)


#: Blocker/terminator options the search considers, smallest first. Smaller polyA
#: signals leak more, so these trade size against dynamic range.
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
    """Enumerate part choices and score each resulting design.

    Small enough to brute force: a handful of promoters times a handful of terminators.
    No need for anything cleverer, and brute force means the Pareto front is exact.
    """
    spec.validate()
    prom_opts = promoters or [
        p for p in PROMOTER_DOWNSIZE_ORDER if p in parts.promoters
    ]
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

        # Dynamic range is the binding constraint for the primary output.
        primary = spec.functions[0].name
        rep = leakage.predict(design, primary, n_samples=n_samples, seed=1)
        dr = rep.dynamic_range_median
        dr = 1e9 if dr == float("inf") else dr

        cands.append(
            Candidate(
                design=design,
                choice=choice,
                circuit_bp=bp,
                dynamic_range=dr,
                fits=verdict.fits,
                headroom=verdict.headroom_bp,
            )
        )

    cands.sort(key=lambda c: c.score)
    return cands


#: Two dynamic ranges within this relative distance count as equivalent. The underlying
#: parameters are order-of-magnitude estimates, so distinguishing 6.0x from 6.1x is
#: false precision -- and acting on it makes the compiler pick needlessly large designs.
DYNAMIC_RANGE_TOLERANCE = 0.10


def _dr_equivalent(a: float, b: float, tol: float = DYNAMIC_RANGE_TOLERANCE) -> bool:
    if a == b:
        return True
    hi = max(a, b)
    return hi > 0 and abs(a - b) / hi <= tol


def best_fitting(
    cands: list[Candidate], min_dynamic_range: float = 10.0
) -> Candidate | None:
    """The smallest design that fits AND clears the dynamic-range floor.

    Among designs whose dynamic ranges are equivalent to within the model's resolution,
    always take the smallest: spending 1300 bp of AAV budget to buy a difference the
    model cannot resolve is a bad trade.
    """
    ok = [c for c in cands if c.fits and c.dynamic_range >= min_dynamic_range]
    if not ok:
        return None
    best_dr = max(c.dynamic_range for c in ok)
    tied = [c for c in ok if _dr_equivalent(c.dynamic_range, best_dr)]
    # Prefer the smallest among the tied set, but do not inflate size to chase DR:
    # anything clearing the floor is acceptable, so the smallest such design wins.
    return min(ok, key=lambda c: (c.circuit_bp, -c.dynamic_range))


def best_effort(cands: list[Candidate]) -> Candidate:
    """Pick when nothing clears the dynamic-range floor.

    Maximise dynamic range, then break ties-within-tolerance by size. Without the
    tolerance this reliably picks the largest design for a vanishing gain.
    """
    fitting = [c for c in cands if c.fits] or cands
    best_dr = max(c.dynamic_range for c in fitting)
    tied = [c for c in fitting if _dr_equivalent(c.dynamic_range, best_dr)]
    return min(tied, key=lambda c: c.circuit_bp)


def pareto_front(cands: list[Candidate]) -> list[Candidate]:
    """Designs not dominated on both size and dynamic range.

    This is the honest output when a circuit is near its packaging limit: there is no
    single best design, there is a curve, and where to sit on it is the user's call.

    Candidates whose dynamic range is equivalent within the model's resolution are
    collapsed to the smallest one, so a front of twelve identical-performing designs
    does not masquerade as a trade-off.
    """
    front: list[Candidate] = []
    for c in sorted(cands, key=lambda x: (x.circuit_bp, -x.dynamic_range)):
        dominated = any(
            o.circuit_bp <= c.circuit_bp
            and (o.dynamic_range >= c.dynamic_range or _dr_equivalent(o.dynamic_range, c.dynamic_range))
            for o in front
        )
        if not dominated:
            front.append(c)
    return front
