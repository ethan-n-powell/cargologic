"""Leakiness, dynamic range, and -- most usefully -- what is actually limiting them.

READ THIS BEFORE QUOTING ANY NUMBER FROM THIS MODULE.

The size verdict in `delivery.py` is arithmetic. This module is a model, and models of
mammalian circuit behaviour have a poor recent track record: a 2026 mammalian circuit CAD
paper found its ODE predictions "did not align well" with its own measurements, and a
2026 ultra-high-throughput study found a plain neural network beat a biophysically
consistent mechanistic model by a wide margin on held-out designs (r^2 ~0.88 vs ~0.37).

So this module is deliberately modest in what it claims and deliberately specific in
what it diagnoses. It propagates three parameters we have order-of-magnitude estimates
for -- excision efficiency, spontaneous excision leak, and terminator read-through --
through the circuit's state graph, reports a distribution rather than a point, and then
runs a sensitivity analysis to say WHICH of the three is limiting you.

That last part is the bit worth having. "Your AND-NOT gate is excision-limited, so a
tighter polyA will not help and a better recombinase will" is actionable and robust even
when the absolute numbers are off, because it depends on the ordering of effects rather
than their magnitudes.

Not modelled: transcriptional squelching, resource competition between transgenes,
retroactivity, position effects, silencing, cell-type differences. Each has broken a
published mammalian design.
"""

from __future__ import annotations

import random
import statistics
from dataclasses import dataclass, field

from .architecture import Cassette, Design


@dataclass
class Params:
    """One draw of the model's parameters."""

    #: Probability each input's recombinase has acted, given that input's state.
    fired: dict[int, float]
    blocker_readthrough: float
    terminator_readthrough: float


@dataclass
class StateExpression:
    minterm: int
    bits: str
    intended: bool
    level_median: float
    level_lo: float
    level_hi: float

    @property
    def intended_str(self) -> str:
        return "ON" if self.intended else "OFF"


@dataclass
class Mechanism:
    """How much the worst OFF state improves if one mechanism were made perfect."""

    name: str
    off_level: float
    improvement: float
    advice: str


@dataclass
class LeakageReport:
    output: str
    states: list[StateExpression]
    on_level_median: float
    off_level_median: float
    dynamic_range_median: float
    dynamic_range_lo: float
    dynamic_range_hi: float
    worst_off_state: str
    weakest_on_state: str
    n_samples: int
    limiting: Mechanism | None = None
    mechanisms: list[Mechanism] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def summary(self) -> str:
        return (
            f"{self.output}: predicted ON/OFF ratio "
            f"{self.dynamic_range_median:.1f}x "
            f"(80% interval {self.dynamic_range_lo:.1f}-{self.dynamic_range_hi:.1f}x)"
        )


def _sample(rng: random.Random, lo: float, hi: float) -> float:
    """Draw from a parameter's reported range.

    Triangular with the mode at the midpoint: we have a literature range and no basis
    for a sharper shape, and a uniform draw would overweight the extremes.
    """
    if hi <= lo:
        return lo
    return rng.triangular(lo, hi, (lo + hi) / 2)


def _draw_params(rng: random.Random, design: Design, minterm: int, cas: Cassette) -> Params:
    n = len(design.inputs)
    fired: dict[int, float] = {}
    for i in range(n):
        rec = design.recombinase_for(i)
        on = bool((minterm >> (n - 1 - i)) & 1)
        fired[i] = _sample(rng, *(rec.excision_efficiency if on else rec.leak))
    return Params(
        fired=fired,
        blocker_readthrough=_sample(rng, *cas.blocker.readthrough),
        terminator_readthrough=_sample(rng, *cas.terminator.readthrough),
    )


def _midpoint_params(design: Design, minterm: int, cas: Cassette) -> Params:
    n = len(design.inputs)
    fired: dict[int, float] = {}
    for i in range(n):
        rec = design.recombinase_for(i)
        on = bool((minterm >> (n - 1 - i)) & 1)
        fired[i] = rec.eff_mid if on else rec.leak_mid
    return Params(
        fired=fired,
        blocker_readthrough=sum(cas.blocker.readthrough) / 2,
        terminator_readthrough=sum(cas.terminator.readthrough) / 2,
    )


def _cassette_reach(cas: Cassette, p: Params) -> float:
    """Probability transcription from this promoter reaches the payload position.

    Each positive literal is a blocker that must be excised. An unexcised blocker still
    passes `blocker_readthrough` of polymerases, which is the floor you cannot design
    away without a tighter terminator.
    """
    out = 1.0
    for i in cas.positive_literals:
        f = p.fired[i]
        out *= f + (1.0 - f) * p.blocker_readthrough
    return out


def _payload_present(neg_literals: list[int], p: Params) -> float:
    """Probability the payload has NOT been excised by a negative-literal recombinase.

    This is where NOT gates live or die. If excision is only 85% efficient, 15% of cells
    keep the payload and the OFF state is 15% ON -- regardless of how good your
    terminators are.
    """
    out = 1.0
    for j in neg_literals:
        out *= 1.0 - p.fired[j]
    return out


def _level(design: Design, relevant: list[Cassette], minterm: int, p: Params) -> float:
    """Predicted output level for one state under one parameter draw."""
    shared = design.strategy == "tandem_promoters"

    if shared:
        neg = relevant[0].negative_literals
        present = _payload_present(neg, p)
        miss = 1.0
        strength = 0.0
        for c in relevant:
            miss *= 1.0 - _cassette_reach(c, p)
            strength = max(strength, c.promoter.strength)
        level = present * (1.0 - miss) * strength
        floor_strength = strength
    else:
        level = 0.0
        floor_strength = 0.0
        for c in relevant:
            reach = _cassette_reach(c, p)
            present = _payload_present(c.negative_literals, p)
            level += reach * present * c.promoter.strength
            floor_strength = max(floor_strength, c.promoter.strength)

    # Read-through from an upstream promoter into a payload that is still present.
    # Scaled by promoter strength, because that is where the polymerases come from --
    # otherwise a strong promoter would appear to improve the ON/OFF ratio for free.
    any_present = max(
        (_payload_present(c.negative_literals, p) for c in relevant), default=1.0
    )
    level += p.terminator_readthrough * 0.02 * floor_strength * any_present
    return level


def _sensitivity(
    design: Design, relevant: list[Cassette], worst_off_minterm: int
) -> list[Mechanism]:
    """Which mechanism, made perfect, most improves the worst OFF state?

    Robust to the absolute parameter values being wrong, because it compares effects
    against each other rather than against an absolute scale.
    """
    cas = relevant[0]
    base_p = _midpoint_params(design, worst_off_minterm, cas)
    base = _level(design, relevant, worst_off_minterm, base_p)
    if base <= 0:
        return []

    n = len(design.inputs)
    out: list[Mechanism] = []

    # (a) Perfect excision: every recombinase whose input is ON fires in every cell.
    p = Params(dict(base_p.fired), base_p.blocker_readthrough, base_p.terminator_readthrough)
    for i in range(n):
        if bool((worst_off_minterm >> (n - 1 - i)) & 1):
            p.fired[i] = 1.0
    lvl = _level(design, relevant, worst_off_minterm, p)
    out.append(
        Mechanism(
            "incomplete excision",
            lvl,
            (base - lvl) / base,
            "Use a more efficient recombinase, raise its dose or expression time, or "
            "switch to an irreversible serine integrase. A tighter polyA will not help.",
        )
    )

    # (b) No spontaneous excision in the absence of input.
    p = Params(dict(base_p.fired), base_p.blocker_readthrough, base_p.terminator_readthrough)
    for i in range(n):
        if not bool((worst_off_minterm >> (n - 1 - i)) & 1):
            p.fired[i] = 0.0
    lvl = _level(design, relevant, worst_off_minterm, p)
    out.append(
        Mechanism(
            "spontaneous recombinase activity",
            lvl,
            (base - lvl) / base,
            "Tighten input control: inducible or split recombinases, lower basal "
            "expression, degron-tagged recombinase.",
        )
    )

    # (c) Perfect termination.
    p = Params(dict(base_p.fired), 0.0, 0.0)
    lvl = _level(design, relevant, worst_off_minterm, p)
    out.append(
        Mechanism(
            "terminator read-through",
            lvl,
            (base - lvl) / base,
            "Use a tighter or tandem polyA as the blocker. Costs bp -- see the Pareto "
            "front.",
        )
    )

    out.sort(key=lambda m: -m.improvement)
    return out


def predict(
    design: Design,
    output: str,
    n_samples: int = 2000,
    seed: int = 0,
) -> LeakageReport:
    """Monte-Carlo the circuit across all input states, then diagnose the limiter."""
    rng = random.Random(seed)
    n = len(design.inputs)

    relevant = [c for c in design.cassettes if output in c.outputs]
    if not relevant:
        states = [
            StateExpression(m, format(m, f"0{n}b"), False, 0.0, 0.0, 0.0)
            for m in range(1 << n)
        ]
        return LeakageReport(
            output=output,
            states=states,
            on_level_median=0.0,
            off_level_median=0.0,
            dynamic_range_median=1.0,
            dynamic_range_lo=1.0,
            dynamic_range_hi=1.0,
            worst_off_state="-",
            weakest_on_state="-",
            n_samples=0,
            warnings=["This output is never expressed in any input state."],
        )

    intended = {m: any(c.implicant.covers(m) for c in relevant) for m in range(1 << n)}
    per_state: dict[int, list[float]] = {m: [] for m in range(1 << n)}
    dyn_samples: list[float] = []

    for _ in range(n_samples):
        levels: dict[int, float] = {}
        for m in range(1 << n):
            p = _draw_params(rng, design, m, relevant[0])
            lvl = _level(design, relevant, m, p)
            levels[m] = lvl
            per_state[m].append(lvl)

        ons = [levels[m] for m in range(1 << n) if intended[m]]
        offs = [levels[m] for m in range(1 << n) if not intended[m]]
        if ons and offs:
            hi_off = max(offs)
            dyn_samples.append(min(ons) / hi_off if hi_off > 0 else float("inf"))
        elif ons:
            dyn_samples.append(float("inf"))

    def pct(xs: list[float], q: float) -> float:
        if not xs:
            return 0.0
        s = sorted(xs)
        k = min(len(s) - 1, max(0, int(q * (len(s) - 1))))
        return s[k]

    states = [
        StateExpression(
            minterm=m,
            bits=format(m, f"0{n}b"),
            intended=intended[m],
            level_median=statistics.median(per_state[m]),
            level_lo=pct(per_state[m], 0.10),
            level_hi=pct(per_state[m], 0.90),
        )
        for m in range(1 << n)
    ]

    on_states = [s for s in states if s.intended]
    off_states = [s for s in states if not s.intended]
    on_med = min((s.level_median for s in on_states), default=0.0)
    off_med = max((s.level_median for s in off_states), default=0.0)

    finite = [d for d in dyn_samples if d != float("inf")]
    dyn_med = statistics.median(finite) if finite else float("inf")
    dyn_lo = pct(finite, 0.10) if finite else float("inf")
    dyn_hi = pct(finite, 0.90) if finite else float("inf")

    worst = max(off_states, key=lambda s: s.level_median) if off_states else None
    mechanisms = _sensitivity(design, relevant, worst.minterm) if worst else []
    limiting = mechanisms[0] if mechanisms and mechanisms[0].improvement > 0.05 else None

    warnings: list[str] = []
    if finite and dyn_med < 10:
        warnings.append(
            f"Predicted dynamic range is only ~{dyn_med:.1f}x. Below roughly 10x, "
            f"flow-cytometry populations overlap and the gate will not read as digital."
        )
    if worst and on_med > 0 and worst.level_median > 0.1 * on_med:
        frac = 100 * worst.level_median / on_med
        if limiting:
            warnings.append(
                f"State {worst.bits} should be OFF but sits at {frac:.0f}% of the "
                f"weakest ON state. Dominant cause: {limiting.name} "
                f"({100 * limiting.improvement:.0f}% of the excess). {limiting.advice}"
            )
        else:
            warnings.append(
                f"State {worst.bits} should be OFF but sits at {frac:.0f}% of the "
                f"weakest ON state."
            )
    max_blockers = max((len(c.positive_literals) for c in relevant), default=0)
    if max_blockers >= 3:
        warnings.append(
            f"One cassette carries {max_blockers} blockers in series. All must be "
            f"excised for the term to fire, so the ON level is a product of three "
            f"efficiencies -- expect a weak ON, not merely a clean OFF."
        )
    n_neg = max((len(c.negative_literals) for c in relevant), default=0)
    if n_neg >= 1 and limiting and limiting.name == "incomplete excision":
        warnings.append(
            "This is the characteristic weakness of NOT logic built on excision: the "
            "OFF state requires the payload to be physically removed, so residual "
            "unexcised cells read as ON. It is a population-heterogeneity problem, not "
            "a transcriptional one, and sorting or an irreversible integrase helps "
            "where better terminators do not."
        )

    return LeakageReport(
        output=output,
        states=states,
        on_level_median=on_med,
        off_level_median=off_med,
        dynamic_range_median=dyn_med,
        dynamic_range_lo=dyn_lo,
        dynamic_range_hi=dyn_hi,
        worst_off_state=worst.bits if worst else "-",
        weakest_on_state=min(on_states, key=lambda s: s.level_median).bits
        if on_states
        else "-",
        n_samples=n_samples,
        limiting=limiting,
        mechanisms=mechanisms,
        warnings=warnings,
    )
