"""Packaging-capacity and titre models for delivery vectors."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .parts import AAV_ITR_BP


@dataclass
class DeliveryVerdict:
    vector: str
    fits: bool
    circuit_bp: int
    budget_bp: int | None
    headroom_bp: int | None
    n_vectors: int
    relative_titre: float | None = None
    severity: str = "ok"  # ok | tight | over | infeasible
    notes: list[str] = field(default_factory=list)

    @property
    def utilisation(self) -> float | None:
        if not self.budget_bp:
            return None
        return self.circuit_bp / self.budget_bp

    def summary(self) -> str:
        if self.budget_bp is None:
            base = f"{self.vector}: {self.circuit_bp} bp, no hard packaging limit"
        else:
            pct = 100 * (self.utilisation or 0)
            base = (
                f"{self.vector}: {self.circuit_bp} bp of {self.budget_bp} bp "
                f"({pct:.0f}% used, {self.headroom_bp:+d} bp headroom)"
            )
        if self.relative_titre is not None:
            base += f", predicted titre {self.relative_titre:.2f}x of a minimal vector"
        return base


@dataclass
class Vector:
    name: str
    capacity_bp: int | None
    overhead_bp: int
    n_vectors: int = 1
    #: Fraction of capacity above which packaging degrades.
    tight_fraction: float = 0.92
    note: str = ""
    #: Extra sequence each vector needs when the circuit is split across vectors.
    split_overhead_bp: int = 0
    can_carry_dna_circuit: bool = True

    def budget(self) -> int | None:
        if self.capacity_bp is None:
            return None
        per = self.capacity_bp - self.overhead_bp - self.split_overhead_bp
        return per * self.n_vectors

    def evaluate(self, circuit_bp: int) -> DeliveryVerdict:
        notes = [self.note] if self.note else []

        if not self.can_carry_dna_circuit:
            return DeliveryVerdict(
                vector=self.name,
                fits=False,
                circuit_bp=circuit_bp,
                budget_bp=None,
                headroom_bp=None,
                n_vectors=self.n_vectors,
                severity="infeasible",
                notes=notes
                + [
                    "This modality cannot carry a DNA excision circuit at all. Excision "
                    "logic stores its state in DNA, and mRNA is transient and never "
                    "enters the nucleus as a recombination substrate. Use it to supply "
                    "the recombinase inputs, with the circuit delivered as DNA.",
                ],
            )

        budget = self.budget()
        if budget is None:
            return DeliveryVerdict(
                vector=self.name,
                fits=True,
                circuit_bp=circuit_bp,
                budget_bp=None,
                headroom_bp=None,
                n_vectors=self.n_vectors,
                severity="ok",
                notes=notes,
            )

        headroom = budget - circuit_bp
        fits = headroom >= 0
        if not fits:
            severity = "over"
            notes.append(
                f"Over budget by {-headroom} bp. Oversized AAV genomes are packaged "
                "with 5'-truncated, heterogeneous inserts rather than simply at lower "
                "titre, so this is a correctness problem and not a yield problem."
                if self.name.startswith("AAV")
                else f"Over budget by {-headroom} bp."
            )
        elif circuit_bp > self.tight_fraction * budget:
            severity = "tight"
            notes.append(
                f"Within {100 * (1 - self.tight_fraction):.0f}% of the limit. Real "
                "constructs pick up spacers and scars during cloning; leave margin."
            )
        else:
            severity = "ok"

        if self.n_vectors > 1 and fits:
            notes.append(
                f"Requires co-delivery of {self.n_vectors} vectors to the same cell. "
                "Reconstitution efficiency is the product of per-vector transduction, "
                "so effective efficiency drops roughly as the square at matched dose."
            )

        return DeliveryVerdict(
            vector=self.name,
            fits=fits,
            circuit_bp=circuit_bp,
            budget_bp=budget,
            headroom_bp=headroom,
            n_vectors=self.n_vectors,
            severity=severity,
            notes=notes,
        )


#: LTRs, psi, RRE and cPPT/FLAP. WPRE is optional and not counted.
LENTI_BACKBONE_BP = 2300

LENTI_REFERENCE_INSERT_BP = 4400

# Anchored on the roughly two orders of magnitude lost between a 4.4 kb and an 11.1 kb
# insert. An estimate for ranking vectors, not a calibrated prediction.
LENTI_LOG10_PER_KB = 2.0 / ((11100 - 4400) / 1000)


def lentiviral_relative_titre(insert_bp: int) -> float:
    """Titre relative to a 4.4 kb insert; clamped at 1.0 below that size."""
    excess_kb = max(0.0, (insert_bp - LENTI_REFERENCE_INSERT_BP) / 1000.0)
    return float(10 ** (-LENTI_LOG10_PER_KB * excess_kb))


@dataclass
class LentivirusVector(Vector):
    def evaluate(self, circuit_bp: int) -> DeliveryVerdict:
        v = super().evaluate(circuit_bp)
        v.relative_titre = lentiviral_relative_titre(circuit_bp)
        v.notes.append(
            f"Packaged genome is {circuit_bp + LENTI_BACKBONE_BP} bp including "
            f"{LENTI_BACKBONE_BP} bp of required vector elements."
        )
        if v.relative_titre < 0.1:
            v.severity = "tight"
            v.notes.append(
                f"Predicted titre is {v.relative_titre:.3f}x a minimal vector, i.e. "
                f"roughly {-math.log10(v.relative_titre):.1f} orders of magnitude down. "
                "Expect to need concentration, and expect multiplicity-of-infection "
                "and copy-number control to get harder."
            )
        return v


def vector_catalogue() -> dict[str, Vector]:
    aav_capacity = 4700
    return {
        "aav_single": Vector(
            name="AAV (single vector)",
            capacity_bp=aav_capacity,
            overhead_bp=2 * AAV_ITR_BP,
            n_vectors=1,
            tight_fraction=0.92,
            note="Wild-type AAV genome is 4.7 kb; usable cargo between the ITRs is "
            "about 4.4 kb.",
        ),
        "aav_dual": Vector(
            name="AAV (dual vector, trans-splicing)",
            capacity_bp=aav_capacity,
            overhead_bp=2 * AAV_ITR_BP,
            n_vectors=2,
            tight_fraction=0.88,
            split_overhead_bp=250,
            note="Split across two AAVs reconstituting in the nucleus. Costs a splice "
            "donor/acceptor or split-intein pair per vector and needs both to reach the "
            "same cell.",
        ),
        "aav_triple": Vector(
            name="AAV (triple vector)",
            capacity_bp=aav_capacity,
            overhead_bp=2 * AAV_ITR_BP,
            n_vectors=3,
            tight_fraction=0.85,
            split_overhead_bp=250,
            note="Rarely practical. Three-way co-transduction efficiency is poor; "
            "include it only to show what the design would need.",
        ),
        "lentivirus": LentivirusVector(
            name="Lentivirus",
            capacity_bp=None,
            overhead_bp=LENTI_BACKBONE_BP,
            n_vectors=1,
            note="No hard packaging ceiling, but titre falls steeply with genome size.",
        ),
        "plasmid": Vector(
            name="Plasmid (transient transfection)",
            capacity_bp=None,
            overhead_bp=0,
            n_vectors=1,
            note="What essentially all published mammalian circuit characterisation "
            "uses, HEK293 included. Convenient, and the least predictive of how the "
            "circuit behaves once integrated at single copy.",
        ),
        "landing_pad": Vector(
            name="Genomic landing pad (integrase docking)",
            capacity_bp=None,
            overhead_bp=0,
            n_vectors=1,
            note="Single-copy, defined locus, stable. The most predictable context and "
            "the one where published part parameters are most likely to transfer.",
        ),
        "mrna_lnp": Vector(
            name="mRNA / LNP",
            capacity_bp=None,
            overhead_bp=0,
            n_vectors=1,
            can_carry_dna_circuit=False,
            note="Transient protein delivery.",
        ),
    }


DEFAULT_VECTORS = ["aav_single", "aav_dual", "lentivirus", "landing_pad", "plasmid"]


def evaluate_all(
    circuit_bp: int, vectors: list[str] | None = None
) -> list[DeliveryVerdict]:
    cat = vector_catalogue()
    out = []
    for n in vectors or DEFAULT_VECTORS:
        if n not in cat:
            raise KeyError(f"unknown vector {n!r}; choose from {sorted(cat)}")
        out.append(cat[n].evaluate(circuit_bp))
    return out


def minimum_vectors_needed(circuit_bp: int) -> int:
    """How many AAVs this design would need, ignoring whether that is practical."""
    single = 4700 - 2 * AAV_ITR_BP
    if circuit_bp <= single:
        return 1
    return max(2, math.ceil(circuit_bp / (single - 250)))
