"""Part library for recombinase-based mammalian logic circuits."""

from __future__ import annotations

from dataclasses import dataclass, field, replace


@dataclass(frozen=True)
class RecombinaseSystem:
    name: str
    site_name: str
    site_bp: int
    cds_bp: int
    excision_efficiency: tuple[float, float]
    leak: tuple[float, float]  # excision fraction with no recombinase supplied
    cross_reacts: frozenset[str] = frozenset()
    family: str = "tyrosine"
    note: str = ""

    @property
    def eff_mid(self) -> float:
        return sum(self.excision_efficiency) / 2

    @property
    def leak_mid(self) -> float:
        return sum(self.leak) / 2


@dataclass(frozen=True)
class Promoter:
    name: str
    bp: int
    strength: float  # arbitrary units, relative to EF1a full = 1.0
    note: str = ""


@dataclass(frozen=True)
class Terminator:
    name: str
    bp: int
    readthrough: tuple[float, float]  # fraction of polymerases getting past it
    note: str = ""


@dataclass(frozen=True)
class Payload:
    """Cargo size, coding sequence only, no promoter/polyA."""

    name: str
    bp: int
    note: str = ""


@dataclass(frozen=True)
class Linker:
    """Ribosome-skipping or cleavable linker for co-expressing payloads."""

    name: str
    bp: int
    note: str = ""


# Efficiency and leak values are literature estimates with wide ranges, for ranking
# designs rather than predicting numbers. Override them with measured values.
RECOMBINASES: dict[str, RecombinaseSystem] = {
    "Cre": RecombinaseSystem(
        name="Cre",
        site_name="loxP",
        site_bp=34,
        cds_bp=1032,
        excision_efficiency=(0.85, 0.99),
        leak=(0.001, 0.05),
        family="tyrosine",
        note="Workhorse. Highest efficiency, but also the most reported leak; "
        "Cre toxicity and pseudo-loxP activity in mammalian genomes are real.",
    ),
    "FlpO": RecombinaseSystem(
        name="FlpO",
        site_name="FRT",
        site_bp=34,
        cds_bp=1272,
        excision_efficiency=(0.70, 0.95),
        leak=(0.001, 0.03),
        family="tyrosine",
        note="Mammalian codon-optimised Flp. Somewhat less efficient than Cre. "
        "Orthogonal to Cre -- the classic second input in BLADE-style devices.",
    ),
    "Dre": RecombinaseSystem(
        name="Dre",
        site_name="rox",
        site_bp=32,
        cds_bp=1035,
        excision_efficiency=(0.70, 0.95),
        leak=(0.001, 0.03),
        family="tyrosine",
        note="Orthogonal to Cre and Flp. Good third input.",
    ),
    "VCre": RecombinaseSystem(
        name="VCre",
        site_name="VloxP",
        site_bp=34,
        cds_bp=1035,
        excision_efficiency=(0.60, 0.90),
        leak=(0.001, 0.03),
        family="tyrosine",
        note="Orthogonal to Cre/Flp/Dre. Less widely characterised; treat its "
        "parameters as weakly constrained.",
    ),
    "Bxb1": RecombinaseSystem(
        name="Bxb1",
        site_name="attB/attP",
        site_bp=88,  # attB 38 + attP 50, counted as the pair
        cds_bp=1500,
        excision_efficiency=(0.75, 0.95),
        leak=(0.0005, 0.01),
        family="serine",
        note="Best-performing serine integrase in human cells; irreversible. Site pair "
        "is bulkier than loxP (~88 bp vs 68 bp per pair) -- this costs you in AAV.",
    ),
    "phiC31": RecombinaseSystem(
        name="phiC31",
        site_name="attB/attP",
        site_bp=73,  # attB 34 + attP 39
        cds_bp=1812,
        excision_efficiency=(0.40, 0.80),
        leak=(0.0005, 0.01),
        family="serine",
        note="Irreversible. Lower efficiency than Bxb1 in human cells and a large CDS. "
        "Has documented pseudo-site activity in the human genome.",
    ),
}

# Pairs with published cross-reactivity concerns: attB/attP scaffolds can cross-react
# unless the sites come from a validated orthogonal set.
_KNOWN_CROSSTALK: list[tuple[str, str]] = [
    ("Bxb1", "phiC31"),
]


def crosstalk_risk(a: str, b: str) -> bool:
    """True if this recombinase pair has published orthogonality concerns."""
    return (a, b) in _KNOWN_CROSSTALK or (b, a) in _KNOWN_CROSSTALK


# Preference order when the compiler picks recombinases.
DEFAULT_RECOMBINASE_ORDER = ["Cre", "FlpO", "Dre", "VCre", "Bxb1", "phiC31"]


PROMOTERS: dict[str, Promoter] = {
    "minCMV": Promoter("minCMV", 100, 0.15, "Minimal CMV core. Tiny and weak -- the "
                       "right default when you are fighting a packaging limit."),
    "SV40": Promoter("SV40", 350, 0.35, "Compact, moderate strength."),
    "EFS": Promoter("EFS", 212, 0.55, "Intron-less EF1a core. Excellent size/strength "
                    "trade-off and resists silencing better than CMV."),
    "hPGK": Promoter("hPGK", 500, 0.45, "Modest, steady, low silencing."),
    "CMV": Promoter("CMV", 600, 1.10, "Very strong in HEK293. Prone to silencing in "
                    "primary cells and in vivo -- a common reason circuits that work "
                    "in HEK fail in T cells."),
    "UBC": Promoter("UBC", 400, 0.40, "Human ubiquitin C. Reliable, moderate."),
    "EF1a": Promoter("EF1a", 1200, 1.00, "Full EF1a with intron. Strong and durable, "
                     "but 1.2 kb is a quarter of your AAV budget."),
    "CAG": Promoter("CAG", 1700, 1.40, "Strongest common option. Almost never fits in "
                    "AAV alongside a real payload."),
}

DEFAULT_PROMOTER = "EFS"
# Order tried when downsizing to meet a budget.
PROMOTER_DOWNSIZE_ORDER = ["EF1a", "CMV", "hPGK", "UBC", "SV40", "EFS", "minCMV"]


TERMINATORS: dict[str, Terminator] = {
    "SPA": Terminator("SPA", 49, (0.02, 0.10), "Synthetic polyA. Smallest option; "
                      "leakier, which directly raises your circuit's OFF floor."),
    "SV40pA": Terminator("SV40pA", 135, (0.01, 0.05), "Compact and well behaved."),
    "bGHpA": Terminator("bGHpA", 225, (0.005, 0.03), "Bovine growth hormone polyA. "
                        "The usual default; tighter than SV40."),
    "hGHpA": Terminator("hGHpA", 480, (0.003, 0.02), "Tightest common option, at a "
                        "real size cost."),
    "tandem": Terminator("tandem", 450, (0.001, 0.008), "Two polyA signals in series. "
                         "Use when OFF-state leak is what is killing your dynamic "
                         "range and you can afford the bp."),
}

DEFAULT_TERMINATOR = "bGHpA"
# Terminator used as the gated-OFF blocker.
DEFAULT_BLOCKER = "bGHpA"


PAYLOADS: dict[str, Payload] = {
    "EGFP": Payload("EGFP", 720, "Reporter."),
    "mCherry": Payload("mCherry", 711, "Reporter."),
    "mKate2": Payload("mKate2", 696, "Reporter used in the COMET part library."),
    "luciferase": Payload("luciferase", 1653, "Firefly luc2."),
    "CAR_2g": Payload("CAR_2g", 1500, "Second-generation CAR: scFv-CD8a hinge/TM-CD28-"
                      "CD3zeta. The canonical cell-therapy payload."),
    "CAR_3g": Payload("CAR_3g", 1800, "Third-generation CAR with an extra costim domain."),
    "iCasp9": Payload("iCasp9", 1245, "Inducible caspase-9 safety switch."),
    "HSV_TK": Payload("HSV_TK", 1131, "Ganciclovir-inducible suicide gene."),
    "tEGFR": Payload("tEGFR", 1100, "Truncated EGFR selection/ablation marker."),
    "IL12": Payload("IL12", 1600, "IL-12 p35+p40, linked."),
    "IL2": Payload("IL2", 462, "Interleukin-2."),
    "SpCas9": Payload("SpCas9", 4104, "Streptococcus pyogenes Cas9. Does not fit in a "
                      "single AAV with any promoter -- this is why the field uses "
                      "split-intein dual AAV."),
    "dCas9_VPR": Payload("dCas9_VPR", 6000, "Catalytically dead Cas9 fused to VPR."),
    "synNotch": Payload("synNotch", 2400, "Synthetic Notch receptor."),
}

LINKERS: dict[str, Linker] = {
    "P2A": Linker("P2A", 66, "Highest skipping efficiency of the common 2A peptides."),
    "T2A": Linker("T2A", 54, "Smallest common 2A."),
    "E2A": Linker("E2A", 72, ""),
    "F2A": Linker("F2A", 72, ""),
    "IRES": Linker("IRES", 580, "EMCV IRES. Large, and downstream expression is much "
                   "weaker than upstream -- 2A is usually the better choice."),
}

DEFAULT_LINKER = "P2A"


# AAV inverted terminal repeat, one at each end of the genome.
AAV_ITR_BP = 130

KOZAK_BP = 10

# Spacer/scar left by assembly between adjacent elements.
JUNCTION_BP = 6


@dataclass
class Parts:
    """A resolved, possibly user-customised part library."""

    recombinases: dict[str, RecombinaseSystem] = field(
        default_factory=lambda: dict(RECOMBINASES)
    )
    promoters: dict[str, Promoter] = field(default_factory=lambda: dict(PROMOTERS))
    terminators: dict[str, Terminator] = field(
        default_factory=lambda: dict(TERMINATORS)
    )
    payloads: dict[str, Payload] = field(default_factory=lambda: dict(PAYLOADS))
    linkers: dict[str, Linker] = field(default_factory=lambda: dict(LINKERS))

    def override_recombinase(self, name: str, **kwargs) -> None:
        if name not in self.recombinases:
            raise KeyError(f"unknown recombinase {name!r}")
        self.recombinases[name] = replace(self.recombinases[name], **kwargs)

    def override_promoter(self, name: str, **kwargs) -> None:
        if name not in self.promoters:
            raise KeyError(f"unknown promoter {name!r}")
        self.promoters[name] = replace(self.promoters[name], **kwargs)

    def override_terminator(self, name: str, **kwargs) -> None:
        if name not in self.terminators:
            raise KeyError(f"unknown terminator {name!r}")
        self.terminators[name] = replace(self.terminators[name], **kwargs)

    def add_payload(self, name: str, bp: int, note: str = "") -> Payload:
        self.payloads[name] = Payload(name, int(bp), note)
        return self.payloads[name]

    def resolve_payload(self, spec: str | int | Payload) -> Payload:
        """Accept a library name, a raw bp count, or a Payload."""
        if isinstance(spec, Payload):
            return spec
        if isinstance(spec, int):
            return Payload(f"payload_{spec}bp", spec, "user-specified size")
        if isinstance(spec, str):
            if spec in self.payloads:
                return self.payloads[spec]
            if ":" in spec:  # "name:1500"
                nm, _, num = spec.partition(":")
                return Payload(nm, int(num), "user-specified size")
            if spec.isdigit():
                return Payload(f"payload_{spec}bp", int(spec), "user-specified size")
            raise KeyError(
                f"unknown payload {spec!r}; pass a bp count or one of "
                f"{sorted(self.payloads)}"
            )
        raise TypeError(f"cannot resolve payload from {spec!r}")

    def orthogonal_set(self, n: int, order: list[str] | None = None) -> list[str]:
        """Pick `n` mutually orthogonal recombinases, preferring well-characterised ones."""
        chosen: list[str] = []
        for name in order or DEFAULT_RECOMBINASE_ORDER:
            if name not in self.recombinases:
                continue
            if any(crosstalk_risk(name, c) for c in chosen):
                continue
            chosen.append(name)
            if len(chosen) == n:
                return chosen
        raise ValueError(
            f"cannot find {n} mutually orthogonal recombinases in the library "
            f"(found {len(chosen)}: {chosen}). Circuits with more than ~4 inputs need "
            f"a validated orthogonal integrase set that this library does not ship."
        )
