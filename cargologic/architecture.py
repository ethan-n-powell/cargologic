"""DNA architecture: turning a logic cover into a physical layout with a length.

HOW THE LOGIC IS REALISED IN DNA
--------------------------------
A single-layer excision circuit implements one product term per cassette:

    [promoter] <a> [blocker polyA] <a>  <b| [payload] [polyA] |b>

  * Each POSITIVE literal (input a must be TRUE) becomes a "blocker": a polyA signal
    flanked by recombinase-a sites, sitting between the promoter and the payload.
    Recombinase a excises the blocker, unblocking transcription. Several blockers in
    series give you AND over positive literals -- every one must be excised.

  * Each NEGATIVE literal (input b must be FALSE) flanks the payload with
    recombinase-b sites. If b fires, the payload is excised and the promoter runs into
    the trailing polyA. Nesting several of these gives NOT b1 AND NOT b2 ...

  * OR across product terms needs one cassette per term. This is the expensive part and
    it is why `logic.minimise` matters: each extra term is another promoter, payload and
    polyA, typically well over a kilobase.

STRATEGIES
----------
`independent`       One full transcription unit per product term. Always valid.
`tandem_promoters`  Valid only when every term shares the same negative-literal set.
                    Gated promoters sit in tandem upstream of ONE shared payload, saving
                    (k-1) copies of payload+polyA. Big win when it applies.

WHAT THIS MODULE DOES AND DOES NOT CLAIM
----------------------------------------
Lengths are arithmetic over real part sizes and are trustworthy. The claim that a given
layout *functions* rests on published single-layer recombinase devices behaving as
advertised, which is well supported for 2-3 inputs and progressively less so beyond.
This module emits an annotated feature table, not nucleotides: you take the actual part
sequences from your own plasmids or Addgene. Emitting invented sequence would look more
impressive and be less useful.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .logic import Cover, Implicant
from .parts import (
    JUNCTION_BP,
    KOZAK_BP,
    Linker,
    Parts,
    Payload,
    Promoter,
    RecombinaseSystem,
    Terminator,
)


@dataclass
class Element:
    """One annotated span of the design."""

    kind: str  # promoter | site | terminator | payload | linker | kozak | itr | spacer
    name: str
    bp: int
    note: str = ""
    # Which input index this element serves, where meaningful.
    input_index: int | None = None

    def __str__(self) -> str:
        return f"{self.name}({self.bp})"


@dataclass
class Cassette:
    """One product term, realised as DNA."""

    outputs: list[str]
    payloads: list[Payload]
    promoter: Promoter
    implicant: Implicant
    positive_literals: list[int]
    negative_literals: list[int]
    blocker: Terminator
    terminator: Terminator
    elements: list[Element] = field(default_factory=list)
    carries_payload: bool = True
    linker: Linker | None = None

    @property
    def bp(self) -> int:
        return sum(e.bp for e in self.elements)

    @property
    def label(self) -> str:
        return "+".join(self.outputs)


@dataclass
class Design:
    """A complete circuit layout."""

    inputs: list[str]
    recombinases: list[RecombinaseSystem]  # index-aligned with `inputs`
    cassettes: list[Cassette]
    strategy: str
    shared_payload_elements: list[Element] = field(default_factory=list)
    recombinase_elements: list[Element] = field(default_factory=list)
    include_recombinases: bool = False
    notes: list[str] = field(default_factory=list)

    @property
    def elements(self) -> list[Element]:
        out: list[Element] = []
        for c in self.cassettes:
            out.extend(c.elements)
        out.extend(self.shared_payload_elements)
        out.extend(self.recombinase_elements)
        return out

    @property
    def circuit_bp(self) -> int:
        """Length of the logic circuit itself, excluding recombinase expression."""
        return sum(c.bp for c in self.cassettes) + sum(
            e.bp for e in self.shared_payload_elements
        )

    @property
    def recombinase_bp(self) -> int:
        return sum(e.bp for e in self.recombinase_elements)

    @property
    def total_bp(self) -> int:
        return self.circuit_bp + self.recombinase_bp

    @property
    def n_terms(self) -> int:
        return len(self.cassettes)

    def recombinase_for(self, input_index: int) -> RecombinaseSystem:
        return self.recombinases[input_index]

    # ---------------------------------------------------------------- reporting

    def feature_table(self) -> list[dict]:
        """Coordinates for every element, 1-based inclusive, ready for GenBank."""
        rows = []
        pos = 1
        for e in self.elements:
            rows.append(
                {
                    "start": pos,
                    "end": pos + e.bp - 1,
                    "bp": e.bp,
                    "kind": e.kind,
                    "name": e.name,
                    "input": self.inputs[e.input_index]
                    if e.input_index is not None
                    else "",
                    "note": e.note,
                }
            )
            pos += e.bp
        return rows

    def to_genbank(self, name: str = "cargologic_design") -> str:
        """A minimal GenBank-format feature table with N-masked sequence.

        The features and coordinates are real; the sequence is deliberately N, because
        the part sequences belong to your plasmids and not to this tool. Drop in the
        real elements in your cloning software and the coordinates will hold.
        """
        total = self.total_bp
        lines = [
            f"LOCUS       {name[:16]:<16} {total} bp    DNA     linear   SYN",
            f"DEFINITION  cargologic design: {len(self.cassettes)} cassette(s), "
            f"strategy={self.strategy}.",
            "ACCESSION   .",
            "KEYWORDS    recombinase logic; synthetic circuit.",
            "FEATURES             Location/Qualifiers",
        ]
        for row in self.feature_table():
            loc = f"{row['start']}..{row['end']}"
            key = {
                "promoter": "promoter",
                "site": "misc_recomb",
                "terminator": "terminator",
                "payload": "CDS",
                "linker": "misc_feature",
                "kozak": "misc_feature",
                "itr": "repeat_region",
                "spacer": "misc_feature",
            }.get(row["kind"], "misc_feature")
            lines.append(f"     {key:<15} {loc}")
            lines.append(f"                     /label=\"{row['name']}\"")
            if row["input"]:
                lines.append(f"                     /note=\"input {row['input']}\"")
            if row["note"]:
                lines.append(f"                     /note=\"{row['note']}\"")
        lines.append("ORIGIN")
        seq = "n" * total
        for i in range(0, total, 60):
            chunk = seq[i : i + 60]
            groups = " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10))
            lines.append(f"{i + 1:>9} {groups}")
        lines.append("//")
        return "\n".join(lines)

    def ascii_map(self, width: int = 78) -> str:
        """A compact linear schematic, one line per cassette."""
        out = []
        for ci, c in enumerate(self.cassettes, 1):
            pieces = []
            for e in c.elements:
                if e.kind == "site":
                    pieces.append(f"<{e.name}>")
                elif e.kind == "promoter":
                    pieces.append(f"[{e.name}>")
                elif e.kind == "terminator":
                    pieces.append(f"|{e.name}|")
                elif e.kind == "payload":
                    pieces.append(f"({e.name})")
                elif e.kind in {"kozak", "spacer"}:
                    continue
                else:
                    pieces.append(f"-{e.name}-")
            term = c.implicant.to_str(self.inputs)
            out.append(f"  cassette {ci}  [{term}]  {c.bp} bp")
            line = "    " + "-".join(pieces)
            while len(line) > width:
                out.append(line[:width])
                line = "      " + line[width:]
            out.append(line)
        if self.shared_payload_elements:
            pieces = []
            for e in self.shared_payload_elements:
                if e.kind == "site":
                    pieces.append(f"<{e.name}>")
                elif e.kind == "payload":
                    pieces.append(f"({e.name})")
                elif e.kind == "terminator":
                    pieces.append(f"|{e.name}|")
                elif e.kind == "linker":
                    pieces.append(f"-{e.name}-")
            out.append(
                f"  shared payload  {sum(e.bp for e in self.shared_payload_elements)} bp"
            )
            out.append("    " + "-".join(pieces))
        if self.recombinase_elements:
            out.append(f"  recombinase expression  {self.recombinase_bp} bp")
            out.append(
                "    "
                + "-".join(
                    f"({e.name})" if e.kind == "payload" else f"[{e.name}>"
                    for e in self.recombinase_elements
                    if e.kind in {"payload", "promoter"}
                )
            )
        return "\n".join(out)


# --------------------------------------------------------------------------------------
# Layout construction
# --------------------------------------------------------------------------------------


def _sites_for(rec: RecombinaseSystem, input_index: int, role: str) -> tuple[Element, Element]:
    """A matched pair of recombinase sites."""
    left = Element(
        "site", rec.site_name, rec.site_bp // 2 + rec.site_bp % 2,
        f"{role}, 5' site", input_index,
    )
    right = Element(
        "site", rec.site_name, rec.site_bp // 2, f"{role}, 3' site", input_index
    )
    return left, right


def build_cassette(
    implicant: Implicant,
    outputs: list[str],
    payloads: list[Payload],
    promoter: Promoter,
    blocker: Terminator,
    terminator: Terminator,
    recombinases: list[RecombinaseSystem],
    linker: Linker | None,
    carries_payload: bool = True,
) -> Cassette:
    """Lay out one product term."""
    pos_lits = implicant.positive_literals()
    neg_lits = implicant.negative_literals()
    els: list[Element] = []

    els.append(Element("promoter", promoter.name, promoter.bp, "drives this term"))

    # AND over positive literals: a blocker per literal, all must be excised.
    for i in pos_lits:
        rec = recombinases[i]
        l, r = _sites_for(rec, i, f"blocker gate for {rec.name}")
        els.append(l)
        els.append(
            Element(
                "terminator",
                blocker.name,
                blocker.bp,
                f"blocks transcription until {rec.name} fires",
                i,
            )
        )
        els.append(r)
        els.append(Element("spacer", "junction", JUNCTION_BP))

    if carries_payload:
        # NOT over negative literals: nested site pairs flanking the payload.
        opens: list[Element] = []
        closes: list[Element] = []
        for i in neg_lits:
            rec = recombinases[i]
            l, r = _sites_for(rec, i, f"payload excision by {rec.name}")
            opens.append(l)
            closes.insert(0, r)
        els.extend(opens)

        for pi, pay in enumerate(payloads):
            els.append(Element("kozak", "Kozak", KOZAK_BP))
            els.append(
                Element("payload", pay.name, pay.bp, outputs[pi] if pi < len(outputs) else "")
            )
            if linker is not None and pi < len(payloads) - 1:
                els.append(Element("linker", linker.name, linker.bp, "co-expression"))

        els.extend(closes)
        els.append(Element("spacer", "junction", JUNCTION_BP))

    els.append(
        Element("terminator", terminator.name, terminator.bp, "transcription stop")
    )
    els.append(Element("spacer", "junction", JUNCTION_BP))

    return Cassette(
        outputs=list(outputs),
        payloads=list(payloads),
        promoter=promoter,
        implicant=implicant,
        positive_literals=pos_lits,
        negative_literals=neg_lits,
        blocker=blocker,
        terminator=terminator,
        elements=els,
        carries_payload=carries_payload,
        linker=linker,
    )


def can_share_payload(cover: Cover) -> bool:
    """True when every product term has the same negative-literal set.

    That is the condition for parking one payload copy downstream of several gated
    promoters: the sites that would excise the payload have to be identical for all
    terms, or one term's NOT would silence another term's ON state.
    """
    if cover.n_terms <= 1:
        return True
    sets = {tuple(sorted(t.negative_literals())) for t in cover.terms}
    return len(sets) == 1


def build_shared_payload_tail(
    neg_lits: list[int],
    outputs: list[str],
    payloads: list[Payload],
    terminator: Terminator,
    recombinases: list[RecombinaseSystem],
    linker: Linker | None,
) -> list[Element]:
    """The single payload copy shared by tandem gated promoters."""
    els: list[Element] = []
    opens: list[Element] = []
    closes: list[Element] = []
    for i in neg_lits:
        rec = recombinases[i]
        l, r = _sites_for(rec, i, f"payload excision by {rec.name}")
        opens.append(l)
        closes.insert(0, r)
    els.extend(opens)
    for pi, pay in enumerate(payloads):
        els.append(Element("kozak", "Kozak", KOZAK_BP))
        els.append(
            Element("payload", pay.name, pay.bp, outputs[pi] if pi < len(outputs) else "")
        )
        if linker is not None and pi < len(payloads) - 1:
            els.append(Element("linker", linker.name, linker.bp, "co-expression"))
    els.extend(closes)
    els.append(Element("spacer", "junction", JUNCTION_BP))
    els.append(Element("terminator", terminator.name, terminator.bp, "transcription stop"))
    return els


def build_recombinase_elements(
    recombinases: list[RecombinaseSystem],
    promoter: Promoter,
    terminator: Terminator,
    linker: Linker,
) -> list[Element]:
    """Self-contained recombinase expression, when the circuit must carry its own inputs.

    Usually you do NOT want this: in published devices the recombinases are supplied
    separately (co-transfection, a second vector, or driven by sensor modules). Including
    them adds 1-2 kb each, which blows an AAV budget immediately -- and the tool says so
    rather than quietly producing something unbuildable.
    """
    els: list[Element] = [
        Element("promoter", promoter.name, promoter.bp, "drives recombinases")
    ]
    for i, rec in enumerate(recombinases):
        els.append(Element("kozak", "Kozak", KOZAK_BP))
        els.append(
            Element("payload", rec.name, rec.cds_bp, f"{rec.name} recombinase CDS")
        )
        if i < len(recombinases) - 1:
            els.append(Element("linker", linker.name, linker.bp, "co-expression"))
    els.append(Element("spacer", "junction", JUNCTION_BP))
    els.append(Element("terminator", terminator.name, terminator.bp, ""))
    return els
