"""Turns a logic cover into a physical DNA layout with a length."""

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
    # Set by layouts whose element order is not cassette-by-cassette, such as the
    # BLADE template where one promoter sits upstream of every address.
    ordered_elements: list[Element] | None = None

    @property
    def elements(self) -> list[Element]:
        if self.ordered_elements is not None:
            return self.ordered_elements + self.recombinase_elements
        return (
            [e for c in self.cassettes for e in c.elements]
            + self.shared_payload_elements
            + self.recombinase_elements
        )

    @property
    def circuit_bp(self) -> int:
        """Length of the logic circuit itself, excluding recombinase expression."""
        if self.ordered_elements is not None:
            return sum(e.bp for e in self.ordered_elements)
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

    def feature_table(self) -> list[dict]:
        """Coordinates for every element, 1-based inclusive."""
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
        """GenBank feature table with real coordinates over an N-masked sequence."""
        total = self.total_bp
        lines = [
            f"LOCUS       {name[:16]:<16} {total} bp    DNA     linear   SYN",
            f"DEFINITION  cargologic design: {len(self.cassettes)} cassette(s), "
            f"strategy={self.strategy}.",
            "ACCESSION   .",
            "KEYWORDS    recombinase logic; synthetic circuit.",
            "FEATURES             Location/Qualifiers",
        ]
        feature_keys = {
            "promoter": "promoter",
            "site": "misc_recomb",
            "terminator": "terminator",
            "payload": "CDS",
            "linker": "misc_feature",
            "kozak": "misc_feature",
            "itr": "repeat_region",
            "spacer": "misc_feature",
        }
        for row in self.feature_table():
            key = feature_keys.get(row["kind"], "misc_feature")
            lines.append(f"     {key:<15} {row['start']}..{row['end']}")
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
        fmt = {
            "site": "<{}>", "promoter": "[{}>", "terminator": "|{}|",
            "payload": "({})", "linker": "-{}-",
        }
        out = []
        for ci, c in enumerate(self.cassettes, 1):
            pieces = [
                fmt.get(e.kind, "-{}-").format(e.name)
                for e in c.elements
                if e.kind not in {"kozak", "spacer"}
            ]
            out.append(
                f"  cassette {ci}  [{c.implicant.to_str(self.inputs)}]  {c.bp} bp"
            )
            line = "    " + "-".join(pieces)
            while len(line) > width:
                out.append(line[:width])
                line = "      " + line[width:]
            out.append(line)
        if self.shared_payload_elements:
            shared_kinds = {"site", "payload", "terminator", "linker"}
            out.append(
                f"  shared payload  {sum(e.bp for e in self.shared_payload_elements)} bp"
            )
            out.append(
                "    "
                + "-".join(
                    fmt[e.kind].format(e.name)
                    for e in self.shared_payload_elements
                    if e.kind in shared_kinds
                )
            )
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


def _sites_for(rec: RecombinaseSystem, input_index: int, role: str) -> tuple[Element, Element]:
    half = rec.site_bp // 2
    return (
        Element("site", rec.site_name, half + rec.site_bp % 2, f"{role}, 5' site", input_index),
        Element("site", rec.site_name, half, f"{role}, 3' site", input_index),
    )


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
    els: list[Element] = [
        Element("promoter", promoter.name, promoter.bp, "drives this term")
    ]

    # AND over positive literals: every blocker must be excised to transcribe.
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
        pairs: list[tuple[Element, Element]] = []
        for i in neg_lits:
            rec = recombinases[i]
            pairs.append(_sites_for(rec, i, f"payload excision by {rec.name}"))
        els.extend(l for l, _ in pairs)

        for pi, pay in enumerate(payloads):
            els.append(Element("kozak", "Kozak", KOZAK_BP))
            els.append(
                Element("payload", pay.name, pay.bp, outputs[pi] if pi < len(outputs) else "")
            )
            if linker is not None and pi < len(payloads) - 1:
                els.append(Element("linker", linker.name, linker.bp, "co-expression"))

        els.extend(r for _, r in reversed(pairs))
        els.append(Element("spacer", "junction", JUNCTION_BP))

    els.append(Element("terminator", terminator.name, terminator.bp, "transcription stop"))
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
    """True when every product term has the same negative-literal set."""
    if cover.n_terms <= 1:
        return True
    return len({tuple(sorted(t.negative_literals())) for t in cover.terms}) == 1


def build_shared_payload_tail(
    neg_lits: list[int],
    outputs: list[str],
    payloads: list[Payload],
    terminator: Terminator,
    recombinases: list[RecombinaseSystem],
    linker: Linker | None,
) -> list[Element]:
    """The single payload copy shared by tandem gated promoters."""
    pairs: list[tuple[Element, Element]] = []
    for i in neg_lits:
        rec = recombinases[i]
        pairs.append(_sites_for(rec, i, f"payload excision by {rec.name}"))
    els: list[Element] = [l for l, _ in pairs]
    for pi, pay in enumerate(payloads):
        els.append(Element("kozak", "Kozak", KOZAK_BP))
        els.append(
            Element("payload", pay.name, pay.bp, outputs[pi] if pi < len(outputs) else "")
        )
        if linker is not None and pi < len(payloads) - 1:
            els.append(Element("linker", linker.name, linker.bp, "co-expression"))
    els.extend(r for _, r in reversed(pairs))
    els.append(Element("spacer", "junction", JUNCTION_BP))
    els.append(Element("terminator", terminator.name, terminator.bp, "transcription stop"))
    return els


def build_recombinase_elements(
    recombinases: list[RecombinaseSystem],
    promoter: Promoter,
    terminator: Terminator,
    linker: Linker,
) -> list[Element]:
    """Recombinase expression unit, for circuits that must carry their own inputs.

    Published devices supply the recombinases separately; each CDS added here costs
    1-2 kb, which exhausts an AAV budget on its own.
    """
    els: list[Element] = [
        Element("promoter", promoter.name, promoter.bp, "drives recombinases")
    ]
    for i, rec in enumerate(recombinases):
        els.append(Element("kozak", "Kozak", KOZAK_BP))
        els.append(Element("payload", rec.name, rec.cds_bp, f"{rec.name} recombinase CDS"))
        if i < len(recombinases) - 1:
            els.append(Element("linker", linker.name, linker.bp, "co-expression"))
    els.append(Element("spacer", "junction", JUNCTION_BP))
    els.append(Element("terminator", terminator.name, terminator.bp, ""))
    return els


def build_blade(
    n_inputs: int,
    on_minterms: set[int],
    outputs: list[str],
    payloads: list[Payload],
    promoter: Promoter,
    terminator: Terminator,
    recombinases: list[RecombinaseSystem],
    linker: Linker | None,
) -> tuple[list[Cassette], list[Element]]:
    """Lay out the BLADE template: one promoter upstream of 2^N addresses.

    Weinberg et al. 2017 state the structure directly: a single transcriptional layer
    with one promoter upstream of up to 2^N addresses surrounded by recombination
    sites, and 2^N-1 recombination site pairs. Every address is physically present;
    an OFF address holds a terminator alone.

    Wins over product terms when the cover cannot factor, parity being the clearest
    case, because it pays one promoter rather than one per term.
    """
    n_addresses = 1 << n_inputs
    full_mask = n_addresses - 1
    ordered: list[Element] = [
        Element("promoter", promoter.name, promoter.bp, "shared by every address")
    ]
    cassettes: list[Cassette] = []

    for address in range(n_addresses):
        bits = format(address, f"0{n_inputs}b")
        if address > 0:
            # 2^N-1 site pairs in total, one boundary between consecutive addresses.
            i = (address - 1) % n_inputs
            rec = recombinases[i]
            left, right = _sites_for(rec, i, f"address {bits} boundary")
            ordered.append(left)
            ordered.append(right)

        if address in on_minterms:
            els: list[Element] = []
            for pi, pay in enumerate(payloads):
                els.append(Element("kozak", "Kozak", KOZAK_BP))
                els.append(
                    Element(
                        "payload",
                        pay.name,
                        pay.bp,
                        outputs[pi] if pi < len(outputs) else "",
                    )
                )
                if linker is not None and pi < len(payloads) - 1:
                    els.append(Element("linker", linker.name, linker.bp, "co-expression"))
            els.append(
                Element("terminator", terminator.name, terminator.bp, f"address {bits}")
            )
            ordered.extend(els)
            implicant = Implicant(bits=address, mask=full_mask, n=n_inputs)
            cassettes.append(
                Cassette(
                    outputs=list(outputs),
                    payloads=list(payloads),
                    promoter=promoter,
                    implicant=implicant,
                    positive_literals=implicant.positive_literals(),
                    negative_literals=implicant.negative_literals(),
                    blocker=terminator,
                    terminator=terminator,
                    elements=els,
                    carries_payload=True,
                    linker=linker,
                )
            )
        else:
            ordered.append(
                Element(
                    "terminator",
                    terminator.name,
                    terminator.bp,
                    f"address {bits} is OFF",
                )
            )
        ordered.append(Element("spacer", "junction", JUNCTION_BP))

    return cassettes, ordered
