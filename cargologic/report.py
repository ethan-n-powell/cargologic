"""Rendering the verdict.

The report is the product. It leads with the answer, separates what is arithmetic from
what is modelled, and names the specific state that will fail rather than handing back a
single score.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from .compiler import CompileResult

_SEV = {"ok": "FITS", "tight": "FITS (TIGHT)", "over": "DOES NOT FIT", "infeasible": "NOT APPLICABLE"}


def _rule(width: int, ch: str = "-") -> str:
    return ch * width


def _header(title: str, width: int) -> str:
    return f"\n{title}\n{_rule(min(width, len(title)))}"


def render(result: "CompileResult", width: int = 80) -> str:
    out: list[str] = []
    d = result.design
    tv = result.target_verdict

    # ---------------------------------------------------------------- headline
    out.append(_rule(width, "="))
    verdict = _SEV.get(tv.severity, tv.severity.upper())
    out.append(f"  {verdict}  -  {tv.vector}")
    out.append(f"  {d.total_bp} bp circuit, {d.n_terms} cassette(s)")
    if tv.budget_bp:
        out.append(
            f"  {tv.headroom_bp:+d} bp against a {tv.budget_bp} bp budget "
            f"({100 * (tv.utilisation or 0):.0f}% used)"
        )
    primary = result.spec.functions[0].name
    if primary in result.leakage:
        rep = result.leakage[primary]
        if rep.dynamic_range_median == float("inf"):
            out.append("  predicted ON/OFF: no leak path modelled")
        else:
            out.append(
                f"  predicted ON/OFF ratio ~{rep.dynamic_range_median:.0f}x "
                f"(80% interval {rep.dynamic_range_lo:.0f}-{rep.dynamic_range_hi:.0f}x)"
            )
    out.append(_rule(width, "="))

    # ---------------------------------------------------------------- spec
    out.append(_header("SPECIFICATION", width))
    out.append(f"inputs      {', '.join(result.spec.inputs)}")
    for f in result.spec.functions:
        pay = result.spec.payloads[f.name]
        out.append(
            f"output {f.name:<4} {f.source or 'table'}   "
            f"-> {pay.name} ({pay.bp} bp)"
        )
    out.append("")
    out.append("recombinase assignment")
    for i, nm in enumerate(result.spec.inputs):
        rec = d.recombinases[i]
        out.append(
            f"  {nm} -> {rec.name:<7} sites {rec.site_name:<10} "
            f"{rec.site_bp} bp/pair, {rec.family}"
        )

    # truth table
    for f in result.spec.functions:
        out.append("")
        hdr = "  " + " ".join(result.spec.inputs) + f"  |  {f.name}"
        out.append(hdr)
        out.append("  " + _rule(len(hdr) - 2))
        for bits, val in f.truth_rows():
            spaced = " ".join(bits)
            out.append(f"  {spaced}  |  {'1' if val else '0'}")

    # ---------------------------------------------------------------- logic
    out.append(_header("LOGIC MINIMISATION", width))
    out.append(
        "Each product term becomes one physical cassette, so term count is the\n"
        "dominant size lever -- far bigger than any part substitution."
    )
    for name, cov in result.covers.items():
        out.append("")
        out.append(f"  {name} = {cov.to_str(result.spec.inputs)}")
        out.append(
            f"    {cov.n_terms} term(s), {cov.n_positive_literals} positive literal(s) "
            f"-> {cov.n_positive_literals} blocker(s), "
            f"{cov.n_literals - cov.n_positive_literals} negative literal(s)"
        )
        for t in cov.terms:
            pos = [result.spec.inputs[i] for i in t.positive_literals()]
            neg = [result.spec.inputs[i] for i in t.negative_literals()]
            bits = []
            if pos:
                bits.append(f"blockers gated by {', '.join(pos)}")
            if neg:
                bits.append(f"payload excised by {', '.join(neg)}")
            out.append(
                f"      [{t.pattern()}] {t.to_str(result.spec.inputs):<22} "
                + ("; ".join(bits) if bits else "ungated")
            )

    # ---------------------------------------------------------------- layout
    out.append(_header("ARCHITECTURE", width))
    out.append(f"strategy    {d.strategy}")
    out.append("")
    out.append(d.ascii_map(width=width - 2))
    out.append("")
    out.append("  legend  [promoter>  <recombinase site>  |terminator|  (payload)")

    out.append("")
    out.append("size breakdown")
    by_kind: dict[str, int] = {}
    for e in d.elements:
        by_kind[e.kind] = by_kind.get(e.kind, 0) + e.bp
    for kind, bp in sorted(by_kind.items(), key=lambda kv: -kv[1]):
        pct = 100 * bp / d.total_bp if d.total_bp else 0
        bar = "#" * int(pct / 2.5)
        out.append(f"  {kind:<12} {bp:>6} bp  {pct:>5.1f}%  {bar}")
    out.append(f"  {'TOTAL':<12} {d.total_bp:>6} bp")
    if d.include_recombinases:
        out.append(
            f"  (of which {d.recombinase_bp} bp is recombinase expression)"
        )

    out.append("")
    out.append(f"parts       promoter {result.choice.promoter}, "
               f"blocker {result.choice.blocker}, terminator {result.choice.terminator}")
    if result.searched:
        out.append(
            f"            selected by searching {len(result.candidates)} part "
            f"combinations"
        )

    # ---------------------------------------------------------------- delivery
    out.append(_header("DELIVERY", width))
    out.append("This section is arithmetic over physical limits. Trust it.")
    out.append("")
    for v in result.all_verdicts:
        mark = {"ok": "  ok  ", "tight": " tight", "over": " OVER ", "infeasible": "  n/a "}.get(
            v.severity, "      "
        )
        out.append(f"  [{mark}] {v.summary()}")
    out.append("")
    for v in result.all_verdicts:
        if v.notes and v.severity in {"over", "tight", "infeasible"}:
            out.append(f"  {v.vector}")
            for nt in v.notes:
                for line in _wrap(nt, width - 6):
                    out.append(f"      {line}")

    # ---------------------------------------------------------------- behaviour
    out.append(_header("PREDICTED BEHAVIOUR (MODEL -- SEE CAVEATS)", width))
    for name, rep in result.leakage.items():
        out.append("")
        out.append(f"  {name}")
        hdr = "    state  intended   predicted level (median [10-90%])"
        out.append(hdr)
        out.append("    " + _rule(len(hdr) - 4))
        peak = max((s.level_median for s in rep.states), default=1.0) or 1.0
        for s in rep.states:
            bar = "#" * int(24 * s.level_median / peak)
            flag = ""
            if not s.intended and s.level_median > 0.1 * rep.on_level_median:
                flag = "  <-- leaky"
            out.append(
                f"    {' '.join(s.bits)}   {s.intended_str:<8}  "
                f"{s.level_median:6.3f} [{s.level_lo:5.3f}-{s.level_hi:5.3f}]  "
                f"{bar}{flag}"
            )
        if rep.dynamic_range_median != float("inf"):
            out.append(
                f"    ON/OFF  {rep.dynamic_range_median:.1f}x  "
                f"(80% interval {rep.dynamic_range_lo:.1f}-{rep.dynamic_range_hi:.1f}x)"
            )
            out.append(
                f"    weakest ON state {rep.weakest_on_state}, "
                f"leakiest OFF state {rep.worst_off_state}"
            )

        if rep.mechanisms:
            out.append("")
            out.append(
                f"    what is limiting state {rep.worst_off_state} "
                f"(share of the excess each mechanism accounts for)"
            )
            for m in rep.mechanisms:
                bar = "#" * int(28 * max(0.0, m.improvement))
                out.append(
                    f"      {m.name:<34} {100 * m.improvement:>5.1f}%  {bar}"
                )
            if rep.limiting:
                out.append("")
                for line in _wrap(f"-> {rep.limiting.advice}", width - 8):
                    out.append(f"      {line}")
            out.append("")
            for line in _wrap(
                "This ranking is more trustworthy than the absolute levels above: it "
                "compares mechanisms against each other, so it survives the parameter "
                "estimates being off by a factor.",
                width - 8,
            ):
                out.append(f"      {line}")

        for w in rep.warnings:
            out.append("")
            for line in _wrap("! " + w, width - 6):
                out.append(f"    {line}")

    # ---------------------------------------------------------------- pareto
    front = result.pareto
    if front and len(front) > 1:
        out.append(_header("SIZE vs DYNAMIC RANGE (PARETO FRONT)", width))
        out.append(
            "Smaller polyA signals save bp and leak more. There is no single best\n"
            "design near the packaging limit, so here is the trade-off."
        )
        out.append("")
        out.append("     bp   ON/OFF   fits   promoter  blocker   terminator")
        out.append("  " + _rule(58))
        for c in front[:12]:
            dr = "inf" if c.dynamic_range >= 1e9 else f"{c.dynamic_range:.0f}x"
            out.append(
                f"  {c.circuit_bp:>5}  {dr:>6}   {'y' if c.fits else 'n':^4}   "
                f"{c.choice.promoter:<9} {c.choice.blocker:<9} {c.choice.terminator}"
            )

    # ---------------------------------------------------------------- notes
    if d.notes:
        out.append(_header("DESIGN NOTES", width))
        for n in d.notes:
            for line in _wrap("- " + n, width - 2):
                out.append(f"  {line}")

    # ---------------------------------------------------------------- caveats
    out.append(_header("CAVEATS", width))
    for c in CAVEATS:
        for line in _wrap("- " + c, width - 2):
            out.append(f"  {line}")

    return "\n".join(out)


CAVEATS = [
    "Sizes are reliable. Predicted expression levels are not: they come from a model "
    "built on literature parameter ranges, and recent work shows mechanistic models of "
    "mammalian circuits underperform badly on held-out designs. Use the predictions to "
    "rank designs and to find the state that will fail, not as expected measurements.",
    "Published quantitative part data is almost entirely HEK293 transient transfection. "
    "If your target is primary T cells, NK cells or iPSCs, treat every efficiency and "
    "leak number here as a guess, and expect promoter silencing that this tool does not "
    "model.",
    "Not modelled: transcriptional squelching at high activator dose, competition for "
    "shared transcriptional and degradation resources, retroactivity, chromosomal "
    "position effects, and transgene silencing. Each of these has broken a real design.",
    "Recombinase orthogonality is assumed from published pairings. Verify site "
    "cross-reactivity experimentally before committing to a multi-input build.",
    "The emitted GenBank has real coordinates and N-masked sequence. Substitute actual "
    "part sequences from your own plasmids or Addgene; this tool does not invent "
    "sequence.",
    "Tandem-promoter sharing assumes transcription initiating upstream still yields "
    "functional payload mRNA after reading through a downstream promoter. That is "
    "usually true and occasionally is not.",
]


def _wrap(text: str, width: int) -> list[str]:
    words = text.split()
    lines: list[str] = []
    cur = ""
    for w in words:
        if not cur:
            cur = w
        elif len(cur) + 1 + len(w) <= width:
            cur += " " + w
        else:
            lines.append(cur)
            cur = "  " + w if text.startswith(("- ", "! ")) else w
    if cur:
        lines.append(cur)
    return lines
