#!/usr/bin/env python3
"""Retrospective benchmark: does the compiler agree with published circuits?

WHY THIS FILE IS THE POINT OF THE PROJECT
-----------------------------------------
Anyone can write a circuit compiler. Nobody has published an evaluation of whether
model-based mammalian circuit design actually predicts which published circuits worked.
The ~113 circuits in the BLADE paper, the ~189 in CREATE, and the 256 predicted by the
BLADE mechanistic model are a free, already-paid-for test set, and no wet lab is needed
to use them.

The result is publishable whichever way it falls. If the compiler ranks working circuits
above broken ones, that is a validated design tool. If it does not, that is a quantified
negative result on mammalian composability -- which is arguably the more useful paper,
and it is the honest outcome to be prepared for.

WHAT IS SCORED
--------------
1. Size accuracy      Computed length vs the construct length reported in the paper.
                      This checks arithmetic and should be near-exact. If it is not,
                      the part library is wrong before anything else is.
2. Feasibility recall Of circuits that were built and worked, how many does the
                      compiler call feasible? Should be ~100%; anything less is a
                      compiler bug, not a biology finding.
3. Discrimination     Does predicted dynamic range rank working circuits above failed
                      ones? Reported as AUC. 0.5 is chance.
4. Calibration        Spearman correlation between predicted and measured fold-change.
                      Expect this to be the weakest number, and report it anyway.

THE DATA IS NOT IN THIS REPOSITORY
----------------------------------
See README.md for exactly where to get each dataset. `circuits.csv` ships as a header
only. `synthetic_smoke.csv` contains FABRICATED rows whose only purpose is to exercise
this script; never report a number computed from it.

    python validation/benchmark.py validation/circuits.csv
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import cargologic as cl  # noqa: E402
from cargologic.logic import LogicError  # noqa: E402


@dataclass
class Record:
    """One published circuit."""

    circuit_id: str
    source: str
    truth_table: str
    inputs: list[str]
    payload_bp: int
    recombinases: list[str] | None
    outcome: str  # worked | failed | partial | unknown
    reported_fold_change: float | None
    reported_construct_bp: int | None
    vector: str
    notes: str = ""

    @property
    def worked(self) -> bool | None:
        if self.outcome == "worked":
            return True
        if self.outcome == "failed":
            return False
        return None


@dataclass
class Scored:
    record: Record
    predicted_bp: int
    predicted_fold_change: float | None
    feasible: bool
    n_terms: int
    limiting: str | None
    error: str | None = None


def load(path: Path) -> list[Record]:
    rows: list[Record] = []
    with path.open() as fh:
        for raw in csv.DictReader(fh):
            if not raw.get("circuit_id") or raw["circuit_id"].lstrip().startswith("#"):
                continue

            def opt_float(key: str) -> float | None:
                v = (raw.get(key) or "").strip()
                return float(v) if v else None

            def opt_int(key: str) -> int | None:
                v = (raw.get(key) or "").strip()
                return int(v) if v else None

            inputs = [s.strip() for s in (raw.get("inputs") or "").split(";") if s.strip()]
            recs = [s.strip() for s in (raw.get("recombinases") or "").split(";") if s.strip()]
            rows.append(
                Record(
                    circuit_id=raw["circuit_id"].strip(),
                    source=(raw.get("source") or "").strip(),
                    truth_table=(raw.get("truth_table") or "").strip(),
                    inputs=inputs or None,
                    payload_bp=int(raw.get("payload_bp") or 720),
                    recombinases=recs or None,
                    outcome=(raw.get("outcome") or "unknown").strip().lower(),
                    reported_fold_change=opt_float("reported_fold_change"),
                    reported_construct_bp=opt_int("reported_construct_bp"),
                    vector=(raw.get("vector") or "plasmid").strip() or "plasmid",
                    notes=(raw.get("notes") or "").strip(),
                )
            )
    return rows


def score_one(rec: Record, n_samples: int = 600) -> Scored:
    try:
        result = cl.compile_circuit(
            rec.truth_table,
            payload=rec.payload_bp,
            vector=rec.vector,
            inputs=rec.inputs,
            recombinases=rec.recombinases,
            as_table=True,
            optimise=False,
            n_samples=n_samples,
        )
    except (LogicError, cl.SynthesisError, KeyError, ValueError) as exc:
        return Scored(rec, 0, None, False, 0, None, error=str(exc))

    rep = next(iter(result.leakage.values()))
    dr = rep.dynamic_range_median
    return Scored(
        record=rec,
        predicted_bp=result.design.total_bp,
        predicted_fold_change=None if dr == float("inf") else dr,
        feasible=result.feasible,
        n_terms=next(iter(result.covers.values())).n_terms,
        limiting=rep.limiting.name if rep.limiting else None,
    )


# ------------------------------------------------------------------ statistics, no deps


def spearman(xs: list[float], ys: list[float]) -> float | None:
    """Rank correlation. Returns None if there is not enough variation to compute one."""
    if len(xs) < 3:
        return None

    def ranks(vs: list[float]) -> list[float]:
        order = sorted(range(len(vs)), key=lambda i: vs[i])
        out = [0.0] * len(vs)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and vs[order[j + 1]] == vs[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                out[order[k]] = avg
            i = j + 1
        return out

    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def auc(positives: list[float], negatives: list[float]) -> float | None:
    """Mann-Whitney U / ROC AUC. Probability a random positive outranks a random negative."""
    if not positives or not negatives:
        return None
    wins = 0.0
    for p, n in itertools.product(positives, negatives):
        wins += 1.0 if p > n else (0.5 if p == n else 0.0)
    return wins / (len(positives) * len(negatives))


# ------------------------------------------------------------------------------ report


def report(scored: list[Scored], warn_synthetic: bool = False) -> str:
    out: list[str] = []
    ok = [s for s in scored if s.error is None]
    errs = [s for s in scored if s.error is not None]

    out.append("=" * 76)
    out.append("  cargologic retrospective benchmark")
    out.append("=" * 76)
    if warn_synthetic:
        out.append("")
        out.append("  *** SYNTHETIC DATA -- FABRICATED FOR SMOKE TESTING ONLY ***")
        out.append("  *** Do not report any number below as a result.          ***")
    out.append("")
    out.append(f"  circuits loaded      {len(scored)}")
    out.append(f"  compiled             {len(ok)}")
    out.append(f"  failed to compile    {len(errs)}")

    sources = sorted({s.record.source for s in scored if s.record.source})
    if sources:
        out.append(f"  sources              {', '.join(sources)}")

    if errs:
        out.append("")
        out.append("  COMPILE FAILURES (compiler bugs, not biology)")
        for s in errs[:10]:
            out.append(f"    {s.record.circuit_id:<16} {s.error}")

    # ---------------------------------------------------------------- 1. size accuracy
    sized = [s for s in ok if s.record.reported_construct_bp]
    out.append("")
    out.append("-" * 76)
    out.append("  1. SIZE ACCURACY  (arithmetic -- should be near exact)")
    out.append("-" * 76)
    if not sized:
        out.append("  no reported construct lengths in this dataset")
    else:
        errors = [s.predicted_bp - s.record.reported_construct_bp for s in sized]
        abs_err = [abs(e) for e in errors]
        mae = sum(abs_err) / len(abs_err)
        rel = [
            abs(e) / s.record.reported_construct_bp
            for e, s in zip(errors, sized)
        ]
        out.append(f"  n                    {len(sized)}")
        out.append(f"  mean absolute error  {mae:.0f} bp")
        out.append(f"  median rel. error    {sorted(rel)[len(rel) // 2] * 100:.1f}%")
        out.append(f"  worst                {max(abs_err)} bp")
        out.append("")
        out.append(f"  {'circuit':<18}{'reported':>10}{'predicted':>11}{'delta':>9}")
        for s in sorted(sized, key=lambda x: -abs(x.predicted_bp - x.record.reported_construct_bp))[:10]:
            d = s.predicted_bp - s.record.reported_construct_bp
            out.append(
                f"  {s.record.circuit_id:<18}{s.record.reported_construct_bp:>10}"
                f"{s.predicted_bp:>11}{d:>+9}"
            )

    # ------------------------------------------------------- 2. feasibility recall
    worked = [s for s in ok if s.record.worked is True]
    out.append("")
    out.append("-" * 76)
    out.append("  2. FEASIBILITY RECALL  (of circuits that worked, how many do we allow?)")
    out.append("-" * 76)
    if not worked:
        out.append("  no circuits labelled as working")
    else:
        n_ok = sum(1 for s in worked if s.feasible)
        out.append(f"  {n_ok}/{len(worked)} = {100 * n_ok / len(worked):.1f}%")
        missed = [s for s in worked if not s.feasible]
        if missed:
            out.append("")
            out.append("  called INFEASIBLE despite having been built and worked:")
            for s in missed[:10]:
                out.append(
                    f"    {s.record.circuit_id:<16} {s.predicted_bp} bp "
                    f"in {s.record.vector}  <- investigate: likely a part-size error"
                )

    # ------------------------------------------------------------ 3. discrimination
    pos = [
        s.predicted_fold_change
        for s in ok
        if s.record.worked is True and s.predicted_fold_change is not None
    ]
    neg = [
        s.predicted_fold_change
        for s in ok
        if s.record.worked is False and s.predicted_fold_change is not None
    ]
    out.append("")
    out.append("-" * 76)
    out.append("  3. DISCRIMINATION  (does predicted dynamic range separate the two?)")
    out.append("-" * 76)
    a = auc(pos, neg)
    if a is None:
        out.append(
            f"  not computable: {len(pos)} working and {len(neg)} failed circuits with "
            f"a finite prediction. Need both."
        )
    else:
        out.append(f"  n working {len(pos)}, n failed {len(neg)}")
        out.append(f"  AUC                  {a:.3f}   (0.5 = chance)")
        if pos:
            out.append(f"  median predicted DR, worked   {sorted(pos)[len(pos) // 2]:.1f}x")
        if neg:
            out.append(f"  median predicted DR, failed   {sorted(neg)[len(neg) // 2]:.1f}x")
        verdict = (
            "no better than chance" if a < 0.6
            else "weak" if a < 0.7
            else "useful" if a < 0.85
            else "strong"
        )
        out.append(f"  interpretation       {verdict}")

    # --------------------------------------------------------------- 4. calibration
    paired = [
        (s.predicted_fold_change, s.record.reported_fold_change)
        for s in ok
        if s.predicted_fold_change is not None and s.record.reported_fold_change
    ]
    out.append("")
    out.append("-" * 76)
    out.append("  4. CALIBRATION  (predicted vs measured fold-change)")
    out.append("-" * 76)
    if len(paired) < 3:
        out.append(f"  not computable: only {len(paired)} circuits have both values")
    else:
        rho = spearman([p for p, _ in paired], [m for _, m in paired])
        out.append(f"  n                    {len(paired)}")
        out.append(f"  Spearman rho         {rho:.3f}" if rho is not None else "  Spearman rho  n/a")
        ratios = [p / m for p, m in paired if m]
        if ratios:
            ratios.sort()
            out.append(
                f"  predicted/measured   median {ratios[len(ratios) // 2]:.2f}x, "
                f"range {ratios[0]:.2f}-{ratios[-1]:.2f}x"
            )
            out.append(
                "  A median far from 1.0 means the model is biased, which is fixable. "
                "A low rho"
            )
            out.append(
                "  means it is uninformative, which is the finding worth publishing."
            )

    # ------------------------------------------------------------------- mechanisms
    mechs: dict[str, int] = {}
    for s in ok:
        if s.limiting:
            mechs[s.limiting] = mechs.get(s.limiting, 0) + 1
    if mechs:
        out.append("")
        out.append("-" * 76)
        out.append("  LIMITING MECHANISM ACROSS THE SET")
        out.append("-" * 76)
        for name, count in sorted(mechs.items(), key=lambda kv: -kv[1]):
            out.append(f"  {name:<36} {count:>4}  ({100 * count / len(ok):.0f}%)")

    out.append("")
    out.append("=" * 76)
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("csv", type=Path, nargs="?",
                    default=Path(__file__).parent / "circuits.csv")
    ap.add_argument("--samples", type=int, default=600)
    ap.add_argument("--out", type=Path, help="also write per-circuit results here")
    args = ap.parse_args(argv)

    if not args.csv.exists():
        print(f"error: {args.csv} not found", file=sys.stderr)
        return 1

    records = load(args.csv)
    if not records:
        print(
            f"No circuits in {args.csv}.\n\n"
            f"This file ships empty on purpose -- the published circuit data is not\n"
            f"redistributed here. See {args.csv.parent / 'README.md'} for exactly where\n"
            f"to get each dataset, or run against synthetic_smoke.csv to exercise the\n"
            f"harness (those numbers are fabricated and must not be reported).",
            file=sys.stderr,
        )
        return 2

    scored = [score_one(r, args.samples) for r in records]
    synthetic = "synthetic" in args.csv.name.lower() or any(
        "synthetic" in s.record.source.lower() for s in scored
    )
    print(report(scored, warn_synthetic=synthetic))

    if args.out:
        with args.out.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow([
                "circuit_id", "source", "outcome", "n_terms", "predicted_bp",
                "reported_bp", "feasible", "predicted_fold_change",
                "reported_fold_change", "limiting", "error",
            ])
            for s in scored:
                w.writerow([
                    s.record.circuit_id, s.record.source, s.record.outcome,
                    s.n_terms, s.predicted_bp, s.record.reported_construct_bp or "",
                    s.feasible,
                    f"{s.predicted_fold_change:.3f}" if s.predicted_fold_change else "",
                    s.record.reported_fold_change or "", s.limiting or "",
                    s.error or "",
                ])
        print(f"\nwrote {args.out}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
