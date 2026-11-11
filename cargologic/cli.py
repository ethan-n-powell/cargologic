"""Command line interface.

    cargologic "A and not B" --payload CAR_2g --vector aav_single
    cargologic AND3 --payload 1500 --vector aav_dual --genbank out.gb
    cargologic "0110" --table --payload EGFP
    cargologic --list-parts
"""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .compiler import compile_circuit
from .delivery import vector_catalogue
from .logic import LogicError, NAMED_GATES, all_functions
from .parts import PAYLOADS, PROMOTERS, RECOMBINASES, TERMINATORS, Parts
from .synthesis import SynthesisError


def _list_parts() -> str:
    out = ["RECOMBINASES"]
    for k, v in RECOMBINASES.items():
        out.append(
            f"  {k:<8} sites {v.site_name:<10} {v.site_bp:>3} bp/pair  "
            f"CDS {v.cds_bp:>4} bp  eff {v.excision_efficiency[0]:.2f}-"
            f"{v.excision_efficiency[1]:.2f}  leak {v.leak[0]:.4f}-{v.leak[1]:.4f}"
        )
    out.append("\nPROMOTERS")
    for k, v in sorted(PROMOTERS.items(), key=lambda kv: kv[1].bp):
        out.append(f"  {k:<8} {v.bp:>5} bp  relative strength {v.strength:.2f}")
    out.append("\nTERMINATORS / BLOCKERS")
    for k, v in sorted(TERMINATORS.items(), key=lambda kv: kv[1].bp):
        out.append(
            f"  {k:<8} {v.bp:>5} bp  read-through "
            f"{v.readthrough[0]:.3f}-{v.readthrough[1]:.3f}"
        )
    out.append("\nPAYLOADS")
    for k, v in sorted(PAYLOADS.items(), key=lambda kv: kv[1].bp):
        out.append(f"  {k:<12} {v.bp:>5} bp  {v.note}")
    out.append("\nVECTORS")
    for k, v in vector_catalogue().items():
        budget = v.budget()
        out.append(f"  {k:<14} budget {budget if budget else 'unbounded'}")
    out.append("\nNAMED GATES")
    out.append("  " + ", ".join(sorted(NAMED_GATES)))
    return "\n".join(out)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="cargologic",
        description="Compile Boolean logic into a recombinase circuit that fits your vector.",
        epilog='examples:\n'
        '  cargologic "A and not B" --payload CAR_2g\n'
        '  cargologic MAJORITY3 --payload 1500 --vector aav_dual\n'
        '  cargologic "1,2" --table --payload EGFP --json\n',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("logic", nargs="?", help="Boolean expression, named gate, or truth table")
    p.add_argument("--payload", default="EGFP",
                   help="library name, bp count, or name:bp (default EGFP)")
    p.add_argument("--vector", default="aav_single",
                   help="target delivery format (default aav_single)")
    p.add_argument("--table", action="store_true",
                   help="parse LOGIC as a truth table rather than an expression")
    p.add_argument("--inputs", help="comma-separated input names")
    p.add_argument("--recombinases", help="comma-separated recombinases, one per input")
    p.add_argument("--include-recombinases", action="store_true",
                   help="count recombinase CDSs against the packaging budget")
    p.add_argument("--promoter", help="force a promoter instead of searching")
    p.add_argument("--blocker", help="force the blocker polyA")
    p.add_argument("--terminator", help="force the trailing polyA")
    p.add_argument("--no-optimise", action="store_true", help="skip the part search")
    p.add_argument("--min-dynamic-range", type=float, default=10.0,
                   help="dynamic-range floor the search must clear (default 10)")
    p.add_argument("--samples", type=int, default=2000,
                   help="Monte-Carlo samples for the behaviour model (default 2000)")
    p.add_argument("--genbank", metavar="PATH", help="write a GenBank feature table")
    p.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    p.add_argument("--width", type=int, default=80, help="report width")
    p.add_argument("--list-parts", action="store_true", help="show the part library and exit")
    p.add_argument("--all-two-input", action="store_true",
                   help="compile all 16 two-input functions as a self-check")
    p.add_argument("--version", action="version", version=f"cargologic {__version__}")
    return p


def _json_payload(result) -> dict:
    return {
        "inputs": result.spec.inputs,
        "outputs": {
            f.name: {
                "expression": f.source,
                "ones": sorted(f.ones),
                "cover": result.covers[f.name].to_str(result.spec.inputs),
                "n_terms": result.covers[f.name].n_terms,
            }
            for f in result.spec.functions
        },
        "recombinases": [r.name for r in result.design.recombinases],
        "strategy": result.design.strategy,
        "circuit_bp": result.design.total_bp,
        "n_cassettes": result.design.n_terms,
        "parts": {
            "promoter": result.choice.promoter,
            "blocker": result.choice.blocker,
            "terminator": result.choice.terminator,
        },
        "target": {
            "vector": result.target_verdict.vector,
            "fits": result.target_verdict.fits,
            "budget_bp": result.target_verdict.budget_bp,
            "headroom_bp": result.target_verdict.headroom_bp,
            "severity": result.target_verdict.severity,
        },
        "delivery": [
            {
                "vector": v.vector,
                "fits": v.fits,
                "budget_bp": v.budget_bp,
                "headroom_bp": v.headroom_bp,
                "relative_titre": v.relative_titre,
                "severity": v.severity,
            }
            for v in result.all_verdicts
        ],
        "behaviour": {
            name: {
                "dynamic_range_median": (
                    None if rep.dynamic_range_median == float("inf")
                    else rep.dynamic_range_median
                ),
                "dynamic_range_p10": (
                    None if rep.dynamic_range_lo == float("inf") else rep.dynamic_range_lo
                ),
                "dynamic_range_p90": (
                    None if rep.dynamic_range_hi == float("inf") else rep.dynamic_range_hi
                ),
                "weakest_on_state": rep.weakest_on_state,
                "leakiest_off_state": rep.worst_off_state,
                "warnings": rep.warnings,
                "states": [
                    {
                        "bits": s.bits,
                        "intended": s.intended,
                        "level_median": s.level_median,
                        "level_p10": s.level_lo,
                        "level_p90": s.level_hi,
                    }
                    for s in rep.states
                ],
            }
            for name, rep in result.leakage.items()
        },
        "features": result.design.feature_table(),
    }


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.list_parts:
        print(_list_parts())
        return 0

    if args.all_two_input:
        return _self_check(args)

    if not args.logic:
        build_parser().print_help()
        return 2

    inputs = [s.strip() for s in args.inputs.split(",")] if args.inputs else None
    recs = (
        [s.strip() for s in args.recombinases.split(",")] if args.recombinases else None
    )

    try:
        result = compile_circuit(
            args.logic,
            payload=args.payload,
            vector=args.vector,
            inputs=inputs,
            recombinases=recs,
            include_recombinases=args.include_recombinases,
            optimise=not args.no_optimise,
            min_dynamic_range=args.min_dynamic_range,
            as_table=args.table,
            promoter=args.promoter,
            blocker=args.blocker,
            terminator=args.terminator,
            n_samples=args.samples,
        )
    except (LogicError, SynthesisError, KeyError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if args.genbank:
        with open(args.genbank, "w") as fh:
            fh.write(result.to_genbank())
        print(f"wrote {args.genbank}", file=sys.stderr)

    if args.json:
        print(json.dumps(_json_payload(result), indent=2))
    else:
        print(result.report(width=args.width))

    return 0 if result.feasible else 3


def _self_check(args) -> int:
    """Compile every two-input Boolean function.

    BLADE demonstrated all 16 experimentally, so a compiler that cannot emit a design
    for each is wrong. This is a correctness check on the compiler, not a claim about
    biological performance.
    """
    print("Compiling all 16 two-input Boolean functions")
    print(f"payload={args.payload}  vector={args.vector}")
    print()
    print(f"{'mask':>4}  {'function':<30} {'terms':>5} {'bp':>6} {'fits':>5} {'ON/OFF':>8}")
    print("-" * 66)
    failures = 0
    for f in all_functions(2):
        # Specify each function as an explicit output column so constant-0 and
        # constant-1 are expressible (a minterm list cannot represent the empty set).
        column = "".join("1" if m in f.ones else "0" for m in range(4))
        try:
            r = compile_circuit(
                column,
                payload=args.payload,
                vector=args.vector,
                as_table=True,
                inputs=f.inputs,
                n_samples=300,
            )
        except Exception as exc:  # pragma: no cover
            print(f"{f.name:>4}  FAILED: {exc}")
            failures += 1
            continue
        cov = next(iter(r.covers.values()))
        rep = next(iter(r.leakage.values()))
        dr = "inf" if rep.dynamic_range_median == float("inf") else f"{rep.dynamic_range_median:.0f}x"
        print(
            f"{f.name:>4}  {cov.to_str(f.inputs):<30} {cov.n_terms:>5} "
            f"{r.circuit_bp:>6} {'y' if r.feasible else 'n':>5} {dr:>8}"
        )
    print()
    print(f"{16 - failures}/16 compiled" + (" -- ALL OK" if not failures else ""))
    return 1 if failures else 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
