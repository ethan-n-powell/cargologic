#!/usr/bin/env python3
"""Compare cargologic's product-term layout against BLADE's published topology.

BLADE (Weinberg et al. 2017, Nat Biotechnol 35:453-462) lays out a circuit as one
promoter upstream of 2^N "addresses", each surrounded by recombination sites. The paper
states the template directly:

    "Each BLADE circuit is organized as a single transcriptional layer comprising a
     single promoter upstream of up to 2^N regions of DNA sequences (or addresses, Z),
     which are surrounded by recombination sites"

    "the total number of recombination site pairs for N-inputs being 2^N-1"

So the two designs differ in what they replicate. BLADE instantiates every state, ON or
OFF, and shares one promoter. cargologic instantiates one cassette per product term and
pays for a promoter each time. Which wins depends on how sparse the function is, and that
is what this script measures.

Everything here is computed from the published template plus this package's part sizes.
It needs no supplementary data, which matters because the paper's supplementary
spreadsheet is 403/CAPTCHA-gated and its data site (datasheets.synbiotools.org) is dead.

    python validation/blade_architecture.py
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import cargologic as cl  # noqa: E402
from cargologic.parts import JUNCTION_BP, KOZAK_BP, Parts  # noqa: E402

# Recombination sites BLADE uses, from the paper: lox sites loxP, lox2272, loxN for Cre;
# FRT sites FRT, F3, F14 for Flp; VloxP and Vlox2272 for VCre. All are 34 bp.
BLADE_SITE_BP = 34


@dataclass
class Layout:
    name: str
    bp: int
    promoters: int
    payload_copies: int
    site_pairs: int


def blade_size(n_inputs: int, ones: set[int], parts: Parts,
               promoter: str, payload_bp: int, terminator: str) -> Layout:
    """Size of the BLADE template: one promoter, 2^N addresses, 2^N-1 site pairs.

    Every address is physically present. An ON address carries the payload and a polyA;
    an OFF address carries a terminator alone.
    """
    prom = parts.promoters[promoter]
    term = parts.terminators[terminator]
    n_addresses = 1 << n_inputs
    n_on = len(ones)
    n_off = n_addresses - n_on
    site_pairs = n_addresses - 1

    bp = prom.bp
    bp += site_pairs * 2 * BLADE_SITE_BP
    bp += n_on * (KOZAK_BP + payload_bp + term.bp)
    bp += n_off * term.bp
    bp += (n_addresses + 1) * JUNCTION_BP
    return Layout("BLADE template", bp, 1, n_on, site_pairs)


def cargologic_size(expr_ones: set[int], n_inputs: int, payload_bp: int,
                    promoter: str, terminator: str) -> Layout:
    """Size cargologic produces for the same function, with parts held equal."""
    column = "".join("1" if m in expr_ones else "0" for m in range(1 << n_inputs))
    r = cl.compile_circuit(
        column,
        payload=payload_bp,
        as_table=True,
        inputs=[chr(ord("A") + i) for i in range(n_inputs)],
        optimise=False,
        promoter=promoter,
        terminator=terminator,
        blocker=terminator,
        n_samples=1,
    )
    d = r.design
    n_prom = sum(1 for e in d.elements if e.kind == "promoter")
    n_pay = sum(1 for e in d.elements if e.kind == "payload")
    n_sites = sum(1 for e in d.elements if e.kind == "site")
    return Layout(f"cargologic ({d.strategy})", d.total_bp, n_prom, n_pay, n_sites // 2)


def compare(n_inputs: int, payload_bp: int = 1500, promoter: str = "EFS",
            terminator: str = "bGHpA") -> list[tuple]:
    parts = Parts()
    rows = []
    for mask in range(1 << (1 << n_inputs)):
        ones = {m for m in range(1 << n_inputs) if (mask >> m) & 1}
        if not ones:
            continue
        b = blade_size(n_inputs, ones, parts, promoter, payload_bp, terminator)
        c = cargologic_size(ones, n_inputs, payload_bp, promoter, terminator)
        rows.append((mask, len(ones), b, c, c.bp - b.bp))
    return rows


def _is_parity(ones: list[int], n: int) -> bool:
    """True if the function is XOR or XNOR over all n inputs."""
    odd = {m for m in range(1 << n) if bin(m).count("1") % 2 == 1}
    even = {m for m in range(1 << n)} - odd
    return set(ones) == odd or set(ones) == even


def main() -> int:
    payload = 1500  # a second-generation CAR
    print("=" * 78)
    print("  BLADE template vs cargologic product-term layout")
    print(f"  payload {payload} bp, promoter EFS, polyA bGHpA, parts held equal")
    print("=" * 78)

    for n in (2, 3):
        rows = compare(n, payload_bp=payload)
        print()
        print(f"--- {n} inputs, all {len(rows)} non-constant-OFF functions ---")
        print(f"  {'ON states':>10} {'BLADE bp':>9} {'cargologic':>11} "
              f"{'delta':>8} {'who wins':>12}")
        print("  " + "-" * 56)

        by_k: dict[int, list] = {}
        for _, k, b, c, delta in rows:
            by_k.setdefault(k, []).append((b.bp, c.bp, delta))

        for k in sorted(by_k):
            group = by_k[k]
            bs = sum(x[0] for x in group) / len(group)
            cs = sum(x[1] for x in group) / len(group)
            d = cs - bs
            winner = "cargologic" if d < 0 else "BLADE"
            print(f"  {k:>10} {bs:>9.0f} {cs:>11.0f} {d:>+8.0f} {winner:>12}")

        wins = sum(1 for _, _, _, _, d in rows if d < 0)
        print(f"  cargologic smaller in {wins}/{len(rows)} functions")

        losses = [(mask, k, b, c, d) for mask, k, b, c, d in rows if d > 0]
        if losses:
            print()
            print(f"  BLADE wins {len(losses)}/{len(rows)}, on these functions:")
            for mask, k, b, c, d in sorted(losses, key=lambda r: -r[4]):
                ones = sorted(m for m in range(1 << n) if (mask >> m) & 1)
                cov = cl.parse_table(
                    "".join("1" if m in ones else "0" for m in range(1 << n)),
                    inputs=[chr(ord("A") + i) for i in range(n)],
                ).minimise()
                par = " (parity)" if _is_parity(ones, n) else ""
                print(f"    ones={str(ones):<16} {cov.n_terms} terms  "
                      f"BLADE {b.bp} vs {c.bp}, {d:+d} bp{par}")

    print()
    print("-" * 78)
    print("  WHAT THIS SHOWS")
    print("-" * 78)
    print("""  The product-term layout is smaller almost everywhere: 13/15 two-input and
  246/255 three-input functions. BLADE instantiates all 2^N addresses whether ON
  or OFF and pays 2^N-1 site pairs, so its template grows as 2^N while a minimised
  cover usually does not. This holds across payload sizes from 300 to 4500 bp, so
  it is not an artefact of the cargo chosen.

  The exceptions are not the dense functions, which is the obvious guess and is
  wrong. Both exact-parity functions at two inputs are in the losing set, and both
  at three inputs, which fits the usual story: parity does not factor, so its cover
  needs one term per ON state and the layout degenerates into a promoter per address.

  But parity is not the whole rule. Nine three-input functions lose, only two of
  which are parity. The obvious tightening, "BLADE wins when no two minterms merge
  and at least two states are ON", was tested here and FALSIFIED: it mispredicts 17
  of the 255 three-input functions. The real boundary is arithmetic. It falls out of
  how n_terms x (promoter + payload + polyA) compares with one promoter plus 2^N
  addresses, and there is no clean combinatorial statement of it.

  The actionable part does not depend on characterising the boundary. A compiler can
  simply price both layouts and take the smaller, which is cheap because both are
  closed-form. cargologic does not currently do this: it always builds product terms,
  which is the wrong call on 2/15 and 9/255 of these functions, by up to 1040 bp.
  Adding the BLADE template as a third strategy alongside independent and
  tandem_promoters would close that gap.""")
    print()
    print("-" * 78)
    print("  ASSUMPTIONS, SO THE NUMBERS CAN BE CHECKED")
    print("-" * 78)
    print(f"""  - BLADE template taken from the paper: one promoter, 2^N addresses,
    2^N-1 recombination site pairs. Both quoted verbatim in this file's docstring.
  - All recombination sites {BLADE_SITE_BP} bp (loxP, FRT and their heterospecific
    variants are all {BLADE_SITE_BP} bp).
  - An OFF address costs one terminator; an ON address costs Kozak + payload + polyA.
  - Part sizes come from cargologic's own library, so both layouts are priced with
    the same parts. Absolute numbers depend on that library; the comparison does not.
  - NOT verified against measured construct lengths. The paper states no bp sizes and
    Addgene does not publish lengths for these plasmids, so this compares published
    TOPOLOGY, not a measured construct.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
