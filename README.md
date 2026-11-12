# cargologic

**Will your genetic circuit fit in the vector you plan to ship it in?**

You give it a truth table and a delivery format. It gives you a DNA architecture, an
exact length, a packaging verdict, and a prediction of what will break.

```bash
pip install -e .
cargologic "A and not B" --payload CAR_2g --vector aav_single
```

```
================================================================================
  FITS  -  AAV (single vector)
  1862 bp circuit, 1 cassette(s)
  +2578 bp against a 4440 bp budget (42% used)
  predicted ON/OFF ratio ~6x (80% interval 4-9x)
================================================================================
```

## The idea

A recombinase logic circuit realises one product term per physical cassette. So the thing
that decides whether your circuit fits in an AAV is not how complicated the logic looks —
it is **how many product terms survive Boolean minimisation.**

That makes Quine–McCluskey a packaging optimisation, which is a slightly strange sentence
but is exactly what is going on. Shaving a promoter saves a few hundred base pairs.
Eliminating a product term saves a kilobase and change.

The compiler puts that relationship at the centre, then adds the constraint that nothing
else in the field models: **the vector you have to ship it in.**

```
cassette = [promoter] <a> [blocker polyA] <a>  <b| [payload] [polyA] |b>

  positive literal "a must be TRUE"  -> a polyA blocker flanked by a-sites.
                                        Recombinase a excises it, unblocking the payload.
                                        Several in series give AND.
  negative literal "b must be FALSE" -> b-sites flanking the payload itself.
                                        If b fires, the payload is physically removed.
  OR across terms                    -> one cassette per term. This is the expensive part.
```

## What is actually new here

The gap is specific and nine years old.

- **BLADE** (Weinberg et al., *Nat Biotechnol* 2017) built a mammalian recombinase logic
  platform and demonstrated all 16 two-input Boolean functions. The authors noted design
  software was unavailable for mammalian cells and built every gate by hand. BU announced
  software was coming. It never shipped.
- **Cello** does truth-table-to-DNA compilation, for bacteria. **CALIN** does it for
  recombinase logic, also for bacteria.
- **GCAD** (*ACS Synth Biol*, Feb 2026) is genuine mammalian circuit CAD — and it
  *deliberately rejects* Boolean specification, has no recombinases, caps circuits at two
  parts, and does not model delivery at all.
- **Asimov's Kernel** does expression vectors and codon optimisation, not circuit logic.

So: no Boolean compiler for mammalian recombinase logic, and **nothing anywhere that puts
AAV packaging limits inside the design loop.** That second gap is the useful one, because
it is the only part of this problem with no biology risk.

## Two kinds of output, kept apart on purpose

This is the design principle of the whole tool.

**Sizes are arithmetic.** AAV holds 4.7 kb; the ITRs take 260 bp; your cassettes take what
they take. A verdict of "this does not fit" cannot be overturned by anything we later
learn about mammalian gene expression. The report says *trust this section*.

**Expression predictions are a model, and models of mammalian circuits have a poor recent
record.** GCAD's ODE predictions "did not align well" with its own measurements. CLASSIC
(*Nature*, Jan 2026) found a plain neural network beat a biophysically consistent
mechanistic model on held-out designs by r² 0.88 to 0.37 — plain linear regression also
beat the mechanistic model. So the behaviour model here reports an interval, not a number,
and is explicitly labelled as being for ranking.

### The part that survives being wrong

Rather than lean on absolute levels, the tool runs a sensitivity analysis and tells you
**which mechanism is limiting you** — because that ranking depends on the *ordering* of
effects, not their magnitudes, so it survives the parameters being off by a factor.

For `A AND NOT B`:

```
    what is limiting state 11 (share of the excess each mechanism accounts for)
      incomplete excision                100.0%  ############################
      terminator read-through              0.6%
      spontaneous recombinase activity     0.0%

      -> Use a more efficient recombinase, raise its dose or expression time,
      or switch to an irreversible serine integrase. A tighter polyA will not help.
```

That advice is checkable, and it checks out. Sweeping the terminator across all five
options changes the dynamic range by nothing (5.7× throughout) while costing up to 862 bp.
Improving excision efficiency from 70–95% to 90–97% nearly doubles it:

| FlpO excision efficiency | predicted ON/OFF | limiting mechanism |
|---|---|---|
| 0.70–0.95 | 5.7× | incomplete excision |
| 0.90–0.97 | 10.7× | terminator read-through |
| 0.98–0.999 | 10.8× | terminator read-through |

This is the characteristic weakness of NOT logic built on excision, and it is a
population-heterogeneity problem rather than a transcriptional one: if excision is 85%
efficient, 15% of cells keep the payload and read as ON no matter how good your
terminators are. The tool says so instead of blaming leak.

## Usage

```bash
# a named gate, and the self-check across all 16 two-input functions
cargologic AND3 --payload 1500 --vector aav_dual
cargologic --all-two-input --payload CAR_2g

# truth tables: output column, explicit rows, or a minterm list
cargologic "0110" --table --payload EGFP
cargologic "00:0,01:1,10:1,11:0" --table --payload EGFP
cargologic "1,2" --table --payload EGFP

# machine-readable, and a GenBank feature table
cargologic "A & B & !C" --payload CAR_2g --json
cargologic "A & B & !C" --payload CAR_2g --genbank design.gb

# what happens if the circuit must carry its own recombinases
cargologic "A & !B" --payload CAR_2g --include-recombinases

cargologic --list-parts
```

```python
import cargologic as cl

r = cl.compile_circuit("A & B & !C", payload="CAR_2g", vector="aav_single")
print(r.feasible, r.circuit_bp)
print(r.report())

# multi-output: a CAR and a cytokine on different logic
r = cl.compile_circuit(
    {"car": "A & !C", "il12": "A & B & !C"},
    payload={"car": "CAR_2g", "il12": "IL12"},
)

# swap in your own measured parameters -- the shipped ones are estimates
p = cl.Parts()
p.override_recombinase("Bxb1", excision_efficiency=(0.96, 0.99), leak=(1e-4, 1e-3))
r = cl.compile_circuit("A & !B", payload=1500, parts=p)
```

## Optimisations it finds for you

| Move | When | Saving |
|---|---|---|
| Boolean minimisation | always | ~1–2 kb per term eliminated |
| Shared payload behind tandem gated promoters | when all terms share the same negative literals | `(k−1) × (payload + polyA)` |
| BLADE template instead of product terms | when the cover doesn't factor, parity being the clear case | up to ~1 kb, and usually better dynamic range |
| 2A-linked cassette merging | two outputs firing in identical states | one promoter + polyA per merged output |
| Part downsizing search | near the packaging limit | up to ~1.5 kb, traded against leak |

OR and XOR both minimise to two terms, but OR's terms share negative literals so one
payload copy serves both: **2473 bp versus 3544 bp** for the same payload.

### Three layouts, priced against each other

The compiler builds all applicable layouts and takes the best, because which one wins is
not predictable from the logic alone:

- **independent** — one transcription unit per product term. Always valid.
- **tandem_promoters** — gated promoters in tandem sharing one payload. Needs every term
  to have the same negative literals.
- **blade** — the template from Weinberg et al. 2017: one promoter upstream of 2^N
  addresses with 2^N−1 recombination site pairs. Every address exists whether ON or OFF.

`validation/blade_architecture.py` prices the first two against the third across every
function of two and three inputs. Product terms win 13/15 and 245/255, stable across
payloads from 300 to 4500 bp — the BLADE template grows as 2^N while a minimised cover
usually doesn't. The exceptions include both parity functions at each input count, where
the cover can't factor and product terms degenerate into a promoter per address.

Building that comparison is also what caught a real bug: recombinase site pairs were
being counted at half their true length, because the part library mixed per-site and
per-pair figures. Two independent implementations disagreeing is how it surfaced.

Near the limit there is no single best design, so you get the trade-off instead of a
fabricated optimum. For a pure `A AND B`, 117 extra base pairs buys **3.4× the dynamic
range**, by switching layout rather than by buying a tighter polyA:

```
SIZE vs DYNAMIC RANGE (PARETO FRONT)  --  A & B
     bp   ON/OFF   promoter  blocker   terminator   layout
   1917     9.5x   minCMV    SPA       SPA          independent
   2034    32.5x   minCMV    SV40pA    SPA          blade
```

### A design rule falls out of this

Run the same search on `A AND B AND NOT C` and the front **collapses to a single point**.
Nothing you can buy with base pairs improves it.

That is the answer, not a bug: what you can fix with budget depends on **what is actually
holding the OFF state down**, and that differs by gate.

| gate | chosen layout | limiting mechanism | can budget fix it? |
|---|---|---|---|
| `A & B` | blade | spontaneous recombinase activity | yes — 117 bp buys 3.4× |
| `A & B & C` | independent | terminator read-through | yes — tighter polyA works |
| `A \| B` | tandem_promoters | spontaneous recombinase activity | partly |
| `A & B & !C` | independent | incomplete excision | **no — front is one point** |

A gate held by a blocker has a terminator problem, and a tighter polyA fixes it. A gate
with a NOT puts its OFF state behind a physical excision event, so the floor is set by
the fraction of cells that failed to excise, and no terminator reaches that. Spend it on
a better recombinase instead.

Note what the BLADE layout does for `A & B`: it has no blocker at all, so there is no
read-through path to leak, and its OFF states are set purely by recombinase behaviour.
That is why it wins the dynamic range despite being slightly larger.

## Validation

`validation/` holds a retrospective benchmark harness. **It ships without data**, and the
README there says exactly where to get each dataset.

This is the part of the project that matters most. Writing a compiler is a weekend.
Nobody has published an evaluation of whether model-based mammalian circuit design
predicts which published circuits worked — and roughly **300–550 circuits with known
outcomes already exist** across BLADE (113), CREATE (189) and the BLADE mechanistic model
(256). Using them needs no wet lab.

It scores four things: size accuracy (arithmetic, should be near-exact), feasibility
recall, discrimination (AUC on working vs failed), and calibration (Spearman on
fold-change).

The result is publishable either way. If predicted dynamic range ranks working circuits
above broken ones, that is the first validated design tool for mammalian recombinase
logic. If it does not, that is a quantified negative result on mammalian composability
with a public benchmark attached — and given the CLASSIC result above, that is the outcome
to be prepared for. **Pre-register the scoring rule before looking at outcomes.**

```bash
python validation/benchmark.py validation/circuits.csv --out results.csv
python validation/benchmark.py validation/synthetic_smoke.csv   # fabricated; smoke test only
```

## Limitations, stated plainly

- **Everything quantitative in mammalian circuit literature is HEK293 transient
  transfection.** If your target is primary T cells, NK cells or iPSCs, treat every
  efficiency and leak value here as a guess, and expect promoter silencing this tool does
  not model. v1 is honestly a HEK293 tool.
- **Not modelled:** transcriptional squelching, competition for shared transcriptional and
  degradation resources, retroactivity, chromosomal position effects, transgene silencing.
  Each has broken a real published design.
- **Part parameters are literature-informed estimates with ranges, not measurements.**
  Override them. `Parts.override_*` exists for exactly this.
- **No sequence is emitted.** The GenBank output has real coordinates and N-masked
  sequence; take the part sequences from your own plasmids or Addgene. Inventing sequence
  would look more impressive and be less useful.
- **Orthogonality is assumed from published pairings.** Verify cross-reactivity before
  committing to a multi-input build.
- Single transcriptional layer only, which is what makes the near-digital treatment
  defensible. Multi-layer cascades need the continuous modelling this tool avoids.

## Development

```bash
pip install -e ".[dev]"
pytest -q          # 88 tests, no dependencies beyond pytest
```

The test suite checks that minimisation never alters a function across **all 2, 16 and 256
functions of 1, 2 and 3 inputs** plus sampled 4-input functions, that the feature table
tiles the design with no gaps or overlap, that a term costs more than any part swap, and
that the mechanism diagnosis is directionally correct.

MIT licensed.
