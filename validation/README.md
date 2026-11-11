# Retrospective validation: where to get the data

**No published circuit data is redistributed here.** `circuits.csv` ships as a header
only. This file tells you exactly where to get each dataset and how to get it into the
schema.

## Why this matters more than the compiler

Writing a circuit compiler is a weekend of work. What nobody has published is an
evaluation of whether model-based mammalian circuit design predicts which published
circuits actually worked. Several hundred circuits with known outcomes already exist in
the literature, paid for by other people's pipettes, and using them needs no wet lab.

The result is publishable either way:

- **If predicted dynamic range ranks working circuits above broken ones** (AUC well
  above 0.5), you have the first validated design tool for mammalian recombinase logic.
- **If it does not**, you have a quantified negative result on mammalian composability,
  with a free public benchmark attached. Given that a 2026 study found a plain neural
  network beat a biophysically consistent mechanistic model by r² 0.88 to 0.37 on
  held-out circuit designs, this is the outcome to expect and to be prepared to defend.

Either way the benchmark itself is the durable contribution. Models get reimplemented;
whoever defines the yardstick gets cited by everyone who later competes on it.

**Pre-register the scoring rule before you look at the outcomes.** Decide now what AUC
counts as success, and write it down. Otherwise this becomes a tuning exercise.

## Priority order

Work down this list. Stop when you have enough, but do the first one regardless because
it is the only one that validates the size arithmetic.

### 1. BLADE — the core test set

Weinberg BH, Pham NTH, Caraballo LD, Lozanoski T, Engel A, Bhatia S, Wong WW.
"Large-scale design of robust genetic circuits with multiple inputs and outputs for
mammalian cells." *Nature Biotechnology* 35(5):453–462 (2017).
doi [10.1038/nbt.3805](https://doi.org/10.1038/nbt.3805)

- **113 circuits**, HEK293 and Jurkat T cells, single transcriptional layer via
  recombinase excision. Includes all 16 two-input Boolean functions.
- What you need: per-circuit truth table, measured output per input state, construct
  identity. In the supplementary tables and source data.
- Plasmids are on Addgene under the Wong lab, which is the most reliable way to get
  **exact construct lengths** — and those are what validate section 1 of the benchmark.
- This is the single most important dataset here. The authors stated design software was
  unavailable for mammalian cells and built all 16 gates by hand.

### 2. BLADE mechanistic model — 256 predicted circuits

"A mechanistic model of the BLADE platform predicts performance characteristics of 256
different synthetic DNA recombination circuits." *PLoS Computational Biology*.
doi [10.1371/journal.pcbi.1007849](https://doi.org/10.1371/journal.pcbi.1007849)

- Gives you a **published baseline to compare against**, which is worth more than more
  data. If your compiler cannot beat or match this model on the BLADE circuits, say so.
- Its headline finding — performance degrades as recombinase expression falls — is a
  direct prediction you can test against your own excision-efficiency sensitivity.

### 3. CREATE — the BLADE successor

Letendre JH, Weinberg BH, et al., Wong WW. "Programmable Mixed-Signal Biocomputers in
Mammalian Cells." bioRxiv
[10.1101/2022.06.07.495130](https://doi.org/10.1101/2022.06.07.495130)

- **189 two-input three-output circuits**, plus digital-to-analog converters and a 2-bit
  multiplier. Recombinase logic layer linked by programmable transcription factors.
- Multi-output, so it exercises the compiler's output-grouping and 2A sharing.
- I could not confirm a journal version. Check for one before citing.

### 4. CLASSIC — the large modern set

O'Connell R, Rai K, et al., Bashor CJ. "Ultra-high-throughput mapping of genetic design
space." *Nature* 650:1035–1044 (2026).
doi [10.1038/s41586-025-09933-9](https://doi.org/10.1038/s41586-025-09933-9)

- Roughly **10⁵ circuit designs** measured in HEK293T at single copy via a Bxb1 landing
  pad, with nanopore long reads linking full circuit sequence to barcode.
- Transcription-factor circuits rather than recombinase logic, so it does not test this
  compiler directly. It tests the *premise*: it reports AND/OR-like behaviour as KL
  divergence from idealised gates, which you can use as a soft outcome label.
- **Access is the open question.** The preprint said data was available on request with
  accessions to follow at publication. Code is at
  [cbashorlab/WIMPY](https://github.com/cbashorlab/WIMPY) and
  [cbashorlab/CLASSIC](https://github.com/cbashorlab/CLASSIC), whose README offers
  additional datasets on request.
- **Email Bashor before you commit to this project.** One email decides whether your
  modelling ceiling is COMET's 1.5 MB or something a thousand times larger. It is the
  cheapest and highest-value de-risking step available.

### 5. COMET — part parameters, not circuits

Donahue PS, Draut JW, Muldoon JJ, Edelstein HI, Bagheri N, Leonard JN. *Nature
Communications* 11:779 (2020).
doi [10.1038/s41467-019-14147-5](https://doi.org/10.1038/s41467-019-14147-5) ·
[leonardlab/COMET](https://github.com/leonardlab/COMET) ·
Zenodo [10.5281/zenodo.3533966](https://doi.org/10.5281/zenodo.3533966)

- 44 zinc-finger activators, 12 inhibitors, 83 cognate promoters, measured in MEFL with
  fitted three-parameter response functions. Source Data is about 1.5 MB.
- Use it to **replace this library's estimated parameters with measured ones** via
  `Parts.override_*`. That is the single biggest available improvement to the model.
- Read its limitations honestly: the authors state they cannot predict a new activator's
  expression magnitude or its orthogonality in advance, and report one off-target pair at
  75% crosstalk with no sequence explanation.

### 6. GCAD — a fitted mammalian part library, MIT licensed

Dreyer KS, et al., Leonard JN, Mangan NM. "GCAD: A Computational Framework for Mammalian
Genetic Program Computer-Aided Design." *ACS Synthetic Biology* 15(3):1033–1052 (2026).
doi [10.1021/acssynbio.5c00670](https://doi.org/10.1021/acssynbio.5c00670) ·
[Mangan-Group/GraphGA](https://github.com/Mangan-Group/GraphGA)

- Ships `parts.pkl` / `promo.pkl` with ODE parameters already fitted for 12 activators,
  12 repressors and 156 promoters. **The most directly reusable artifact in this list.**
- Also your closest prior art, so read its validation section carefully: its predictions
  "did not align well" with its own measurements until recalibrated after the fact, and
  it is explicit that it does not do Boolean specification, recombinases, or delivery
  constraints. Those three gaps are what this compiler fills.

### 7. Therapeutic-context circuits

- Li HS, Israni DV, et al., Khalil AS. "Multidimensional control of therapeutic human
  cell function with synthetic gene circuits." *Science* 378(6625):1227–1234 (2022).
  doi [10.1126/science.ade0156](https://doi.org/10.1126/science.ade0156).
  Data on Zenodo [7378325](https://zenodo.org/records/7378325) and Dryad
  `10.5061/dryad.s7h44j19v`, about 1.6 GB, mostly `.fcs` with per-figure spreadsheets.
  Primary human T cells. This is functional experiment data, not a part-characterisation
  matrix, so expect to re-derive anything quantitative from flow files.
- Senti Biosciences, *Cell Systems* (April 2026): **more than 60 NOT-gated CAR circuit
  designs** in T and NK cells. A clean held-out discrimination set in exactly the
  therapeutic context this tool targets.
- Roas K, et al., Nissim L. *Nature Communications* 17:7827 (2026).
  doi [10.1038/s41467-026-74408-y](https://doi.org/10.1038/s41467-026-74408-y).
  Plasmids Addgene 251264–251338, RNA-seq GEO `GSE316863`. Characterisation is binary
  ON/OFF, and the paper states it reports no original code.

## Schema

`circuits.csv` columns:

| column | required | meaning |
|---|---|---|
| `circuit_id` | yes | unique label; rows starting `#` are skipped |
| `source` | | short citation key, e.g. `BLADE2017` |
| `truth_table` | yes | output column MSB-first (`0110`), explicit rows (`00:0,01:1,...`), or a minterm list (`1,2`) |
| `inputs` | | `;`-separated input names, e.g. `A;B` |
| `payload_bp` | | payload coding length; defaults to 720 (EGFP) |
| `recombinases` | | `;`-separated, one per input, e.g. `Cre;FlpO` |
| `outcome` | | `worked`, `failed`, `partial`, `unknown` |
| `reported_fold_change` | | measured ON/OFF from the paper |
| `reported_construct_bp` | | construct length, ideally from the Addgene map |
| `vector` | | `plasmid`, `aav_single`, `lentivirus`, `landing_pad`, ... |
| `notes` | | free text |

Define `outcome` once, up front, and apply it mechanically. A defensible rule: `worked`
if the paper's own measured ON/OFF exceeds a stated threshold (say 10×) in the intended
states, `failed` if it does not, `partial` if the paper reports it as functional but
off-specification. **Do not** assign outcomes by reading the narrative.

## Running it

```bash
python validation/benchmark.py validation/circuits.csv --out results.csv
```

To exercise the harness before you have real data:

```bash
python validation/benchmark.py validation/synthetic_smoke.csv
```

`synthetic_smoke.csv` contains **fabricated rows**. It exists only to prove the script
runs. Every number it produces is meaningless and the harness prints a warning saying so.
