# Superfold Architecture

This document is intended to be a quick reference for developers to understand the overall structure of Superfold.

## What It Does

Superfold is a bioinformatics pipeline for RNA secondary structure prediction.
It takes SHAPE-MaP chemical probing data (nucleotide reactivity measurements) as input and uses it to guide thermodynamic folding via the RNAstructure software suite.
The key innovation is a windowed, parallelized approach that makes large RNA folding (>5000 nt) tractable.
It also adds ensemble-level statistics; low Shannon entropy from base-pairing probabilities and low SHAPE reactivity identify high-confidence structures.

**Typical invocation:**
```bash
python Superfold.py myRNA.map --np 8 --foldWindowSize 3000 --partitionWindowSize 1200
```

---

## Repository Layout

```
Superfold/
├── Superfold.py                    # Main pipeline entry point (~1491 lines)
├── RNAtools.py                     # Core data structures: CT, dotPlot (~858 lines)
├── batchSubmit.py                  # Parallel job execution via multiprocessing (~241 lines)
├── pvclient.py                     # PseudoViewer web client for structure rendering (~833 lines)
├── PyCircleCompareSF.py            # Circular secondary structure plots in PostScript (~918 lines)
├── drawArcRibbons_simple.py        # Arc ribbon plots via matplotlib (~451 lines)
├── differenceByWindowSHAPEMAP.py   # Differential SHAPE reactivity analysis (~449 lines)
├── 16S.1m7.map.txt                 # Example SHAPE-MaP data file
├── README                          # Installation and usage documentation
└── COPYING                         # GNU GPL v3 license
```

---

## Data Flow

```
Input: myRNA.map (SHAPE-MaP data file)
  columns: nt_index | SHAPE_reactivity | std_error | nucleotide

  ↓ parseArgs()
  Load map file into shapeMAP object. Parse CLI args, constraints files.

  ↓ genFiles()
  For each window [ntStart, ntEnd]:
    Write windowed .seq, .shape, .con (constraint) files to temp folder.

  ↓ Partition Phase — generateAndRunPartition()
  │
  ├── batchSubmit() → RNAstructure `partition` (parallel, N workers)
  │     Input: .seq + .shape files per window
  │     Output: .pfs (partition function) files
  │
  ├── batchSubmit() → RNAstructure `ProbabilityPlot` (parallel)
  │     Input: .pfs files
  │     Output: .dp (dot plot / base pair probability) files
  │
  └── mainAssemble() → concatonateDP()
        Read all windowed .dp files, trim edges, merge overlapping windows
        by averaging duplicate base pair probabilities.
        Output: merged dotPlot object

  ↓ Fold Phase — generateAndRunFold()
  │
  ├── batchSubmit() → RNAstructure `Fold` (parallel, N workers)
  │     Input: .seq + .shape + .con files per window
  │     Output: .ct (connect table / MFE structure) files
  │
  └── MasterModel_readAndRenumberAll() → MasterModel_findOverlapPairs()
        Read all windowed .ct files, renumber to full-length coordinates,
        find consensus pairs present in majority of overlapping windows.
        Output: consensus CT object

  ↓ Analysis
  │
  ├── calcShannon()       — Shannon entropy per nucleotide from merged dotPlot
  └── mainShannonFunc()  — Identify low-entropy, low-SHAPE (well-structured) regions

  ↓ Visualization
  │
  ├── ensembleRNA_splitPlot() → drawArcRibbons_simple.writeManyPlot()
  │     Arc ribbon plot colored by base pair probability → PDF
  │
  ├── PyCircleCompareSF.makeCircle()
  │     Circular secondary structure plot annotated with SHAPE + Shannon → PostScript
  │
  └── pvclient.python_client()
        Render each low-Shannon region structure via PseudoViewer web service → EPS + XRNA

Output directory: {prefix}_Superfold/
  ├── merged_*.dp           — merged partition function (dot plot)
  ├── merged_*.ct           — consensus MFE structure
  ├── merged_wPK_*.ct       — structure with pseudoknots
  ├── shannonShape_*.pdf    — Shannon entropy + SHAPE reactivity plot
  ├── arcPlot_*.pdf         — arc ribbon secondary structure diagram
  ├── regions/
  │   ├── *.ct              — per-region structures
  │   ├── *.ps              — per-region PostScript renders
  │   └── *.xrna            — per-region XRNA files
  └── log_*.txt             — execution log
```

---

## Module Descriptions

### Superfold.py — Pipeline Orchestration

The main executable. Owns the full pipeline lifecycle: argument parsing, windowed file generation, job dispatch, result assembly, and visualization calls.

**Key components:**

| Function | Role |
|---|---|
| `main()` | Top-level orchestration — calls all phases in order |
| `parseArgs()` | argparse-based CLI; sets window sizes, step sizes, nprocs, SHAPE slope/intercept |
| `runCheck()` | Validates RNAstructure executables (`Fold`, `partition`, `ProbabilityPlot`) are in PATH |
| `shapeMAP` class | Holds parsed SHAPE data: `seq`, `ntNum`, `shape`, `stdErr`, `zfactor` |
| `genFiles()` | Writes windowed `.seq` / `.shape` / `.con` input files for a given `[ntStart, ntEnd]` |
| `generateAndRunPartition()` | Builds partition commands for all windows, dispatches via batchSubmit, assembles merged DP |
| `generateAndRunFold()` | Builds Fold commands for all windows, dispatches via batchSubmit, assembles consensus CT |
| `mainAssemble()` | Reads windowed `.dp` files, trims terminal ends, returns merged dotPlot |
| `concatonateDP()` | Merges duplicate base pairs across windows by averaging `-log10(probability)` |
| `MasterModel_readAndRenumberAll()` | Reads all windowed CT files, renumbers to global coordinates |
| `MasterModel_findOverlapPairs()` | Consensus vote: keeps pairs present in majority of overlapping windows |
| `mainShannonFunc()` | Finds contiguous low-Shannon, low-SHAPE regions for detailed rendering |
| `ensembleRNA_splitPlot()` | Creates arc plot PDF colored by partition function probabilities |
| `progress()` | CLI progress bar |

**Windowing logic:**
- RNA is divided into overlapping windows (`windowSize`, `stepSize`)
- Edge windows at 5' and 3' ends use shorter sizes (+/-100, +/-50 nt offsets) to ensure coverage
- Partition windows default: size=1200, step=100
- Fold windows default: size=3000, step=300

---

### RNAtools.py — Core Data Structures

Defines the two primary data objects used throughout the pipeline.

**`CT` class — secondary structure (connect table format):**
- `num`: list of nucleotide indices (1-based)
- `seq`: list of nucleotide characters
- `ct`: list of pairing partners (0 = unpaired)
- `name`: structure label/title
- Key methods: `readCT()`, `writeCT()`, `pairList()`, `pair2CT()`, `cutCT()`, `contactDistance()`, `extractHelices()`, `extractPK()`, `fillPairs()`, `copy()`

**`dotPlot` class — base pairing probability ensemble:**
- `length`: RNA length
- `name`: label
- `dp`: dict with keys `'i'`, `'j'`, `'logBP'`, `'coverage'`
  - `i`, `j`: nucleotide indices of a potential base pair
  - `logBP`: `-log10(probability)` (smaller = more probable)
  - `coverage`: how many windows contributed to this entry
- Key methods: `readDP()`, `writeDP()`, `pairList()`, `requireProb()`, `trimEnds()`, `calcShannon()`, `averageSlippedBPs()`

**Module-level utilities:** `padCT()`, `readSHAPE()`, `writeSHAPE()`, `readSeq()`

---

### batchSubmit.py — Parallel Job Execution

Wraps Python `multiprocessing` to run many RNAstructure commands concurrently.

- `batchSubmit(jobCommands, nproc=4)` — main entry; submits jobs in batches of `nproc`, polls for completion
- `submitJob(command)` — executes a single shell command via `subprocess`
- `updateJobs(currentJobs, jobTrack)` — polls `AsyncResult` handles to track completion
- Can also be run standalone from CLI with a file of commands

---

### pvclient.py — PseudoViewer Web Client

Contacts the PseudoViewer SOAP web service to render RNA secondary structures as PostScript, then enhances the EPS output with SHAPE/differential reactivity coloring.

- Depends on `httplib2` for HTTP
- `python_client(CT, shape, ...)` — main entry point; takes CT object + SHAPE data, returns annotated EPS + XRNA
- `getEPS(title, startNuc, seq, struct)` — raw SOAP call to PseudoViewer
- `modifyEPS()` — post-processes EPS to inject color strings, titles
- `makeColorStrings()` — dispatches to SHAPE / differential / range color formatters

**Note:** This module has an optional import in Superfold.py (wrapped in `try/except`) so the pipeline continues if PseudoViewer is unavailable.

---

### PyCircleCompareSF.py — Circular Structure Plots

Generates circular RNA secondary structure diagrams in PostScript format. Can compare a predicted structure to a reference structure and color-annotate differences.

- `makeCircle(ct_pred, ct_correct, shannon, shape, ...)` — primary entry; writes PostScript to disk
- `compareRNA(predRNA, correctRNA)` — classifies pairs as correct, missing, extra, or slipped
- Contains a local copy of the `CT` class (independent of RNAtools.py for portability)
- Annotates nucleotides by SHAPE reactivity and Shannon entropy

---

### drawArcRibbons_simple.py — Arc Ribbon Plots

Creates publication-quality arc diagrams of RNA secondary structures using matplotlib Bezier curves.

- `writeManyPlot(outPath, title, seq, pairedNucArr, arcColors, alpha, maxDistance)` — primary interface; draws multiple structure overlays
- `writePlot()` — single or double structure comparison plot
- `plotArcRibbons(pairedNuc, ax, color, alpha)` — draws each helix as a colored arc
- `parseCTFile(lines)` — parses CT file text into seq + pairing arrays
- Output: PDF or PNG

---

### differenceByWindowSHAPEMAP.py — Differential SHAPE Analysis

Standalone utility for comparing SHAPE reactivity between two experimental conditions. Uses windowed scaling and normalization to account for global reactivity differences before computing per-nucleotide differences.

- `scaleSampleReactReport()` — main entry; scales two SHAPE datasets and reports differences
- `scaleShapeDataWindow()` — window-based linear scaling
- `normSimple()` — percentile-based normalization
- `fitLinear()` — linear interpolation for scaling factors
- Not called by the main Superfold pipeline; used independently

---

## External Dependencies

### Required Software (must be in PATH)
| Executable | Purpose |
|---|---|
| `Fold` | Minimum free energy structure prediction (RNAstructure) |
| `partition` | Partition function calculation (RNAstructure) |
| `ProbabilityPlot` | Convert partition function to dot plot (RNAstructure) |

Also requires `DATAPATH` environment variable pointing to RNAstructure thermodynamic data tables.

### Python Libraries
| Library | Use |
|---|---|
| `numpy` | Array operations, statistics, probability calculations |
| `matplotlib` | Arc plots, Shannon entropy plots, PDF/PNG output |
| `httplib2` | HTTP requests to PseudoViewer web service (pvclient only) |
| Standard library | `argparse`, `subprocess`, `multiprocessing`, `os`, `sys`, `re`, `math`, `time`, `copy`, `hashlib`, `shlex` |

---

## Internal Module Dependency Graph

```
Superfold.py
├── RNAtools.py           (CT, dotPlot, padCT, writeSHAPE)
├── batchSubmit.py        (batchSubmit)
├── drawArcRibbons_simple.py (writeManyPlot)
├── PyCircleCompareSF.py  (makeCircle)
└── pvclient.py           (python_client)  [optional — try/except import]

pvclient.py
└── RNAtools.py           (CT class)

PyCircleCompareSF.py
└── [local CT copy — no import of RNAtools]

drawArcRibbons_simple.py
└── [no internal deps]

differenceByWindowSHAPEMAP.py
└── [standalone — no internal deps]

batchSubmit.py
└── [no internal deps]
```

---

## Key Design Patterns

- **Windowed parallelism**: RNA is split into overlapping windows; each window is folded independently via `RNAstructure` executables, then results are merged. This scales to arbitrarily long RNAs.
- **External process dispatch**: All heavy computation is delegated to RNAstructure CLI executables via `subprocess`. Python manages orchestration only.
- **CT + dotPlot as universal currency**: All intermediate structure data is exchanged as `CT` (single structure) or `dotPlot` (ensemble) objects defined in RNAtools.py.
- **Optional visualization**: `pvclient.py` is imported with `try/except`; if PseudoViewer is unreachable the rest of the pipeline continues.
- **Duplicate CT class**: `PyCircleCompareSF.py` contains a local copy of `CT` rather than importing from `RNAtools.py`, making it independently portable but creating a maintenance burden.
