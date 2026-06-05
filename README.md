# SuperFold

*Gregg Rice, 2014 — gmr@unc.edu*

SHAPE-directed RNA secondary structure prediction using windowed partition function and MFE folding.

- See [ARCHITECTURE.md](./ARCHITECTURE.md) for code structure.

---

## Requirements

- Python 2.7
- RNAStructure
- matplotlib
- httplib2 (optional, for structure drawing)

### RNAStructure Setup

Follow the instructions on the [RNAStructure website](https://rna.urmc.rochester.edu/RNAstructure.html) for your platform.

### Python Modules (matplotlib, httplib2)

For each module: download the source, extract, then install:

```bash
python setup.py install --user
```

---

## Usage

```bash
python SuperFold.py RNA.map
```

All flags are optional. See `--help` for details:

```bash
python SuperFold.py --help
```

---

## Input Files

### `.map` — SHAPE Reactivity (required)

Generated automatically by the ShapeMapper pipeline. Columns: nucleotide #, SHAPE reactivity, error, sequence. `T` is converted to `U` automatically.

```
1   0.002512    0.053798    G
2   -0.034906   0.143529    T
3   -0.077852   0.257623    T
```

### `.mapd` — Differential SHAPE (optional)

Columns: nucleotide #, differential SHAPE reactivity, std error, sequence, Z-factor.

```
1   -999.0  -999.0  G   -999.0
2   -0.0124 0.2673  U   -74.244
3    0.0951 0.0833  U   -2.349
```

Create a `.mapd` file with:

```bash
python differenceByWindowSHAPEMAP.py nmia.map 1m6.map nmia-1m6.mapd 25
```

Pass to SuperFold via `--differentialFile`.

### `ssConstraints.txt` — Single-Strand Constraints (optional)

One nucleotide index per line:

```
34
35
36
```

### `ListofPKs_ds.txt` — Pseudoknot Constraints (optional)

Paired nucleotide indices, one pair per line. Used to reassemble pseudoknotted positions in the final step:

```
34 78
35 77
36 76
```

### DMS data (ShapeMapper >= 2.2.0)

If data was generated with `--dms`, run SuperFold with `--DMS` to use fold/partition commands compatible with DMS SM 2.2+ data.

---

## Output

SuperFold creates folders automatically and appends a cryptographic hash to folder/file names to avoid collisions. A log file is written to the results folder.

| Folder       | Contents |
|--------------|----------|
| `partition/` | Intermediate partition function calculations |
| `fold/`      | Intermediate MFE fold calculations |
| `results/`   | Merged partition function and MFE structures (`merged.*`) |
| `regions/`   | Per-region structure files and plots |

**Arc plot key** (base pair probabilities from partition function):

| Color  | Probability (%) |
|--------|-------------|
| Green  | > 80 |
| Blue   | > 30 |
| Yellow | > 10 |
| Gray   | > 3 |

The `ShannonSHAPE` PDF shows Shannon entropy and SHAPE analysis. Region cut sites are written to the log file.

---

## Troubleshooting

If no base pairs appear in the output, try a smaller window size:

```bash
--partitionWindowSize 1000
--foldWindowSize 1000
```

For window sizes under 1000, set `--trimInterior 200` to include interior windows. Note that smaller window sizes bias toward shorter-range interactions.
