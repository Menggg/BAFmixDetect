# BAFmixDetect

A lightweight, dependency-free R implementation for B-allele frequency (BAF) regression-based DNA mixture detection from Illumina SNP array data.

This implementation is an enhanced all-in-one version of the original **bafRegress** workflow, providing native support for Illumina Final Reports and VCF files without requiring external Python preprocessing.

---

# Features

- Pure base R implementation
- No external R packages required
- Supports Illumina FinalReport
- Supports VCF input
- Automatic sample detection
- Optional BAF recomputation
    - normalized X/Y
    - raw intensity
- Multi-sample batch processing
- Conservative built-in SNP quality control
- Automatic logging
- Automatic run management
- Publication-ready summary tables

---

# Quality Control

The pipeline performs conservative SNP filtering before regression.

Default filters include

| Filter | Default |
|---------|---------|
| Missing marker | Remove |
| Missing MAF | Remove |
| Missing BAF | Remove |
| Invalid BAF (<0 or >1) | Remove |
| Duplicate markers | Remove |
| Input CNV records | Remove |
| Non-autosomal markers | Remove |
| Low GenCall Score | 0.15 |
| Low GC Score | 0.15 |
| Low intensity | X+Y ≥ 1 |
| Sample Call Rate | ≥0.90 |

All thresholds are user configurable.

---

# Installation

Clone the repository

```bash
git clone https://github.com/Menggg/BAFmixDetect.git
cd BAFmixDetect
```

No package installation is required.

Only R (≥4.1 recommended) is needed.

---

# Quick Start

```bash
Rscript run_bafRegress.R \
    --final FinalReport.txt \
    --freqfile population_frequency.txt
```

---

# Input

Supported formats

- Illumina Final Report
- VCF

Population frequency file

```
Marker    MAF
rs10001   0.32
rs10002   0.11
...
```

---

# Output

```
baf_out/

    latest/

    logs/

    plots/

    tables/

        estimate_summary.tsv

        qc_summary.tsv

        samples.txt

    run_info.tsv

    cmd.sh
```

---

# QC Summary

Each sample receives a QC report including

- Total SNPs
- Duplicate markers removed
- Invalid BAF removed
- Missing MAF removed
- CNV markers removed
- Low GenCall removed
- Low GC Score removed
- Low intensity removed
- Non-autosomal markers removed
- Final usable SNPs
- Homozygous SNPs
- Call rate

---

# Optional Parameters

## MAF

```bash
--min_maf 0
```

Official bafRegress-compatible default.

---

## GenCall Score

```bash
--min_gencall 0.15
```

---

## GC Score

```bash
--min_gcscore 0.15
```

---

## Intensity

```bash
--min_intensity 1
```

Uses

```
X + Y
```

---

## Sample Call Rate

```bash
--min_callrate 0.90
```

---

## Autosomes only

```bash
--autosome_only TRUE
```

---

## Disable QC

Example

```bash
--autosome_only FALSE

--min_gencall 0

--min_gcscore 0

--min_intensity 0
```

---

# BAF Recalculation

The pipeline can optionally recalculate BAF directly from probe intensities.

Supported modes

```
original

xy_y

xy_x

raw_y

raw_x
```

Example

```bash
--rebaf_mode raw_y
```

---

# Compatibility

The regression model is fully compatible with the original **bafRegress** implementation.

Default MAF behavior follows the original software

```
MAF > 0
```

rather than imposing additional allele-frequency thresholds.

---

# Citation

If you use this software, please cite

> (Manuscript in preparation)

---

---

# License

Copyright (c) 2026 Meng Wang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

MIT License
