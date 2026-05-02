# 🧬 NovaSNP — Novel Variant Calling Pipeline

> A lightweight, educational variant calling pipeline built from scratch in Python — no GATK, no black boxes.

![Python](https://img.shields.io/badge/Python-3.8+-blue?style=flat-square&logo=python)
![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-green?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)

---

## 🔬 What is NovaSNP?

**NovaSNP** is a fully self-contained variant calling pipeline that demonstrates the core concepts behind real-world tools like GATK HaplotypeCaller and FreeBayes — built entirely from first principles in Python.

It simulates the complete NGS data flow:

```
Reference Genome → Read Simulation → Pileup → Variant Calling → QC → VCF Output
```

---

## 🧩 Pipeline Modules

| Module | Description |
|--------|-------------|
| `ReferenceGenome` | Generates a simulated reference genome with 5 annotated gene regions (BRCA1, TP53, KRAS, EGFR, PIK3CA) |
| `NGSSimulator` | Simulates Illumina-style reads with tunable coverage, error rate, and planted true variants |
| `PileupEngine` | Builds per-position base pileups from aligned reads (analogous to `samtools mpileup`) |
| `VariantCaller` | Quality-aware allele counting with SNP detection, depth/AF/BQ filtering |
| `QCReporter` | Computes sensitivity, specificity, mean depth, allele frequency stats |
| `VCFWriter` | Standard VCF 4.2 output with INFO fields |
| `TSVWriter` | Human-readable TSV for downstream analysis |
| `JSONReporter` | JSON QC report for programmatic parsing |

---

## ⚙️ Installation

```bash
git clone https://github.com/kirthi-m-hangallmath/novasnp-pipeline.git
cd novasnp-pipeline
python variant_caller.py
```

No dependencies beyond Python standard library! _(Optional: add numpy/matplotlib for extended visualization)_

---

## 🚀 Usage

```bash
# Basic run (defaults: 30× coverage, 10kb genome)
python variant_caller.py

# Custom parameters
python variant_caller.py \
  --coverage 50 \
  --genome-length 20000 \
  --error-rate 0.005 \
  --min-depth 15 \
  --min-af 0.10 \
  --output-dir my_results
```

### Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--coverage` | 30 | Sequencing depth (×) |
| `--genome-length` | 10000 | Reference genome length (bp) |
| `--error-rate` | 0.01 | Per-base sequencing error rate |
| `--min-depth` | 10 | Minimum depth to call a variant |
| `--min-af` | 0.05 | Minimum allele frequency threshold |
| `--output-dir` | novasnp_output | Output directory |

---

## 📂 Output Files

```
novasnp_output/
├── variants.vcf            # Standard VCF 4.2 format
├── variants_summary.tsv    # Tabular PASS variant summary
└── qc_report.json          # Pipeline QC metrics
```

### VCF Example
```
##fileformat=VCFv4.2
##source=NovaSNP_Pipeline_v1.0
#CHROM  POS   ID  REF ALT QUAL  FILTER  INFO
chr1    1042  .   G   A   35.2  PASS    DP=28;AF=0.214;ANN=Missense;GENE=BRCA1;IMPACT=MODERATE
chr1    2317  .   C   T   33.8  PASS    DP=31;AF=0.387;ANN=Synonymous;GENE=TP53;IMPACT=LOW
```

---

## 🔍 Variant Filtering Logic

```
LowDepth    → total_depth < min_depth
LowAltDepth → alt_depth < 3
LowAF       → allele_freq < min_af
LowQual     → mean_base_quality < 20
PASS        → all thresholds met
```

---

## 📊 Annotation System

Variants are annotated based on gene region and base-change properties:

| Annotation | Impact | Description |
|------------|--------|-------------|
| Missense | MODERATE | Amino acid change (transition in coding region) |
| Nonsense | HIGH | Premature stop codon (transversion in coding region) |
| Synonymous | LOW | Silent mutation (third codon position) |
| Intergenic | MODIFIER | Outside annotated gene regions |

---

## 📈 Interactive Dashboard

Open `dashboard.html` in any browser to explore the pipeline results visually:

- Allele Frequency distribution histogram
- Read Depth distribution
- Filter status donut chart
- Variant impact & annotation breakdown
- Interactive PASS variant table

---

## 🧪 Concepts Demonstrated

- NGS read simulation with realistic quality score modeling
- Pileup construction from aligned reads
- Quality-aware allele frequency estimation
- VCF 4.2 format generation
- Sensitivity analysis against ground truth
- Gene-level variant annotation
- CLI tool design with argparse

---

## 👩‍💻 Author

**Kirthi M Hangallmath**  
Genome Analyst · B.E. Biotechnology (VTU, 2026)  
📧 kirthihanagalmath@gmail.com  
🔗 [LinkedIn](https://linkedin.com/in/kirthi-m-hangallmath)  
📍 Bengaluru, Karnataka

---

## 📄 License

MIT License — free to use, modify, and distribute.
