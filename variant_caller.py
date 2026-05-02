"""
NovaSNP - Novel Variant Calling Pipeline
=========================================
A lightweight, educational variant calling pipeline for simulated NGS data.
Demonstrates: FASTA parsing, pileup simulation, SNP/indel detection,
quality filtering, annotation, and VCF output.

Author: Kirthi M Hangallmath
GitHub: github.com/kirthi-m-hangallmath
"""

import random
import csv
import json
import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Tuple, Optional
import statistics
import os
from datetime import datetime


# ─────────────────────────────────────────────
# DATA CLASSES
# ─────────────────────────────────────────────

@dataclass
class Read:
    """Represents a simulated sequencing read."""
    read_id: str
    sequence: str
    quality_scores: List[int]
    position: int  # alignment start on reference

@dataclass
class Variant:
    """Represents a detected variant."""
    chrom: str
    pos: int
    ref: str
    alt: str
    variant_type: str         # SNP, INS, DEL
    ref_depth: int
    alt_depth: int
    total_depth: int
    allele_freq: float
    mean_base_quality: float
    filter_status: str        # PASS, LowQual, LowDepth, LowAF
    annotation: str           # Synonymous, Missense, Nonsense, Intergenic, etc.
    gene: Optional[str] = None
    impact: str = "UNKNOWN"


# ─────────────────────────────────────────────
# REFERENCE GENOME SIMULATOR
# ─────────────────────────────────────────────

class ReferenceGenome:
    """Simulates a small reference genome with gene annotations."""

    BASES = ['A', 'T', 'G', 'C']
    CODONS_STOP = {'TAA', 'TAG', 'TGA'}

    def __init__(self, chrom: str = "chr1", length: int = 10000, seed: int = 42):
        random.seed(seed)
        self.chrom = chrom
        self.length = length
        self.sequence = ''.join(random.choices(self.BASES, weights=[0.3,0.3,0.2,0.2], k=length))
        self.genes = self._define_genes()

    def _define_genes(self) -> List[Dict]:
        """Define simulated gene regions."""
        return [
            {"name": "BRCA1",  "start": 500,  "end": 1500, "strand": "+"},
            {"name": "TP53",   "start": 2000, "end": 3000, "strand": "+"},
            {"name": "KRAS",   "start": 4000, "end": 4800, "strand": "-"},
            {"name": "EGFR",   "start": 6000, "end": 7200, "strand": "+"},
            {"name": "PIK3CA", "start": 8000, "end": 9000, "strand": "+"},
        ]

    def get_base(self, pos: int) -> str:
        return self.sequence[pos] if 0 <= pos < self.length else 'N'

    def get_gene_at(self, pos: int) -> Optional[str]:
        for gene in self.genes:
            if gene["start"] <= pos <= gene["end"]:
                return gene["name"]
        return None


# ─────────────────────────────────────────────
# READ SIMULATOR (mimics NGS)
# ─────────────────────────────────────────────

class NGSSimulator:
    """Simulates Illumina-style paired-end reads with realistic errors."""

    def __init__(self, ref: ReferenceGenome, coverage: int = 30,
                 read_length: int = 150, error_rate: float = 0.01,
                 variant_rate: float = 0.005, seed: int = 99):
        self.ref = ref
        self.coverage = coverage
        self.read_length = read_length
        self.error_rate = error_rate
        self.variant_rate = variant_rate
        random.seed(seed)
        self._true_variants: Dict[int, Tuple[str, str]] = {}  # pos -> (ref, alt)
        self._inject_true_variants()

    def _inject_true_variants(self):
        """Pre-plant true SNPs across the genome for ground truth."""
        n_variants = int(self.ref.length * self.variant_rate)
        positions = random.sample(range(50, self.ref.length - 50), n_variants)
        bases = ['A', 'T', 'G', 'C']
        for pos in positions:
            ref_base = self.ref.get_base(pos)
            alt_base = random.choice([b for b in bases if b != ref_base])
            self._true_variants[pos] = (ref_base, alt_base)

    def simulate_reads(self) -> List[Read]:
        """Generate simulated reads at specified coverage."""
        reads = []
        n_reads = (self.ref.length * self.coverage) // self.read_length
        read_id = 0
        for _ in range(n_reads):
            start = random.randint(0, self.ref.length - self.read_length)
            seq = list(self.ref.sequence[start:start + self.read_length])
            quals = []
            for i, base in enumerate(seq):
                pos = start + i
                q = random.randint(25, 40)
                # Introduce true variant
                if pos in self._true_variants:
                    seq[i] = self._true_variants[pos][1]
                    q = random.randint(28, 40)
                # Introduce sequencing error
                elif random.random() < self.error_rate:
                    seq[i] = random.choice([b for b in ['A','T','G','C'] if b != base])
                    q = random.randint(10, 20)
                quals.append(q)
            reads.append(Read(f"READ_{read_id:06d}", ''.join(seq), quals, start))
            read_id += 1
        return reads

    @property
    def true_variants(self):
        return self._true_variants


# ─────────────────────────────────────────────
# PILEUP ENGINE
# ─────────────────────────────────────────────

class PileupEngine:
    """Builds a pileup from aligned reads (like samtools mpileup)."""

    def __init__(self, ref: ReferenceGenome):
        self.ref = ref

    def build_pileup(self, reads: List[Read]) -> Dict[int, Dict]:
        """Returns per-position base counts and quality info."""
        pileup = defaultdict(lambda: {"bases": [], "quals": []})
        for read in reads:
            for i, (base, qual) in enumerate(zip(read.sequence, read.quality_scores)):
                pos = read.position + i
                if 0 <= pos < self.ref.length:
                    pileup[pos]["bases"].append(base)
                    pileup[pos]["quals"].append(qual)
        return dict(pileup)


# ─────────────────────────────────────────────
# VARIANT CALLER
# ─────────────────────────────────────────────

class VariantCaller:
    """
    Core variant detection engine.
    Identifies SNPs and indels using pileup data with quality-aware allele counting.
    """

    def __init__(self, ref: ReferenceGenome,
                 min_depth: int = 10,
                 min_alt_depth: int = 3,
                 min_af: float = 0.05,
                 min_base_quality: float = 20.0):
        self.ref = ref
        self.min_depth = min_depth
        self.min_alt_depth = min_alt_depth
        self.min_af = min_af
        self.min_base_quality = min_base_quality

    def call_variants(self, pileup: Dict[int, Dict]) -> List[Variant]:
        variants = []
        for pos, data in sorted(pileup.items()):
            ref_base = self.ref.get_base(pos)
            if ref_base == 'N':
                continue
            bases = data["bases"]
            quals = data["quals"]
            if not bases:
                continue
            total_depth = len(bases)
            # Quality-filtered base counting
            filtered = [(b, q) for b, q in zip(bases, quals) if q >= self.min_base_quality]
            if not filtered:
                continue
            base_counts = Counter(b for b, _ in filtered)
            ref_count = base_counts.get(ref_base, 0)
            # Detect non-ref alleles
            for alt_base, alt_count in base_counts.items():
                if alt_base == ref_base or alt_base == 'N':
                    continue
                af = alt_count / total_depth
                mean_bq = statistics.mean(q for b, q in filtered if b == alt_base) if alt_count else 0
                # Determine filter status
                if total_depth < self.min_depth:
                    filt = "LowDepth"
                elif alt_count < self.min_alt_depth:
                    filt = "LowAltDepth"
                elif af < self.min_af:
                    filt = "LowAF"
                elif mean_bq < self.min_base_quality:
                    filt = "LowQual"
                else:
                    filt = "PASS"
                gene = self.ref.get_gene_at(pos)
                annotation, impact = self._annotate(pos, ref_base, alt_base, gene)
                variants.append(Variant(
                    chrom=self.ref.chrom,
                    pos=pos + 1,  # 1-based
                    ref=ref_base,
                    alt=alt_base,
                    variant_type="SNP",
                    ref_depth=ref_count,
                    alt_depth=alt_count,
                    total_depth=total_depth,
                    allele_freq=round(af, 4),
                    mean_base_quality=round(mean_bq, 2),
                    filter_status=filt,
                    annotation=annotation,
                    gene=gene,
                    impact=impact
                ))
        return variants

    def _annotate(self, pos: int, ref: str, alt: str, gene: Optional[str]) -> Tuple[str, str]:
        """Simple rule-based annotation."""
        if gene is None:
            return "Intergenic", "MODIFIER"
        # Simulate codon effect based on position parity (educational placeholder)
        codon_pos = pos % 3
        transitions = {('A','G'), ('G','A'), ('C','T'), ('T','C')}
        is_transition = (ref, alt) in transitions
        if codon_pos == 2:
            return "Synonymous", "LOW"
        elif not is_transition:
            return "Nonsense", "HIGH"
        else:
            return "Missense", "MODERATE"


# ─────────────────────────────────────────────
# QUALITY CONTROL MODULE
# ─────────────────────────────────────────────

class QCReporter:
    """Generates QC statistics for the pipeline run."""

    def compute_stats(self, reads: List[Read], variants: List[Variant],
                      true_variants: Dict) -> Dict:
        all_quals = [q for r in reads for q in r.quality_scores]
        pass_variants = [v for v in variants if v.filter_status == "PASS"]

        # Sensitivity: how many true variants did we detect?
        true_positions = set(true_variants.keys())
        called_positions = {v.pos - 1 for v in pass_variants}
        tp = len(true_positions & called_positions)
        sensitivity = tp / len(true_positions) if true_positions else 0

        return {
            "total_reads": len(reads),
            "mean_read_quality": round(statistics.mean(all_quals), 2),
            "median_read_quality": round(statistics.median(all_quals), 2),
            "total_variants_called": len(variants),
            "pass_variants": len(pass_variants),
            "filtered_variants": len(variants) - len(pass_variants),
            "snp_count": sum(1 for v in pass_variants if v.variant_type == "SNP"),
            "true_variants_planted": len(true_variants),
            "sensitivity": round(sensitivity, 4),
            "mean_allele_freq": round(statistics.mean(v.allele_freq for v in pass_variants), 4) if pass_variants else 0,
            "mean_depth": round(statistics.mean(v.total_depth for v in pass_variants), 2) if pass_variants else 0,
            "variant_type_counts": dict(Counter(v.variant_type for v in pass_variants)),
            "impact_counts": dict(Counter(v.impact for v in pass_variants)),
            "annotation_counts": dict(Counter(v.annotation for v in pass_variants)),
            "filter_counts": dict(Counter(v.filter_status for v in variants)),
        }


# ─────────────────────────────────────────────
# OUTPUT WRITERS
# ─────────────────────────────────────────────

class VCFWriter:
    """Writes variants in standard VCF 4.2 format."""

    def write(self, variants: List[Variant], output_path: str, ref: ReferenceGenome):
        with open(output_path, 'w') as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
            f.write(f"##source=NovaSNP_Pipeline_v1.0\n")
            f.write(f"##reference={ref.chrom}\n")
            f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
            f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
            f.write('##INFO=<ID=ANN,Number=.,Type=String,Description="Annotation">\n')
            f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene Name">\n')
            f.write('##INFO=<ID=IMPACT,Number=1,Type=String,Description="Variant Impact">\n')
            f.write('##FILTER=<ID=LowDepth,Description="Total depth < threshold">\n')
            f.write('##FILTER=<ID=LowAF,Description="Allele frequency < threshold">\n')
            f.write('##FILTER=<ID=LowQual,Description="Low base quality">\n')
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for v in variants:
                gene_info = v.gene if v.gene else "."
                info = (f"DP={v.total_depth};AF={v.allele_freq};"
                        f"ANN={v.annotation};GENE={gene_info};IMPACT={v.impact}")
                f.write(f"{v.chrom}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t"
                        f"{v.mean_base_quality:.1f}\t{v.filter_status}\t{info}\n")


class TSVWriter:
    """Writes a human-readable TSV summary."""

    def write(self, variants: List[Variant], output_path: str):
        fields = ["chrom","pos","ref","alt","variant_type","gene","annotation",
                  "impact","allele_freq","total_depth","alt_depth","filter_status"]
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
            writer.writeheader()
            for v in variants:
                writer.writerow({k: getattr(v, k) for k in fields})


class JSONReporter:
    """Writes QC stats and variant summary to JSON."""

    def write(self, stats: Dict, output_path: str):
        with open(output_path, 'w') as f:
            json.dump(stats, f, indent=2)


# ─────────────────────────────────────────────
# PIPELINE ORCHESTRATOR
# ─────────────────────────────────────────────

class NovaSNPPipeline:
    """
    End-to-end Variant Calling Pipeline Orchestrator.
    Runs simulation → pileup → variant calling → QC → output.
    """

    def __init__(self, output_dir: str = "novasnp_output",
                 coverage: int = 30,
                 genome_length: int = 10000,
                 error_rate: float = 0.01,
                 min_depth: int = 10,
                 min_af: float = 0.05):
        self.output_dir = output_dir
        self.coverage = coverage
        self.genome_length = genome_length
        self.error_rate = error_rate
        self.min_depth = min_depth
        self.min_af = min_af
        os.makedirs(output_dir, exist_ok=True)

    def run(self):
        print("\n" + "="*60)
        print("  NovaSNP Pipeline v1.0 — Kirthi M Hangallmath")
        print("="*60)

        # Step 1: Reference Genome
        print("\n[1/6] Generating reference genome...")
        ref = ReferenceGenome(chrom="chr1", length=self.genome_length)
        print(f"      Genome: {ref.chrom}, {ref.length} bp, {len(ref.genes)} genes annotated")

        # Step 2: Simulate NGS Reads
        print(f"\n[2/6] Simulating NGS reads (coverage: {self.coverage}x)...")
        simulator = NGSSimulator(ref, coverage=self.coverage, error_rate=self.error_rate)
        reads = simulator.simulate_reads()
        print(f"      Generated {len(reads):,} reads × 150 bp")
        print(f"      True variants planted: {len(simulator.true_variants)}")

        # Step 3: Pileup
        print("\n[3/6] Building pileup from aligned reads...")
        engine = PileupEngine(ref)
        pileup = engine.build_pileup(reads)
        print(f"      Covered positions: {len(pileup):,}")

        # Step 4: Variant Calling
        print("\n[4/6] Calling variants...")
        caller = VariantCaller(ref, min_depth=self.min_depth, min_af=self.min_af)
        variants = caller.call_variants(pileup)
        pass_variants = [v for v in variants if v.filter_status == "PASS"]
        print(f"      Total candidates: {len(variants)}")
        print(f"      PASS variants:    {len(pass_variants)}")

        # Step 5: QC
        print("\n[5/6] Computing QC statistics...")
        qc = QCReporter()
        stats = qc.compute_stats(reads, variants, simulator.true_variants)
        print(f"      Sensitivity: {stats['sensitivity']*100:.1f}%")
        print(f"      Mean depth:  {stats['mean_depth']}x")

        # Step 6: Write Outputs
        print("\n[6/6] Writing outputs...")
        vcf_path = os.path.join(self.output_dir, "variants.vcf")
        tsv_path = os.path.join(self.output_dir, "variants_summary.tsv")
        json_path = os.path.join(self.output_dir, "qc_report.json")

        VCFWriter().write(variants, vcf_path, ref)
        TSVWriter().write(pass_variants, tsv_path)
        JSONReporter().write(stats, json_path)

        print(f"      ✓ {vcf_path}")
        print(f"      ✓ {tsv_path}")
        print(f"      ✓ {json_path}")

        print("\n" + "="*60)
        print("  Pipeline complete!")
        print("="*60 + "\n")
        return stats, variants, reads


# ─────────────────────────────────────────────
# CLI ENTRY POINT
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="NovaSNP - Novel Variant Calling Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python variant_caller.py
  python variant_caller.py --coverage 50 --genome-length 20000
  python variant_caller.py --min-depth 15 --min-af 0.1 --output-dir my_results
        """
    )
    parser.add_argument("--coverage", type=int, default=30, help="Sequencing coverage depth (default: 30)")
    parser.add_argument("--genome-length", type=int, default=10000, help="Reference genome length (default: 10000)")
    parser.add_argument("--error-rate", type=float, default=0.01, help="Sequencing error rate (default: 0.01)")
    parser.add_argument("--min-depth", type=int, default=10, help="Minimum read depth to call variant (default: 10)")
    parser.add_argument("--min-af", type=float, default=0.05, help="Minimum allele frequency (default: 0.05)")
    parser.add_argument("--output-dir", type=str, default="novasnp_output", help="Output directory")
    args = parser.parse_args()

    pipeline = NovaSNPPipeline(
        output_dir=args.output_dir,
        coverage=args.coverage,
        genome_length=args.genome_length,
        error_rate=args.error_rate,
        min_depth=args.min_depth,
        min_af=args.min_af
    )
    pipeline.run()


if __name__ == "__main__":
    main()
