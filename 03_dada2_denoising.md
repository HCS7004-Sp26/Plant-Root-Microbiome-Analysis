# Module 3 — DADA2 Denoising and Feature Table Construction

## From Raw Reads to Amplicon Sequence Variants (ASVs)

---

## Background

### From Reads to ASVs: What DADA2 Actually Does

Raw amplicon reads contain two sources of sequence heterogeneity: (1) true
biological variation among the microorganisms in your sample, and (2)
sequencing errors introduced by the Illumina instrument. Separating these
two signals is the central challenge of amplicon analysis.

**DADA2** (Divisive Amplicon Denoising Algorithm 2) addresses this by
building a statistical error model from your data — learning the specific
error rates of your sequencing run — and then using that model to infer
which sequences are true biological variants and which are instrument
artifacts. The result is a set of **Amplicon Sequence Variants (ASVs)**:
exact biological sequences with sequencing error removed.

The DADA2 workflow has four stages:

```
Paired-end reads
      │
      ▼
1. Quality filtering and truncation
   (remove low-quality bases; truncate at specified positions)
      │
      ▼
2. Error learning
   (estimate per-base substitution error rates from a subset of reads)
      │
      ▼
3. Dereplication + sample inference (denoising)
   (identify true ASVs using the error model;
    pool or pseudo-pool information across samples)
      │
      ▼
4. Paired-end merging and chimera removal
   (merge R1 and R2; remove chimeric sequences)
      │
      ▼
Feature table (ASV × sample count matrix)  +  Representative sequences (ASV FASTA)
```

### ASVs vs. OTUs: A Critical Methodological Note

Edwards et al. (2015) used **OTU clustering at 97% identity** — the standard
at the time. We reanalyze their data with DADA2's ASV approach because:

| Property | OTUs (97%) | ASVs (DADA2) |
| --- | --- | --- |
| Resolution | Groups sequences within 3% | Single-nucleotide resolution |
| Reproducibility | Reference-dependent | Reproducible across studies |
| Sensitivity | Loses rare diversity | Preserves rare variants |
| Chimera handling | Post-hoc filtering | Integrated into model |

This means your results will show more features (ASVs) than the original
paper reported (OTUs), and some numbers will differ. This is biologically
meaningful: ASVs are a more accurate representation of the community.

---

## Learning Objectives

* Select appropriate DADA2 truncation lengths from the demux summary visualization
* Run QIIME2's `dada2 denoise-paired` on the 16-sample subset
* Inspect and interpret the denoising statistics table
* Summarize and explore the resulting feature table
* Understand the feature table as the foundational data object for all downstream analyses

---

## Step 1: Confirm Your Truncation Lengths

Before running DADA2, confirm the truncation lengths you identified from
the `demux-summary.qzv` in Module 2. Record these in your lab notebook.

Key criteria:
1. Truncated reads must overlap by at least **20 bp** after merging
   (V4 amplicon ≈ 253 bp; minimum retained length: 253 + 20 = 273 bp)
2. Truncate where median quality drops below **Q25**

For this tutorial's MiSeq 2×251 bp data, typical values are:

```
--p-trunc-len-f 230   (R1: truncate at position 230)
--p-trunc-len-r 200   (R2: truncate at position 200)
```

> **Adjust these values based on your own quality plot.** The values above
> are a starting point for 2×251 bp data with high-quality R1 and moderate
> quality R2 tail. If your R2 drops earlier, use a shorter truncation to avoid
> retaining low-quality bases that increase chimera formation.

---

## Step 2: Run DADA2 Denoising

DADA2 is the most compute-intensive step in the pipeline. It builds
a per-base error model for each run and then corrects errors across
all samples. We use 16 threads to parallelize the sample-level denoising.

> **Expected runtime:** ~45–90 minutes for 16 samples on 16 cores.

```bash
cat > ${MICROBIOME}/scripts/03a_dada2_denoise.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=dada2_denoise
#SBATCH --output=logs/dada2_denoise_%j.out
#SBATCH --error=logs/dada2_denoise_%j.err
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
SHARED_Q2=/fs/scratch/PAS3260/Microbiome
Q2_CONTAINER=${SHARED_Q2}/Containers/qiime2.sif

# ---- Truncation lengths: adjust based on your demux-summary.qzv ----
TRUNC_F=230
TRUNC_R=200

echo "=== DADA2 Denoising ==="
echo "Started: $(date)"
echo "Truncation: R1=${TRUNC_F}  R2=${TRUNC_R}"
echo "Threads: 16"
echo ""

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs /data/02_qiime2/import/demux-paired-end.qza \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-trunc-len-f ${TRUNC_F} \
    --p-trunc-len-r ${TRUNC_R} \
    --p-n-threads 16 \
    --o-table             /data/02_qiime2/denoising/table.qza \
    --o-representative-sequences /data/02_qiime2/denoising/rep-seqs.qza \
    --o-denoising-stats   /data/02_qiime2/denoising/denoising-stats.qza \
    --verbose

echo ""
echo "=== DADA2 complete ==="
echo "Outputs:"
echo "  table.qza              → ASV × sample count matrix"
echo "  rep-seqs.qza           → representative sequences (one per ASV)"
echo "  denoising-stats.qza   → per-step read retention statistics"
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/03a_dada2_denoise.sh
squeue -u ${USER}
```

---

## Step 3: Visualize Denoising Statistics

The denoising stats artifact tells you how many reads were retained at each
step for each sample. Low retention at any step signals a problem.

```bash
cat > ${MICROBIOME}/scripts/03b_denoising_viz.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=dada2_viz
#SBATCH --output=logs/dada2_viz_%j.out
#SBATCH --error=logs/dada2_viz_%j.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
SHARED_Q2=/fs/scratch/PAS3260/Microbiome
Q2_CONTAINER=${SHARED_Q2}/Containers/qiime2.sif

echo "=== Denoising stats visualization ==="

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime metadata tabulate \
    --m-input-file /data/02_qiime2/denoising/denoising-stats.qza \
    --o-visualization /data/03_visualizations/denoising-stats.qzv

echo "=== Feature table summary ==="

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime feature-table summarize \
    --i-table /data/02_qiime2/denoising/table.qza \
    --m-sample-metadata-file /data/02_qiime2/metadata.tsv \
    --o-visualization /data/03_visualizations/table-summary.qzv

echo "=== Representative sequence summary ==="

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime feature-table tabulate-seqs \
    --i-data /data/02_qiime2/denoising/rep-seqs.qza \
    --o-visualization /data/03_visualizations/rep-seqs-summary.qzv

echo ""
echo "Visualizations written to ${MICROBIOME}/03_visualizations/"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/03b_denoising_viz.sh
squeue -u ${USER}
```

Transfer and open the visualizations:

```bash
# Run on your LOCAL machine
scp <your_username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/Microbiome/03_visualizations/denoising-stats.qzv ~/Desktop/
scp <your_username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/Microbiome/03_visualizations/table-summary.qzv ~/Desktop/
```

---

## Step 4: Interpreting Denoising Statistics

Open `denoising-stats.qzv` at [view.qiime2.org](https://view.qiime2.org).
You will see a table with the following columns:

| Column | What it means | Warning threshold |
| --- | --- | --- |
| `input` | Reads entering DADA2 | — |
| `filtered` | After quality filtering + truncation | < 70% of input → consider relaxing truncation |
| `denoised` | After error model correction | — |
| `merged` | After R1/R2 merging | < 50% of filtered → overlap may be insufficient |
| `non-chimeric` | After chimera removal | < 75% of merged → chimera rate is high |
| `passed-filter` | Final retained reads | < 1,000 per sample → may need to be excluded |

**Typical expected values for this MiSeq 2×251 bp dataset:**

- 85–95% of input reads pass quality filtering (longer reads, excellent quality)
- 90–97% of filtered reads merge successfully (generous overlap at 2×251 bp)
- 90–98% of merged reads are non-chimeric
- Final retention: ~75–85% of input reads

> **If merging efficiency is low** (< 50%): the R1 and R2 truncated reads
> do not overlap sufficiently. Reduce either truncation length to increase
> the remaining sequence length. Remember: `TRUNC_F + TRUNC_R` must exceed
> the amplicon length (253 bp) by at least 20 bp.

> **If chimera rate is high** (> 25% chimeric): this is common in highly
> diverse samples (e.g., bulk soil) and is not necessarily a problem —
> DADA2 is correctly identifying and removing chimeric sequences.

---

## Step 5: Interpreting the Feature Table Summary

Open `table-summary.qzv` at [view.qiime2.org](https://view.qiime2.org).

Key statistics to record for Module 5 (rarefaction depth selection):

1. **Total features (ASVs):** How many unique ASVs were detected across all 16 samples?
2. **Per-sample feature counts:** Minimum, median, and maximum reads per sample
3. **Sample frequency distribution:** Are any samples dramatically under-sampled?

The **interactive sample frequency detail** tab lets you see exactly how
many reads each sample retained. Sort by `compartment` in the metadata
panel to compare across compartments.

> **Expected pattern:** Endosphere samples should have fewer reads than
> bulk soil samples, because the endosphere is less densely colonized and
> yields less microbial DNA relative to plant DNA. This is a real biological
> signal, not a technical artifact — but it means you will need to choose a
> rarefaction depth that retains endosphere samples without discarding too
> many reads from the soil samples. You will make this decision in Module 5.

**Record these values for later use:**

```
Total ASVs:               ___________
Minimum reads per sample: ___________  (sample ID: _________)
Median reads per sample:  ___________
Maximum reads per sample: ___________
```

---

## Expected Files After Module 3

```
${MICROBIOME}/02_qiime2/denoising/
├── table.qza                  ← ASV × sample count matrix (feature table)
├── rep-seqs.qza               ← representative sequences (one FASTA per ASV)
└── denoising-stats.qza        ← per-sample read retention statistics

${MICROBIOME}/03_visualizations/
├── denoising-stats.qzv        ← per-step read retention table
├── table-summary.qzv          ← feature table summary with interactive plots
└── rep-seqs-summary.qzv       ← ASV sequence lengths and frequencies
```

---

## Checkpoint: Before You Proceed

* All three denoising output `.qza` files exist in `02_qiime2/denoising/`
* DADA2 log shows no fatal errors (check `logs/dada2_denoise_*.out`)
* You have opened `denoising-stats.qzv` and confirmed > 60% read retention for all samples
* No samples have fewer than 1,000 non-chimeric reads (if any do, note them for exclusion)
* You have opened `table-summary.qzv` and recorded the minimum, median, and maximum per-sample read counts
* You understand the difference between an ASV and an OTU and can explain why ASVs are methodologically preferable

---

*Previous: [Module 2 → Data and QC](02_data_and_qc.md) | Next: [Module 4 → Taxonomic Classification](04_taxonomy.md)*