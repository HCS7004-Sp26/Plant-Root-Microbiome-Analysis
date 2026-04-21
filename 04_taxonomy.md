# Module 4 — Taxonomic Classification and Community Visualization

## Assigning Identity to ASVs Using the SILVA Reference Database

---

## Background

### What Is Taxonomic Classification?

After DADA2 produces a set of exact ASV sequences, each one is simply a
DNA string — we know it exists in the community, but not which organism
it represents. **Taxonomic classification** maps each ASV sequence to
a position in the tree of life by comparing it against a reference
database of known 16S rRNA sequences.

QIIME2 uses a **naïve Bayes classifier** trained on the SILVA 138 database
(specifically, the V4 region sequences bounded by primers 515F/806R). The
classifier assigns each ASV a taxonomy from domain down to genus (or species
when confidence is sufficient), along with a confidence score.

### The SILVA Database

**SILVA** is the most widely used 16S rRNA reference database in microbiome
research, containing several million high-quality rRNA sequences curated
from public repositories. We use SILVA release 138 because:

- It is the last release with a fully curated taxonomy aligned with NCBI taxonomy
- The QIIME2 pre-trained classifier is available for exactly the 515F/806R V4 amplicon
- It is the reference database most compatible with the original Edwards et al. analysis

The pre-trained classifier (`silva-138-99-515-806-nb-classifier.qza`) is
pre-installed in the shared class directory. **Do not attempt to re-train it.**
Training a SILVA classifier requires ~32 GB of RAM and several hours of compute.

### Why Remove Chloroplasts and Mitochondria?

Rice roots contain plant cells, and those cells contain chloroplasts and
mitochondria. Both organelles have their own 16S-like rRNA genes (because
they descend from endosymbiotic bacteria). The 515F/806R primers amplify
these sequences just as efficiently as bacterial sequences. If not removed,
chloroplast and mitochondrial ASVs will inflate your feature count and
distort diversity metrics. This is especially important for endosphere
samples, where plant-derived sequences can dominate the amplicon pool.

---

## Learning Objectives

* Classify ASVs against the SILVA 138 V4 naïve Bayes classifier
* Visualize taxonomy assignments across samples (taxa bar plot)
* Filter the feature table to remove non-bacterial contaminants (chloroplasts, mitochondria, Eukaryota)
* Build a de novo phylogenetic tree from representative sequences
* Understand how contamination filtering changes the feature table

---

## Step 1: Classify ASVs Against the SILVA Classifier

```bash
cat > ${MICROBIOME}/scripts/04a_classify_taxonomy.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=classify_taxonomy
#SBATCH --output=logs/classify_taxonomy_%j.out
#SBATCH --error=logs/classify_taxonomy_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
SHARED_Q2=/fs/scratch/PAS3260/Team_Project/Containers/QIIME2
Q2_CONTAINER=${SHARED_Q2}/qiime2_amplicon_2024.10.sif
SILVA_CLASSIFIER=${SHARED_Q2}/silva-138-99-515-806-nb-classifier.qza

echo "=== Taxonomic Classification: SILVA 138 naïve Bayes ==="
echo "Started: $(date)"
echo ""

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --bind ${SHARED_Q2}:/shared \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime feature-classifier classify-sklearn \
    --i-classifier /shared/silva-138-99-515-806-nb-classifier.qza \
    --i-reads      /data/02_qiime2/denoising/rep-seqs.qza \
    --p-n-jobs 16 \
    --o-classification /data/02_qiime2/taxonomy/taxonomy.qza \
    --verbose

echo ""
echo "=== Taxonomy visualization ==="

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime metadata tabulate \
    --m-input-file /data/02_qiime2/taxonomy/taxonomy.qza \
    --o-visualization /data/03_visualizations/taxonomy.qzv

echo ""
echo "=== Taxa bar plot (pre-filtering) ==="

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime taxa barplot \
    --i-table    /data/02_qiime2/denoising/table.qza \
    --i-taxonomy /data/02_qiime2/taxonomy/taxonomy.qza \
    --m-metadata-file /data/02_qiime2/metadata.tsv \
    --o-visualization /data/03_visualizations/taxa-barplot-unfiltered.qzv

echo ""
echo "Classification complete."
echo "Outputs: taxonomy.qza, taxonomy.qzv, taxa-barplot-unfiltered.qzv"
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/04a_classify_taxonomy.sh
squeue -u ${USER}
```

> **Expected runtime:** ~60–90 minutes for 16 samples with 16 threads.
> The naïve Bayes classifier is memory-intensive. Do not reduce `--mem` below 48G.

---

## Step 2: Inspect the Unfiltered Taxa Bar Plot

Transfer `taxa-barplot-unfiltered.qzv` to your local computer and open
it in [view.qiime2.org](https://view.qiime2.org).

```bash
# Run on your LOCAL machine
scp <your_username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/Microbiome/03_visualizations/taxa-barplot-unfiltered.qzv ~/Desktop/
```

In the interactive bar plot:
1. Change the **taxonomic level** to Level 2 (phylum) and Level 3 (class)
2. Group samples by `compartment` (use the metadata column selector)
3. Identify any bars labeled `Chloroplast` or `Mitochondria`

**Guiding questions:**

1. In which compartment(s) do you see the highest proportion of
   chloroplast/mitochondrial sequences? Why does this make biological sense?
2. At phylum level, what are the dominant bacterial phyla across all compartments?
3. Can you already see a qualitative difference in community composition
   between bulk soil and endosphere at this level of resolution?

---

## Step 3: Filter Contaminants from the Feature Table

We remove three categories of non-bacterial/non-archaeal sequences
that would otherwise confound diversity calculations:

- `Chloroplast` (plant organelles; classified under Cyanobacteria)
- `Mitochondria` (plant organelles; classified under Proteobacteria)
- `Eukaryota` (any eukaryotic sequences that passed amplification)
- Features with no taxonomic assignment at the phylum level (`p__`)

```bash
cat > ${MICROBIOME}/scripts/04b_filter_contaminants.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=filter_contaminants
#SBATCH --output=logs/filter_contaminants_%j.out
#SBATCH --error=logs/filter_contaminants_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
SHARED_Q2=/fs/scratch/PAS3260/Team_Project/Containers/QIIME2
Q2_CONTAINER=${SHARED_Q2}/qiime2_amplicon_2024.10.sif

echo "=== Filtering chloroplasts, mitochondria, and unassigned features ==="
echo "Started: $(date)"

# ---- Step 1: Remove Chloroplast and Mitochondria ----
apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime taxa filter-table \
    --i-table    /data/02_qiime2/denoising/table.qza \
    --i-taxonomy /data/02_qiime2/taxonomy/taxonomy.qza \
    --p-exclude Mitochondria,Chloroplast,Eukaryota \
    --o-filtered-table /data/02_qiime2/taxonomy/table-no-contam.qza

# ---- Step 2: Also filter the representative sequences to match ----
apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime taxa filter-seqs \
    --i-sequences /data/02_qiime2/denoising/rep-seqs.qza \
    --i-taxonomy  /data/02_qiime2/taxonomy/taxonomy.qza \
    --p-exclude Mitochondria,Chloroplast,Eukaryota \
    --o-filtered-sequences /data/02_qiime2/taxonomy/rep-seqs-no-contam.qza

# ---- Step 3: Summarize the filtered table ----
apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime feature-table summarize \
    --i-table /data/02_qiime2/taxonomy/table-no-contam.qza \
    --m-sample-metadata-file /data/02_qiime2/metadata.tsv \
    --o-visualization /data/03_visualizations/table-no-contam-summary.qzv

# ---- Step 4: Updated taxa barplot (post-filtering) ----
apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime taxa barplot \
    --i-table    /data/02_qiime2/taxonomy/table-no-contam.qza \
    --i-taxonomy /data/02_qiime2/taxonomy/taxonomy.qza \
    --m-metadata-file /data/02_qiime2/metadata.tsv \
    --o-visualization /data/03_visualizations/taxa-barplot-filtered.qzv

echo ""
echo "Filtering complete."
echo "Clean table: table-no-contam.qza"
echo "Clean rep-seqs: rep-seqs-no-contam.qza"
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/04b_filter_contaminants.sh
squeue -u ${USER}
```

---

## Step 4: Build a Phylogenetic Tree

Beta diversity metrics that incorporate evolutionary relatedness
(weighted and unweighted UniFrac) require a phylogenetic tree.
We build one de novo from the filtered representative sequences
using MAFFT alignment followed by FastTree.

```bash
cat > ${MICROBIOME}/scripts/04c_phylogeny.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=phylogeny
#SBATCH --output=logs/phylogeny_%j.out
#SBATCH --error=logs/phylogeny_%j.err
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
SHARED_Q2=/fs/scratch/PAS3260/Team_Project/Containers/QIIME2
Q2_CONTAINER=${SHARED_Q2}/qiime2_amplicon_2024.10.sif

echo "=== Phylogenetic tree construction ==="
echo "Pipeline: MAFFT alignment → masking → FastTree → midpoint rooting"
echo "Started: $(date)"

# QIIME2 convenience wrapper: align, mask, build tree, root — all in one
apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences /data/02_qiime2/taxonomy/rep-seqs-no-contam.qza \
    --p-n-threads 8 \
    --o-alignment         /data/02_qiime2/phylogeny/aligned-rep-seqs.qza \
    --o-masked-alignment  /data/02_qiime2/phylogeny/masked-aligned-rep-seqs.qza \
    --o-tree              /data/02_qiime2/phylogeny/unrooted-tree.qza \
    --o-rooted-tree       /data/02_qiime2/phylogeny/rooted-tree.qza \
    --verbose

echo ""
echo "Phylogenetic tree complete."
echo "Rooted tree: ${MICROBIOME}/02_qiime2/phylogeny/rooted-tree.qza"
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/04c_phylogeny.sh
squeue -u ${USER}
```

> **Why midpoint rooting?** A rooted tree is required for UniFrac distance
> calculations. We use midpoint rooting (placing the root at the midpoint
> of the longest branch) because we do not have an outgroup sequence in
> this dataset. This is the standard approach for amplicon data.

---

## Step 5: Explore the Filtered Taxa Bar Plot

After the filtering job completes, download and open `taxa-barplot-filtered.qzv`:

```bash
# Run on your LOCAL machine
scp <your_username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/Microbiome/03_visualizations/taxa-barplot-filtered.qzv ~/Desktop/
```

**Guided exploration task:**

At taxonomic Level 2 (phylum), group samples by `compartment`.
Fill in the table below with the visually dominant phyla in each compartment:

| Compartment | Dominant phyla (estimated %) |
| --- | --- |
| Bulk Soil | |
| Rhizosphere | |
| Rhizoplane | |
| Endosphere | |

> **Expected pattern from Edwards et al.:** Bulk soil and rhizosphere are
> dominated by Acidobacteria, Proteobacteria, Actinobacteria, and Verrucomicrobia.
> The endosphere shows a marked shift toward Proteobacteria (especially
> Gammaproteobacteria) and a dramatic reduction in Acidobacteria —
> reflecting the selective filter imposed by the plant.

**Does your analysis recapitulate this pattern?** Record your observations
before proceeding to quantitative diversity analysis in Module 5.

---

## Expected Files After Module 4

```
${MICROBIOME}/02_qiime2/taxonomy/
├── taxonomy.qza                  ← ASV taxonomy assignments
├── table-no-contam.qza           ← filtered feature table (no chloroplast/mito)
└── rep-seqs-no-contam.qza        ← filtered representative sequences

${MICROBIOME}/02_qiime2/phylogeny/
├── aligned-rep-seqs.qza          ← MAFFT multiple sequence alignment
├── masked-aligned-rep-seqs.qza   ← masked alignment (removes hypervariable sites)
├── unrooted-tree.qza             ← FastTree unrooted phylogeny
└── rooted-tree.qza               ← midpoint-rooted tree (used for UniFrac)

${MICROBIOME}/03_visualizations/
├── taxonomy.qzv                        ← taxonomy classification table
├── taxa-barplot-unfiltered.qzv         ← community composition (all ASVs)
├── taxa-barplot-filtered.qzv           ← community composition (bacteria/archaea only)
└── table-no-contam-summary.qzv         ← filtered table summary
```

---

## Checkpoint: Before You Proceed

* `taxonomy.qza`, `table-no-contam.qza`, and `rep-seqs-no-contam.qza` all exist
* `rooted-tree.qza` exists in `02_qiime2/phylogeny/`
* You have viewed `taxa-barplot-filtered.qzv` and recorded dominant phyla per compartment
* You can describe the biological reason for removing chloroplast and mitochondrial sequences
* You have compared the unfiltered and filtered bar plots and noted the change in the endosphere samples
* You have updated your notes on the **minimum per-sample feature count** from `table-no-contam-summary.qzv` (needed for rarefaction depth in Module 5)

---

*Previous: [Module 3 → DADA2 Denoising](03_dada2_denoising.md) | Next: [Module 5 → Diversity Analysis](05_diversity.md)*
