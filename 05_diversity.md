# Module 5 — Diversity Analysis

## Alpha Diversity, Beta Diversity, UniFrac, and Statistical Testing

---

## Background

This module is the quantitative heart of the tutorial. We move from
visualizing *who is present* (Module 4) to measuring *how diverse* each
community is (alpha diversity) and *how different* communities are from
each other (beta diversity). Both frameworks are needed to formally test
the three hypotheses introduced in the README.

### Alpha Diversity: How Rich Is a Single Community?

Alpha diversity describes the diversity *within* a single sample.
Several metrics capture different aspects of this:

| Metric | What it measures | Sensitive to |
| --- | --- | --- |
| **Observed Features** | Raw ASV count | Rare species; sequencing depth |
| **Shannon Entropy** | Richness + evenness; penalizes dominance | Abundant species more than rare |
| **Faith's Phylogenetic Diversity (PD)** | Total branch length of the community's phylogenetic tree | Both richness and evolutionary breadth |
| **Pielou's Evenness** | How equitably abundance is distributed | Dominance by a few taxa |

For the compartment-gradient hypothesis, we predict:
> **H1:** Faith's PD and Shannon entropy decrease monotonically from
> bulk soil → rhizosphere → rhizoplane → endosphere.

This prediction follows directly from the filtering model: each successive
compartment imposes a stricter biochemical and physical filter on the
community, reducing both the number of taxa and (by removing deep-branching
phyla like Acidobacteria) the phylogenetic breadth of the survivors.

### Beta Diversity: How Different Are Two Communities?

Beta diversity measures compositional differences *between* samples.
Different metrics capture different aspects of dissimilarity:

| Metric | What it measures | Requires tree? |
| --- | --- | --- |
| **Bray-Curtis dissimilarity** | Abundance-weighted compositional difference | No |
| **Jaccard distance** | Presence/absence compositional difference | No |
| **Unweighted UniFrac** | Phylogenetic presence/absence difference | Yes |
| **Weighted UniFrac** | Phylogenetic abundance-weighted difference | Yes |

For the compartment-gradient hypothesis, we predict:
> **H2:** Samples cluster by compartment in PCoA ordination, and
> compartment explains a significant portion of total community
> variance (PERMANOVA p < 0.05).

The UniFrac metrics are particularly important here because they
capture evolutionary divergence, not just taxonomic overlap. If host
selection enriches for specific lineages, communities should cluster
in phylogenetic space even more strongly than in taxonomy space.

### Rarefaction: Equalizing Sequencing Depth

Different samples contain different numbers of reads — some through
biology (endosphere has less microbial biomass) and some through
sequencing variation. Alpha and beta diversity metrics are sensitive
to sequencing depth: a deeply sequenced sample will appear more diverse
simply because you sampled more of it.

**Rarefaction** addresses this by randomly subsampling each sample to
a common depth (`--p-sampling-depth`). Samples with fewer reads than
the chosen depth are discarded. Choosing the depth requires a tradeoff:

- **Deeper depth** → more diversity captured per sample, but more samples lost
- **Shallower depth** → all samples retained, but rare taxa undersampled

Use the feature table summary from Module 3 and 4 to choose a depth
that retains all 16 samples while being as deep as possible. A good
rule of thumb is to rarefy to ≥80% of your minimum-depth sample's count.

---

## Learning Objectives

* Choose an appropriate rarefaction depth and justify it
* Compute core diversity metrics (alpha + beta) in a single QIIME2 command
* Produce and interpret alpha diversity box plots with statistical tests
* Produce and interpret PCoA ordination plots (Emperor)
* Run PERMANOVA to statistically test whether compartment explains community composition
* Interpret results in the context of the Edwards et al. compartment-gradient model

---

## Step 1: Choose Your Rarefaction Depth

From your `table-no-contam-summary.qzv` (Module 4), record:

```
Minimum reads per sample (after contamination filtering): ___________
```

Choose a rarefaction depth that:
1. Retains all 16 samples (depth ≤ minimum per-sample read count)
2. Is as deep as possible given constraint 1
3. Is a round number for readability (e.g., 1000, 2000, 5000)

> **For this tutorial's greenhouse M104 samples:** Endosphere samples
> typically have lower read counts after contamination filtering, because
> plant organellar sequences dominate the amplicon pool from root interior
> tissue. Check your `table-no-contam-summary.qzv` from Module 4 to find
> the actual minimum. For the MiSeq 2×251 bp runs in this tutorial, a depth
> of **3,000–5,000 reads** often retains all samples. You may need to
> exclude individual samples with very low counts. Always document your
> choice and reasoning.

Run the rarefaction curve visualization to help choose the depth:

```bash
cat > ${MICROBIOME}/scripts/05a_rarefaction.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=rarefaction
#SBATCH --output=logs/rarefaction_%j.out
#SBATCH --error=logs/rarefaction_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
SHARED_Q2=/fs/scratch/PAS3260/Microbiome
Q2_CONTAINER=${SHARED_Q2}/Containers/qiime2.sif

echo "=== Alpha rarefaction curves ==="
echo "Started: $(date)"

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime diversity alpha-rarefaction \
    --i-table      /data/02_qiime2/taxonomy/table-no-contam.qza \
    --i-phylogeny  /data/02_qiime2/phylogeny/rooted-tree.qza \
    --p-max-depth 10000 \
    --p-steps 20 \
    --m-metadata-file /data/02_qiime2/metadata.tsv \
    --o-visualization /data/03_visualizations/alpha-rarefaction.qzv

echo "Rarefaction curve visualization: alpha-rarefaction.qzv"
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/05a_rarefaction.sh
```

Open `alpha-rarefaction.qzv` and look for the depth at which rarefaction
curves plateau. Your rarefaction depth should be at or beyond this plateau —
meaning deeper sequencing would not substantially increase the observed
diversity. If endosphere curves plateau earlier than bulk soil curves,
this is already a biological signal consistent with H1.

---

## Step 2: Core Diversity Metrics

The `core-metrics-phylogenetic` pipeline rarefies your table, then
computes alpha and beta diversity metrics in a single command, outputting
a comprehensive set of artifacts.

> **Replace `SAMPLING_DEPTH` (line 192) below with your chosen value from Step 1.**

```bash
cat > ${MICROBIOME}/scripts/05b_core_diversity.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=core_diversity
#SBATCH --output=logs/core_diversity_%j.out
#SBATCH --error=logs/core_diversity_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
SHARED_Q2=/fs/scratch/PAS3260/Microbiome
Q2_CONTAINER=${SHARED_Q2}/Containers/qiime2.sif

# ---- Set your rarefaction depth here ----
SAMPLING_DEPTH=10000    # ← replace with your chosen depth

echo "=== Core diversity metrics ==="
echo "Sampling depth: ${SAMPLING_DEPTH}"
echo "Started: $(date)"

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny  /data/02_qiime2/phylogeny/rooted-tree.qza \
    --i-table      /data/02_qiime2/taxonomy/table-no-contam.qza \
    --p-sampling-depth ${SAMPLING_DEPTH} \
    --m-metadata-file  /data/02_qiime2/metadata.tsv \
    --p-n-jobs-or-threads 8 \
    --output-dir   /data/02_qiime2/diversity/core-metrics-results

echo ""
echo "=== Core diversity outputs (in core-metrics-results/) ==="
ls ${MICROBIOME}/02_qiime2/diversity/core-metrics-results/
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/05b_core_diversity.sh
squeue -u ${USER}
```

The output directory will contain:

```
core-metrics-results/
├── rarefied_table.qza                   ← rarefied feature table
├── observed_features_vector.qza         ← alpha: Observed ASVs
├── shannon_vector.qza                   ← alpha: Shannon entropy
├── evenness_vector.qza                  ← alpha: Pielou's evenness
├── faith_pd_vector.qza                  ← alpha: Faith's PD
├── bray_curtis_distance_matrix.qza      ← beta: Bray-Curtis
├── jaccard_distance_matrix.qza          ← beta: Jaccard
├── weighted_unifrac_distance_matrix.qza ← beta: weighted UniFrac
├── unweighted_unifrac_distance_matrix.qza ← beta: unweighted UniFrac
├── bray_curtis_pcoa_results.qza         ← PCoA of Bray-Curtis
├── jaccard_pcoa_results.qza             ← PCoA of Jaccard
├── weighted_unifrac_pcoa_results.qza    ← PCoA of weighted UniFrac
├── unweighted_unifrac_pcoa_results.qza  ← PCoA of unweighted UniFrac
└── *_emperor.qzv                        ← interactive 3D PCoA visualizations
```

---

## Step 3: Alpha Diversity Statistical Testing

Test H1 (alpha diversity decreases from soil to endosphere) using
the Kruskal-Wallis test — a non-parametric test appropriate for small
sample sizes (n=4 per group) without assuming normal distributions.

```bash
cat > ${MICROBIOME}/scripts/05c_alpha_stats.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=alpha_stats
#SBATCH --output=logs/alpha_stats_%j.out
#SBATCH --error=logs/alpha_stats_%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
SHARED_Q2=/fs/scratch/PAS3260/Microbiome
Q2_CONTAINER=${SHARED_Q2}/Containers/qiime2.sif

DIVERSITY_RESULTS=${MICROBIOME}/02_qiime2/diversity/core-metrics-results

echo "=== Alpha diversity significance testing (Kruskal-Wallis) ==="
echo "Started: $(date)"

for METRIC in faith_pd shannon observed_features evenness; do
  echo "--- Testing: ${METRIC} ---"
  apptainer exec \
    --bind ${MICROBIOME}:/data \
    --env MPLCONFIGDIR=/tmp \
    ${Q2_CONTAINER} \
    qiime diversity alpha-group-significance \
      --i-alpha-diversity /data/02_qiime2/diversity/core-metrics-results/${METRIC}_vector.qza \
      --m-metadata-file   /data/02_qiime2/metadata.tsv \
      --o-visualization   /data/03_visualizations/alpha-${METRIC}-significance.qzv
done

echo ""
echo "Alpha significance visualizations written to 03_visualizations/"
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/05c_alpha_stats.sh
```

---

## Step 4: Beta Diversity Statistical Testing (PERMANOVA)

Test H2 (compartment explains significant variation in community composition)
using PERMANOVA (Permutational Multivariate Analysis of Variance). This test
determines what fraction of total beta diversity variance is attributable to
the `compartment` grouping factor, and whether this is more than expected by chance.

```bash
cat > ${MICROBIOME}/scripts/05d_beta_stats.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=beta_stats
#SBATCH --output=logs/beta_stats_%j.out
#SBATCH --error=logs/beta_stats_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
SHARED_Q2=/fs/scratch/PAS3260/Microbiome
Q2_CONTAINER=${SHARED_Q2}/Containers/qiime2.sif

echo "=== Beta diversity significance testing (PERMANOVA) ==="
echo "Grouping variable: compartment"
echo "Started: $(date)"

for METRIC in bray_curtis weighted_unifrac unweighted_unifrac jaccard; do
  echo "--- PERMANOVA: ${METRIC} ---"
  apptainer exec \
    --bind ${MICROBIOME}:/data \
    --env MPLCONFIGDIR=/tmp \
    ${Q2_CONTAINER} \
    qiime diversity beta-group-significance \
      --i-distance-matrix /data/02_qiime2/diversity/core-metrics-results/${METRIC}_distance_matrix.qza \
      --m-metadata-file   /data/02_qiime2/metadata.tsv \
      --m-metadata-column compartment \
      --p-method permanova \
      --p-pairwise \
      --p-permutations 999 \
      --o-visualization /data/03_visualizations/beta-${METRIC}-permanova.qzv
done

echo ""
echo "PERMANOVA visualizations written to 03_visualizations/"
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/05d_beta_stats.sh
```

---

## Step 5: Interpreting Your Results

### Interpreting Alpha Diversity Box Plots

Open each `alpha-*-significance.qzv` file. Each shows box plots of the
diversity metric for each compartment, and a table of Kruskal-Wallis
H-statistics and p-values.

Fill in the table below:

| Metric | BulkSoil median | Endosphere median | Kruskal-Wallis p |
| --- | --- | --- | --- |
| Faith's PD | | | |
| Shannon Entropy | | | |
| Observed Features | | | |
| Pielou's Evenness | | | |

**Guiding questions:**

1. Is the decrease from bulk soil to endosphere statistically significant for
   Faith's PD? For Shannon entropy? If they disagree, what does this suggest about
   the nature of the filtering — does host selection preferentially remove
   phylogenetically distinct taxa, or is it purely numerical?

2. Is evenness higher or lower in the endosphere compared to bulk soil? A lower
   evenness in the endosphere would suggest that one or a few taxa dominate the
   endosphere community, consistent with a very tight host filter.

### Interpreting PCoA Plots (Emperor)

Open `bray_curtis_emperor.qzv` in [view.qiime2.org](https://view.qiime2.org).
The Emperor interface displays an interactive 3D ordination. Color samples by
`compartment`.

**What to look for:**

- Do samples from the same compartment cluster together?
- Is there a visible gradient from bulk soil to endosphere along a principal
  coordinate axis?
- How much variance is explained by PC1 and PC2?

Repeat with `unweighted_unifrac_emperor.qzv` and `weighted_unifrac_emperor.qzv`.

**Guiding question:** Which beta diversity metric shows the strongest
compartment clustering — Bray-Curtis, unweighted UniFrac, or weighted UniFrac?
What does this imply about whether the difference between compartments is
driven by *which taxa are present* (unweighted) or *how abundant they are*
(weighted)?

### Interpreting PERMANOVA Results

Open each `beta-*-permanova.qzv`. Look at:

- **Pseudo-F statistic:** Larger values indicate stronger compartment signal
- **p-value (permuted):** Should be < 0.05 for a significant effect
- **R² value:** Fraction of total variance explained by compartment
- **Pairwise tests:** Which compartment pairs are most different?

Fill in the summary table:

| Distance metric | R² | p-value | Most different compartment pair |
| --- | --- | --- | --- |
| Bray-Curtis | | | |
| Weighted UniFrac | | | |
| Unweighted UniFrac | | | |

> **Expected result from Edwards et al.:** Compartment is consistently
> the strongest predictor of community composition, explaining 40–60% of
> total beta diversity variance. The bulk soil–endosphere comparison is
> always the most dissimilar pair. This result strongly supports H2.

---

## Expected Files After Module 5

```
${MICROBIOME}/02_qiime2/diversity/
└── core-metrics-results/
    ├── rarefied_table.qza
    ├── *_vector.qza                          ← 4 alpha diversity metrics
    ├── *_distance_matrix.qza                 ← 4 beta diversity distance matrices
    ├── *_pcoa_results.qza                    ← 4 PCoA ordinations
    └── *_emperor.qzv                         ← 4 interactive 3D ordination plots

${MICROBIOME}/03_visualizations/
├── alpha-rarefaction.qzv
├── alpha-faith_pd-significance.qzv
├── alpha-shannon-significance.qzv
├── alpha-observed_features-significance.qzv
├── alpha-evenness-significance.qzv
├── beta-bray_curtis-permanova.qzv
├── beta-weighted_unifrac-permanova.qzv
├── beta-unweighted_unifrac-permanova.qzv
└── beta-jaccard-permanova.qzv
```

---

## Checkpoint: Before You Proceed

* All core-metrics-results artifacts exist and the output directory contains at least 14 `.qza` files and 4 `.qzv` Emperor files
* You have opened and interpreted at least one alpha diversity significance visualization and recorded Kruskal-Wallis p-values
* You have opened the Bray-Curtis Emperor ordination and confirmed that samples cluster by compartment
* You have a PERMANOVA R² and p-value for at least one distance metric
* You can explain the difference between weighted and unweighted UniFrac and predict which would be more sensitive to the endosphere filter
* You have documented your rarefaction depth choice and the biological justification for it

---

*Previous: [Module 4 → Taxonomic Classification](04_taxonomy.md) | Next: [Module 6 → Differential Abundance](06_differential_abundance.md)*