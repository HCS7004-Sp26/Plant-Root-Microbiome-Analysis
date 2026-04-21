# Module 6 — Differential Abundance and Ecological Synthesis

## Identifying Host-Enriched and Host-Excluded Taxa with ANCOM-BC

---

## Background

### From "Who Is There?" to "Who Changed?"

Diversity analyses (Module 5) told us *how much* communities differ across
compartments. This module asks the next logical question: **which specific
taxa are responsible for those differences?** This is the domain of
**differential abundance analysis** — identifying taxa that are statistically
enriched or depleted in one condition relative to another.

### Why Differential Abundance in Microbiome Data Is Hard

Microbiome count data has two properties that make standard approaches
(like a simple t-test on raw counts) statistically invalid:

1. **Compositionality:** Amplicon sequencing gives you *relative* abundances
   (proportions), not absolute counts. If one taxon becomes more abundant,
   all others appear less abundant by arithmetic — even if their absolute
   numbers didn't change. This creates **false negative correlations** between
   taxa that can produce spurious results.

2. **Sparsity:** Most ASVs are present in only a subset of samples, with
   many zero counts. Standard parametric tests assume non-sparse data.

### ANCOM-BC: The Modern Solution

**ANCOM-BC** (Analysis of Composition of Microbiomes with Bias Correction)
addresses compositionality by:

1. Estimating and correcting for sample-specific **sampling fraction bias**
   (i.e., the fact that different samples have different total microbial loads)
2. Modeling log-transformed counts with a linear model
3. Applying multiple testing correction (Holm-Bonferroni by default)

The result is a set of taxa with statistically supported log-fold changes
between groups, controlling for false discovery. ANCOM-BC is currently the
most statistically rigorous and widely recommended method for compositional
differential abundance testing in microbiome studies.

### The Biological Question for This Module

We are testing **H3:** Specific bacterial taxa are significantly enriched
in the endosphere relative to bulk soil, consistent with active host selection
for specific microbial partners. We expect:

- **Enriched in endosphere:** Proteobacteria (especially Betaproteobacteria
  and Gammaproteobacteria), Actinobacteria
- **Depleted in endosphere:** Acidobacteria, Verrucomicrobia (these
  are competitive in bulk soil but cannot colonize root interiors)

These predictions come directly from Edwards et al. 2015 and are
well-established in the plant microbiome literature.

---

## Learning Objectives

* Run ANCOM-BC differential abundance analysis comparing endosphere vs. bulk soil
* Interpret log-fold change values and ANCOM-BC q-values
* Visualize differentially abundant taxa using bar plots and volcano-style summaries
* Synthesize results from all modules to evaluate the three core hypotheses
* Articulate limitations of the analysis and directions for future work

---

## Step 1: Run ANCOM-BC

We run two comparisons: (1) all compartments simultaneously (omnibus test
for any compartment effect) and (2) a targeted endosphere vs. bulk soil
pairwise comparison. The `compartment` column in the metadata is the
formula variable. The default reference level for QIIME2's ANCOM-BC is
the first factor level alphabetically — for `compartment`, this will be
`BulkSoil`.

> **Note:** ANCOM-BC requires non-rarefied count data — use the filtered
> table (`table-no-contam.qza`), **not** the rarefied table. ANCOM-BC
> accounts for differential sampling depth internally through its bias
> correction step.

```bash
cat > ${MICROBIOME}/scripts/06a_ancombc.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=ancombc
#SBATCH --output=logs/ancombc_%j.out
#SBATCH --error=logs/ancombc_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
SHARED_Q2=/fs/scratch/PAS3260/Microbiome
Q2_CONTAINER=${SHARED_Q2}/Containers/qiime2.sif

echo "=== ANCOM-BC: Differential abundance across compartments ==="
echo "Formula: compartment (reference level: BulkSoil)"
echo "Input: table-no-contam.qza (non-rarefied counts)"
echo "Started: $(date)"

# ---- ANCOM-BC: omnibus test across all compartments ----
apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime composition ancombc \
    --i-table /data/02_qiime2/taxonomy/table-no-contam.qza \
    --m-metadata-file /data/02_qiime2/metadata.tsv \
    --p-formula 'compartment' \
    --p-reference-levels 'compartment::BulkSoil' \
    --o-differentials /data/02_qiime2/differential_abundance/ancombc-differentials.qza \
    --verbose

echo ""
echo "=== Tabulate ANCOM-BC results ==="

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime composition tabulate \
    --i-data /data/02_qiime2/differential_abundance/ancombc-differentials.qza \
    --o-visualization /data/03_visualizations/ancombc-differentials.qzv

echo ""
echo "=== ANCOM-BC bar chart: taxonomy-aggregated differentials ==="

# Add taxonomy to the differentials for interpretable bar chart
apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime composition da-barplot \
    --i-data /data/02_qiime2/differential_abundance/ancombc-differentials.qza \
    --p-significance-threshold 0.05 \
    --p-level-delimiter ';' \
    --o-visualization /data/03_visualizations/ancombc-da-barplot.qzv

echo ""
echo "ANCOM-BC complete."
echo "Outputs: ancombc-differentials.qza, ancombc-differentials.qzv, ancombc-da-barplot.qzv"
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/06a_ancombc.sh
squeue -u ${USER}
```

---

## Step 2: Visualize and Interpret ANCOM-BC Results

Transfer the visualizations to your local computer:

```bash
# Run on your LOCAL machine
scp <your_username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/Microbiome/03_visualizations/ancombc-differentials.qzv ~/Desktop/
scp <your_username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/Microbiome/03_visualizations/ancombc-da-barplot.qzv ~/Desktop/
```

### Interpreting `ancombc-differentials.qzv`

This table shows, for each ASV and each compartment coefficient
(Rhizosphere, Rhizoplane, Endosphere — all relative to BulkSoil):

| Column | Meaning |
| --- | --- |
| `lfc` | Log-fold change vs. BulkSoil (positive = enriched in that compartment) |
| `se` | Standard error of the log-fold change estimate |
| `W` | Wald statistic |
| `p_val` | Raw p-value |
| `q_val` | FDR-corrected q-value (Holm-Bonferroni) |
| `passed_ss` | Whether sensitivity score threshold was passed |

**Focus on the Endosphere coefficient columns.** Look for:
- ASVs with large positive `lfc` (endosphere-enriched)
- ASVs with large negative `lfc` (endosphere-depleted)
- Both with `q_val` < 0.05

### Interpreting `ancombc-da-barplot.qzv`

The differential abundance bar chart shows which taxa are significantly
enriched (positive bars, shown in orange/red) or depleted (negative bars,
shown in blue/purple) in each compartment relative to bulk soil, with
features sorted by effect size. This is the most visually intuitive
summary of H3.

**Guided analysis:**

Fill in the table with the top 5 enriched and depleted taxa in the
endosphere (use the highest-confidence features with q_val < 0.05):

| Direction | Taxonomy | log-fold change | q-value |
| --- | --- | --- | --- |
| Enriched | | | |
| Enriched | | | |
| Enriched | | | |
| Depleted | | | |
| Depleted | | | |

**Guiding questions:**

1. At what phylum or class level do the endosphere-enriched taxa cluster?
   Is this consistent with the Proteobacteria enrichment reported by Edwards et al.?

2. Are any of the endosphere-depleted taxa affiliated with Acidobacteria or
   Verrucomicrobia? These phyla are known to be dominant in bulk agricultural
   soils but largely excluded from root interiors.

3. Note that ANCOM-BC reports log-fold changes on the *compositional* scale —
   an enriched taxon in the endosphere could be enriched because it is genuinely
   more abundant, or because other taxa are depleted. How would you design a
   follow-up experiment to distinguish these interpretations?

---

## Step 3: Cross-Module Synthesis

You have now completed the full analysis pipeline. Use the table below
to compile evidence for and against each hypothesis:

### H1: Alpha Diversity Decreases from Bulk Soil to Endosphere

| Evidence | Supports H1? | Notes |
| --- | --- | --- |
| Faith's PD Kruskal-Wallis result | | |
| Shannon entropy trend across compartments | | |
| Rarefaction curves plateau earlier in endosphere | | |
| **Overall verdict:** | | |

### H2: Beta Diversity Differs Significantly Among Compartments

| Evidence | Supports H2? | Notes |
| --- | --- | --- |
| PCoA ordination: samples cluster by compartment | | |
| PERMANOVA R² and p-value | | |
| Endosphere–BulkSoil pairwise PERMANOVA | | |
| **Overall verdict:** | | |

### H3: Specific Taxa Are Enriched in the Endosphere

| Evidence | Supports H3? | Notes |
| --- | --- | --- |
| Significant ANCOM-BC results (q < 0.05) | | |
| Endosphere-enriched taxa: phylum-level pattern | | |
| Endosphere-depleted taxa: phylum-level pattern | | |
| Consistency with taxa bar plot (Module 4) | | |
| **Overall verdict:** | | |

---

## Step 4: Limitations and Future Directions

### Limitations of This Analysis

1. **Marker gene vs. function:** 16S rRNA sequencing tells you *who is there*
   but not *what they are doing*. Two Proteobacteria ASVs classified at the
   same genus may have very different metabolic capabilities. Shotgun
   metagenomics, metatranscriptomics, or isolate culturing would be needed
   to link community membership to function.

2. **Cross-sectional design:** This is a single time point (vegetative stage).
   The dynamics of microbiome assembly during plant development — from seed
   germination through senescence — are invisible in this analysis.

3. **Compositional data constraints:** ANCOM-BC corrects for compositionality,
   but cannot determine whether enrichment reflects *absolute* increases in
   a taxon or *relative* increases due to depletion of competitors.

4. **Sample size:** With n=4 replicates per compartment, statistical power
   is limited. Some biologically real differences may not reach significance.

5. **V4 resolution limit:** The V4 region resolves most taxa to genus level
   but rarely to species. Host-specific selection may operate at the strain
   level, which is invisible to 16S rRNA analysis.

### Future Directions

- **Shotgun metagenomics** of the endosphere samples to profile functional gene
  content and identify metabolic pathways enriched by host selection
- **Comparative analysis across cultivars and conditions:** Does M104 vs.
  Nipponbare vs. IR50 select for different endosphere communities? Does the
  greenhouse community differ from the field-grown community using the same
  soil? Both questions are directly testable with the full Edwards et al.
  dataset (you used only M104 greenhouse samples here)
- **Developmental time series:** How does the endosphere community shift from
  vegetative stage through flowering to grain fill?
- **Microbial isolation and inoculation experiments:** Can endosphere-enriched
  taxa (identified here) promote rice growth when inoculated into sterile plants
  grown in defined conditions?

---

## Expected Files After Module 6

```
${MICROBIOME}/02_qiime2/differential_abundance/
└── ancombc-differentials.qza          ← ANCOM-BC result artifact

${MICROBIOME}/03_visualizations/
├── ancombc-differentials.qzv           ← tabular ANCOM-BC results
└── ancombc-da-barplot.qzv              ← differential abundance bar chart
```

---

## Final Checkpoint

* `ancombc-differentials.qza` and both `.qzv` files exist
* You have opened `ancombc-da-barplot.qzv` and identified at least 3 enriched and 3 depleted taxa in the endosphere
* You have completed the cross-module synthesis table (H1, H2, H3)
* You can articulate at least two biological limitations of the 16S rRNA amplicon approach
* You can describe one follow-up experiment that would extend the findings from this tutorial

---

## Complete Tutorial File Tree

```
${MICROBIOME}/
├── 00_raw_reads/          ← 32 FASTQ.gz files
├── 01_qc/                 ← FastQC + MultiQC reports
├── 02_qiime2/
│   ├── metadata.tsv
│   ├── import/            ← pe-manifest.tsv, demux-paired-end.qza
│   ├── denoising/         ← table.qza, rep-seqs.qza, denoising-stats.qza
│   ├── taxonomy/          ← taxonomy.qza, table-no-contam.qza, rep-seqs-no-contam.qza
│   ├── phylogeny/         ← aligned seqs, tree artifacts
│   ├── diversity/         ← core-metrics-results/
│   └── differential_abundance/ ← ancombc-differentials.qza
├── 03_visualizations/     ← all .qzv files for view.qiime2.org
├── scripts/               ← 06a–06a SLURM scripts
└── logs/                  ← SLURM stdout/stderr
```

---

*Previous: [Module 5 → Diversity Analysis](05_diversity.md)*

---

*Tutorial complete. Data source: Edwards et al. (2015) PNAS 112:E911–E920,
BioProject PRJNA255789. Pipeline: QIIME2 2024.5 amplicon distribution.*