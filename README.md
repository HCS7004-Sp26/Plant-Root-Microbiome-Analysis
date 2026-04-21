# HCS 7004 — Genome Analytics

## Tutorial: Plant Root Microbiome Analysis

### *Oryza sativa* (Rice) Rhizosphere and Endosphere Communities — From Raw 16S Reads to Ecological Inference

---

## Overview

In the preceding tutorial series you assembled, predicted, and functionally
annotated a eukaryotic genome using short- and long-read sequencing data.
That work was fundamentally about a *single organism*. This tutorial expands
your analytical scope to *communities*: the thousands of microbial species
that colonize plant root systems and collectively shape plant health,
nutrient acquisition, and stress tolerance.

We use data from **Edwards et al. (2015)** — *"Structure, variation, and
assembly of the root-associated microbiomes of rice"* (*PNAS* 112:E911–E920,
[doi:10.1073/pnas.1414592112](https://doi.org/10.1073/pnas.1414592112)) —
one of the most cited studies in the plant microbiome field. This paper
characterizes the bacterial communities associated with four spatially
distinct root compartments of field-grown rice and addresses a central
ecological hypothesis: **that the plant root progressively filters and
enriches specific bacterial taxa from the surrounding soil, imposing
increasing selective pressure as communities move from bulk soil into
the root interior**.

All 16S rRNA amplicon data (Illumina, V4 region, primers 515F/806R) are
publicly available under **BioProject PRJNA238564**. The analysis pipeline
runs entirely through **QIIME2** in an Apptainer container on the Ohio
Supercomputer Center (OSC).

> **Working directory:** Each student works in their own directory:
> `/fs/scratch/PAS3260/${user_name}/Microbiome`
> Replace `${user_name}` with your OSC username in every command and script.
> The first thing you do in Module 1 is set this variable.

---

## The Central Scientific Questions

This tutorial is organized around three testable hypotheses drawn directly
from Edwards et al.:

| Hypothesis | Ecological concept | Analysis module |
| --- | --- | --- |
| **H1.** Alpha diversity (richness and evenness) decreases progressively from bulk soil through the endosphere | Species filtering / host selection | Module 5 |
| **H2.** Beta diversity (community composition) differs significantly among root compartments | Community assembly / niche differentiation | Module 5 |
| **H3.** Specific bacterial taxa are significantly enriched or depleted in the endosphere relative to bulk soil | Differential abundance / host enrichment | Module 6 |

These are not retrospective questions — you will test them with real data
using the same statistical framework as the original paper.

---

## Tutorial Map

| Module | File | Topic | Est. time |
| --- | --- | --- | --- |
| 0 | `README.md` | Introduction, biological context, pipeline overview | — |
| 1 | `01_setup.md` | **Directory structure, containers, and shared resources** | 30 min |
| 2 | `02_data_and_qc.md` | SRA download, read QC, and QIIME2 import | 45 min |
| 3 | `03_dada2_denoising.md` | DADA2 denoising, feature table, and representative sequences | 60 min |
| 4 | `04_taxonomy.md` | Taxonomic classification with SILVA and community visualization | 45 min |
| 5 | `05_diversity.md` | Alpha and beta diversity, UniFrac, PERMANOVA | 60 min |
| 6 | `06_differential_abundance.md` | ANCOM-BC differential abundance and ecological synthesis | 60 min |

> **Total active work:** ~5–6 hours + SLURM queue time.
> All bioinformatics tools run through **Apptainer containers**.
> The QIIME2 container and the SILVA classifier are pre-installed in a
> shared class directory. Students pull the remaining utility containers
> in Module 1.

---

## Learning Objectives

By the end of this tutorial you will be able to:

1. **Explain** the amplicon sequencing (16S rRNA) workflow and how it differs from whole-genome or transcriptomic sequencing
2. **Execute** the QIIME2 CLI pipeline end-to-end on a real plant microbiome dataset on OSC
3. **Define and calculate** alpha diversity metrics (Shannon entropy, Faith's phylogenetic diversity, Observed Features) and interpret their ecological meaning
4. **Construct and interpret** beta diversity ordinations (PCoA of Bray-Curtis, weighted and unweighted UniFrac distances) and apply PERMANOVA to test group differences
5. **Identify** differentially abundant taxa between root compartments using ANCOM-BC and interpret results in the context of host selection
6. **Critically evaluate** the compartment-gradient model of microbiome assembly proposed by Edwards et al. using your own analysis of their data

---

## Biological and Ecological Background

### The Rice Root Microbiome: A Spatial Gradient

When a rice plant grows in soil, it does not passively accept whatever
microbes are present. It actively shapes the surrounding microbial
community through root exudates, mucilage secretion, and cell wall
chemistry. Edwards et al. conceptualize this as a **series of filters**:

```
                Billions of soil bacteria
                    ▼  (first filter)
           ┌─────────────────────────┐
           │       Rhizosphere       │  ← soil immediately surrounding roots
           │   soil enriched by      │    ~1–2 mm zone; shaped by exudates
           │   root exudates         │
           └─────────────┬───────────┘
                         ▼  (second filter)
           ┌─────────────────────────┐
           │       Rhizoplane        │  ← root surface
           │   root surface biofilm  │    attached to epidermal cells
           └─────────────┬───────────┘
                         ▼  (third filter: tightest)
           ┌─────────────────────────┐
           │       Endosphere        │  ← root interior
           │   inside root tissue    │    inside living plant cells
           └─────────────────────────┘
```

Each filter is selective: the community that enters each compartment is a
*subset* of the community in the compartment above it, enriched for taxa
that can tolerate or exploit that particular chemical and physical
environment. By the time you reach the endosphere, community diversity is
sharply reduced and the remaining taxa are highly consistent across
replicate plants — a signature of *deterministic host selection* rather
than random dispersal.

### The QIIME2 Amplicon Pipeline

The 16S rRNA gene contains nine hypervariable regions (V1–V9). Edwards et al.
amplified the **V4 region** using primers 515F/806R, a standard choice that
balances taxonomic resolution with Illumina read length. QIIME2 processes
the resulting paired-end reads through the following conceptual stages:

```
Raw paired-end FASTQ reads (515F/806R V4 amplicons)
             │
             ▼
     Quality control (FastQC / MultiQC)
             │
             ▼
     Import into QIIME2 artifact system (.qza)
             │
             ▼
     DADA2 denoising
     (error correction → exact amplicon sequence variants, ASVs)
             │ feature table + representative sequences
             ▼
     Taxonomic classification
     (SILVA 138 naïve Bayes classifier)
             │ taxonomy assignments
             ▼
     Contamination filtering
     (remove chloroplasts, mitochondria)
             │ clean feature table
             ▼
     Phylogenetic placement
     (MAFFT alignment → FastTree → rooted tree)
             │
      ┌──────┴──────┐
      ▼             ▼
  Alpha         Beta diversity
  diversity     (Bray-Curtis, weighted
  (Shannon,      and unweighted UniFrac,
  Faith's PD,    PCoA, PERMANOVA)
  Observed Feat.)
      └──────┬──────┘
             ▼
   Differential abundance
   (ANCOM-BC: compartment comparisons)
```

### ASVs vs. OTUs: Why It Matters

Earlier microbiome studies used **operational taxonomic units (OTUs)** —
sequences clustered at 97% identity. QIIME2's DADA2 plugin instead
produces **amplicon sequence variants (ASVs)**, which are exact biological
sequences corrected for sequencing error. ASVs have higher resolution,
are reproducible across studies, and allow you to track individual
sequence variants rather than fuzzy clusters.

> **Note on the tutorial dataset:** Edwards et al. (2015) originally
> used OTU-based analysis. We reanalyze their raw reads with the modern
> ASV approach via DADA2, which produces a more refined and reproducible
> feature table. Minor differences from the published figures are
> expected and are themselves a valuable teaching point about methodological
> evolution.

---

## Tutorial Dataset: Sample Subset

BioProject PRJNA238564 contains ~600 samples spanning three cultivars,
four compartments, three developmental stages, and two field sites.
For this tutorial we use a curated **pedagogical subset** that preserves
the central scientific question while remaining computationally tractable:

| Factor | Selection | Rationale |
| --- | --- | --- |
| Cultivar | Nipponbare | *O. sativa* japonica reference cultivar |
| Field site | Sacramento | Uniform soil background |
| Developmental stage | Vegetative | Single time point; simplest design |
| Compartments | Bulk soil, Rhizosphere, Rhizoplane, Endosphere | Full spatial gradient |
| Replicates | 4 per compartment | Sufficient power for PERMANOVA |

**Total: 16 samples** — enough to fully test all three hypotheses, manageable
in a single SLURM session per module.

The curated accession list and metadata table are pre-installed in the
shared class directory. You will use them starting in Module 2.

---

## Shared Resources for This Tutorial

The QIIME2 container, the SILVA 138 classifier, and the tutorial metadata
are pre-installed in a shared class directory and do **not** need to be
downloaded by individual students:

```
/fs/scratch/PAS3260/Team_Project/Containers/QIIME2/
├── qiime2.sif    ← QIIME2 2024.5 Apptainer container
├── silva-138-99-515-806-nb-classifier.qza   ← pre-trained V4 classifier
└── tutorial_metadata/
    ├── accessions_tutorial_subset.txt       ← 16 SRR accessions
    └── metadata.tsv                         ← sample metadata (QIIME2-formatted)
```

---

## Pipeline Tool Summary

| Tool | Version | Purpose | Source |
| --- | --- | --- | --- |
| QIIME2 (amplicon distribution) | 2024.5 | Full amplicon analysis pipeline | quay.io/qiime2/amplicon |
| sra-tools | 3.2.1 | Download reads from NCBI SRA | Seqera community |
| FastQC | 0.12.1 | Per-read quality assessment | Seqera community |
| MultiQC | 1.25.1 | Aggregate QC report | Seqera community |

---

## Pre-Tutorial Discussion Questions

Before starting Module 1, consider the following with your group:

1. The 16S rRNA gene is a *marker gene* — it does not tell you what a
   bacterium *does*, only what it *is* (taxonomically). What are the
   practical advantages and limitations of this approach compared to
   shotgun metagenomics for studying the rice root microbiome?

2. The Edwards et al. design is observational, not experimental — plants
   were grown in field conditions, not in sterile pots with defined
   microbiomes. How does this affect your ability to draw causal inferences
   about host selection? What experimental design would be needed to
   establish causality?

3. If host selection is the dominant force structuring the endosphere
   community, you would predict that endosphere communities are *more
   similar to each other* (lower beta diversity within the endosphere)
   than bulk soil communities. Why? What ecological null model does this
   contrast with?

---

*Proceed to:* **[Module 1 → Setup](01_setup.md)**
