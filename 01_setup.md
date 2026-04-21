# Module 1 — Environment Setup

## Directory Structure, Containers, and Shared Resources

---

## Background

The QIIME2 artifact system keeps track of every analysis step through
`.qza` (QIIME2 Artifact) and `.qzv` (QIIME2 Visualization) files. Each
artifact is a zip archive that contains the data, its provenance history,
and its semantic type — a record of exactly how it was produced. This means
your analysis is self-documenting by design. A well-organized directory
structure complements this system by keeping raw data, intermediate
artifacts, and visualizations clearly separated.

> **Relationship to the previous tutorial:** This is a standalone tutorial.
> You do not need any files from the Genome Analytics or RNA-seq tutorials.
> All data are downloaded fresh from NCBI SRA in Module 2.

---

## Learning Objectives

* Set the `user_name` variable that customizes every path in this tutorial
* Create the microbiome project directory tree
* Pull utility containers (sra-tools, FastQC, MultiQC) to your personal
  containers directory
* Verify access to the shared QIIME2 container and SILVA classifier
* Set persistent environment variables for the project

---

## Step 1: Set Your Username Variable

> **Do this first, every session.** All paths in this tutorial use
> `${user_name}` to point to your personal working directory.

```bash
# Replace "myusername" with your actual OSC username
export user_name=Jonathan

# Verify
echo "user_name    = ${user_name}"
echo "Working root = /fs/scratch/PAS3260/${user_name}/Microbiome"
```

---

## Step 2: Define Path Variables

```bash
# Personal working directory and containers
export MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
export CONTAINERS=/fs/scratch/PAS3260/${user_name}/Microbiome/containers

# Shared QIIME2 resources (pre-installed, read-only)
export SHARED_Q2=/fs/scratch/PAS3260/Microbiome
export Q2_CONTAINER=${SHARED_Q2}/Containers/qiime2.sif
export SILVA_CLASSIFIER=${SHARED_Q2}/Classifiers/silva-138-99-nb-classifier.qza 
export TUTORIAL_META=${SHARED_Q2}/tutorial_metadata
```

Confirm the shared resources are accessible:

```bash
ls -lh ${SHARED_Q2}/
# Expected output:
# Classifiers
# Containers
# tutorial_metadata/

ls -lh ${TUTORIAL_META}/
# Expected output:
# accessions_tutorial_subset.txt
# metadata.tsv
```

---

## Step 3: Create the Project Directory Structure

```bash
mkdir -p \
  ${MICROBIOME}/00_raw_reads \
  ${MICROBIOME}/01_qc \
  ${MICROBIOME}/02_qiime2/import \
  ${MICROBIOME}/02_qiime2/denoising \
  ${MICROBIOME}/02_qiime2/taxonomy \
  ${MICROBIOME}/02_qiime2/phylogeny \
  ${MICROBIOME}/02_qiime2/diversity \
  ${MICROBIOME}/02_qiime2/differential_abundance \
  ${MICROBIOME}/03_visualizations \
  ${MICROBIOME}/logs \
  ${MICROBIOME}/scripts \
  ${CONTAINERS}

tree -d ${MICROBIOME}
```

### Annotated directory tree

```
/fs/scratch/PAS3260/${user_name}/Microbiome/
├── 00_raw_reads/            ← downloaded FASTQ files from SRA (16 samples)
├── 01_qc/                   ← FastQC reports and MultiQC summary
├── 02_qiime2/
│   ├── import/              ← demux artifact (.qza) + PE manifest file
│   ├── denoising/           ← DADA2 output: table, rep-seqs, denoising stats
│   ├── taxonomy/            ← taxonomy classifications, filtered table
│   ├── phylogeny/           ← aligned rep-seqs, masked alignment, tree artifacts
│   ├── diversity/           ← rarefaction, alpha/beta diversity artifacts and stats
│   └── differential_abundance/ ← ANCOM-BC output
├── 03_visualizations/       ← exported .qzv files for QIIME2 View
├── containers/              ← Apptainer .sif images (utility tools)
├── logs/                    ← SLURM stdout/stderr for every job
└── scripts/                 ← all SLURM batch scripts
```

> **Note:** `.qza` artifacts are the currency of QIIME2 — they carry both
> data and provenance. The `03_visualizations/` folder is for `.qzv` files
> that you will open in your browser at [view.qiime2.org](https://view.qiime2.org).

---

## Step 4: Persist Environment Variables

Write all path variables to `~/.bashrc` so they survive log-out.

```bash
cat >> ~/.bashrc << EOF

# ---- HCS 7004 Microbiome tutorial ----
export user_name=${user_name}
export MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
export CONTAINERS=/fs/scratch/PAS3260/${user_name}/Microbiome/containers
export SHARED_Q2=/fs/scratch/PAS3260/Microbiome
export Q2_CONTAINER=\${SHARED_Q2}/Containers/qiime2.sif
export SILVA_CLASSIFIER=\${SHARED_Q2}/Classifiers/silva-138-99-nb-classifier.qza
export TUTORIAL_META=\${SHARED_Q2}/tutorial_metadata
EOF

source ~/.bashrc

# Verify
echo "user_name          = ${user_name}"
echo "MICROBIOME         = ${MICROBIOME}"
echo "Q2_CONTAINER       = ${Q2_CONTAINER}"
echo "SILVA_CLASSIFIER   = ${SILVA_CLASSIFIER}"
```

---

## Step 5: Pull Utility Containers

The QIIME2 container is a shared class resource — **do not pull it**.
You only need to pull sra-tools, FastQC, and MultiQC.

```bash
cat > ${MICROBIOME}/scripts/01_pull_containers.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=pull_containers
#SBATCH --output=logs/pull_containers_%j.out
#SBATCH --error=logs/pull_containers_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
CONTAINERS=${MICROBIOME}/containers

echo "=== Pulling utility containers ==="
echo "Started: $(date)"

# ---- sra-tools (SRA download) ----
apptainer pull \
  ${CONTAINERS}/sratools.sif \
  oras://community.wave.seqera.io/library/sra-tools:3.2.1--846898724ee33c64

# ---- FastQC (per-read quality assessment) ----
apptainer pull \
  ${CONTAINERS}/fastqc.sif \
  oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960

# ---- MultiQC (aggregate QC report) ----
apptainer pull \
  ${CONTAINERS}/multiqc.sif \
  oras://community.wave.seqera.io/library/multiqc:1.34--4fc8657c816047c0

echo ""
echo "=== All utility containers pulled ==="
echo "Finished: $(date)"
ls -lh ${CONTAINERS}/
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/01_pull_containers.sh
squeue -u ${USER}
```

> **Note:** Container tags change as new builds are released. If any
> `oras://` URI returns a 404, visit <https://seqera.io/containers/> and
> search for the tool name to find the current tag.

---

## Step 6: Verify All Containers

Once the pull job completes (~20–30 min), verify each container:

```bash
echo "=== Shared QIIME2 container ==="
apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime info 2>&1 | head -6

echo ""
echo "=== Utility containers ==="
echo -n "sra-tools:  "
apptainer exec ${CONTAINERS}/sratools.sif \
  fastq-dump --version 2>&1 | grep "fastq-dump"

echo -n "FastQC:     "
apptainer exec ${CONTAINERS}/fastqc.sif \
  fastqc --version 2>&1

echo -n "MultiQC:    "
apptainer exec ${CONTAINERS}/multiqc.sif \
  multiqc --version 2>&1

echo ""
echo "Verification complete."
```

Expected output from `qiime info`:
```
System versions
Python version: 3.9.x
QIIME 2 release: 2026.1
QIIME 2 version: 2026.1.x
...
```

---

## How QIIME2 Is Run in This Tutorial

Every QIIME2 command uses the following Apptainer execution pattern,
binding your project directory into the container and setting a writable
temp directory for matplotlib configuration:

```bash
apptainer exec \
  --bind ${MICROBIOME}:/data \
  --bind ${SHARED_Q2}:/shared \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime <plugin> <method> \
    --i-... /data/... \
    --o-... /data/... \
    ...
```

| Bind mount | Purpose |
| --- | --- |
| `${MICROBIOME}:/data` | Your project data, accessible as `/data` inside container |
| `${SHARED_Q2}:/shared` | Shared SILVA classifier and metadata, as `/shared` |
| `MPLCONFIGDIR=/tmp` | Prevents matplotlib config errors in read-only home dirs |

> **QIIME2 path convention:** All input/output paths in QIIME2 commands
> inside the container use `/data/` prefixes (e.g., `/data/02_qiime2/import/`).
> On the host, these map to `${MICROBIOME}/02_qiime2/import/`.

---

## SLURM Header Template for This Tutorial

Every script uses the following header. Always `cd ${MICROBIOME}` before
calling `sbatch` so SLURM logs land in `${MICROBIOME}/logs/`.

```bash
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=descriptive_name
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=HH:MM:SS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=N
#SBATCH --mem=XG

set -euo pipefail

# ---- Path variables (hardcoded — do not use $user_name inside scripts) ----
user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
CONTAINERS=${MICROBIOME}/containers
SHARED_Q2=/fs/scratch/PAS3260/Microbiome
Q2_CONTAINER=${SHARED_Q2}/Containers/qiime2.sif
SILVA_CLASSIFIER=${SHARED_Q2}/Classifiers/silva-138-99-nb-classifier.qza
TUTORIAL_META=${SHARED_Q2}/tutorial_metadata
```

---

## Checkpoint: Before You Proceed

* `echo ${user_name}` prints your OSC username (not empty)
* `echo ${MICROBIOME}` prints `/fs/scratch/PAS3260/<your_username>/Microbiome`
* `tree -d ${MICROBIOME}` shows the full directory structure
* `ls ${SHARED_Q2}/` lists `Classifiers/`, `Containers/`, and `tutorial_metadata/`
* `ls ${TUTORIAL_META}/` lists `accessions_tutorial_subset.txt` and `metadata.tsv`
* All 3 utility container `.sif` files exist in `${CONTAINERS}/`
* QIIME2 container prints version string without errors (Step 6)
* All utility containers print version strings (Step 6)

---

*Previous: [00 README / Overview](README.md) | Next: [Module 2 → Data Download and QC](02_data_and_qc.md)*