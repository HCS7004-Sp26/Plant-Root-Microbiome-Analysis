# Module 2 — Data Download, Quality Control, and QIIME2 Import

## Obtaining the Edwards et al. (2015) Tutorial Subset from NCBI SRA

---

## Background

The Edwards et al. dataset (BioProject **PRJNA238564**) contains ~600 samples
spanning three rice cultivars, four root compartments, three developmental
stages, and two California field sites. For this tutorial we work with a
curated pedagogical subset: **16 samples** representing the Nipponbare
cultivar at the vegetative stage from the Sacramento field site, with four
replicates per compartment (bulk soil, rhizosphere, rhizoplane, endosphere).
This design isolates the **compartment gradient** — the central scientific
variable — while keeping computational demands feasible.

### Why This Subset?

| Full dataset | Tutorial subset | Reason |
| --- | --- | --- |
| ~600 samples | 16 samples | Computationally tractable on OSC |
| 3 cultivars | Nipponbare only | Focus on compartment gradient hypothesis |
| 2 field sites | Sacramento only | Remove site-to-site variation as a confounder |
| 3 time points | Vegetative only | Single time point; cleanest design |
| 4 compartments | All 4 | Preserve the full spatial gradient |

### About the Metadata File

QIIME2 requires a **sample metadata file** — a tab-separated table where
the first column is `sample-id` (matching your read file names) and
subsequent columns are sample attributes used in statistical tests and
visualizations. The metadata file for this tutorial is pre-installed in
the shared class directory:

```
${TUTORIAL_META}/metadata.tsv
```

Inspect it now to understand the experimental design:

```bash
cat ${TUTORIAL_META}/metadata.tsv
```

You should see columns including `sample-id`, `compartment`, `cultivar`,
`field-site`, `developmental-stage`, and `replicate`. The `compartment`
column (values: `BulkSoil`, `Rhizosphere`, `Rhizoplane`, `Endosphere`) is
the primary grouping variable for all downstream analyses.

### About the Accession File

```
${TUTORIAL_META}/accessions_tutorial_subset.txt
```

This file lists 16 SRR accessions — one per line — corresponding to the
16 samples in the tutorial subset. You will use it to drive the SRA
download in Step 2.

---

## Learning Objectives

* Download 16S rRNA paired-end reads from NCBI SRA using sra-tools
* Assess read quality with FastQC and aggregate results with MultiQC
* Interpret amplicon QC metrics (sequence length distribution, per-base quality, adapter content)
* Construct a QIIME2 paired-end manifest file
* Import reads into the QIIME2 artifact system
* Summarize and visualize the demultiplexed artifact to guide DADA2 truncation

---

## Step 1: Copy the Metadata and Accession Files to Your Working Directory

```bash
# Copy shared resources to your project directory
cp ${TUTORIAL_META}/metadata.tsv         ${MICROBIOME}/02_qiime2/
cp ${TUTORIAL_META}/accessions_tutorial_subset.txt  ${MICROBIOME}/

# Inspect both files
echo "=== metadata.tsv (first 5 lines) ==="
head -5 ${MICROBIOME}/02_qiime2/metadata.tsv

echo ""
echo "=== accession list ==="
cat ${MICROBIOME}/accessions_tutorial_subset.txt
echo "Total accessions: $(wc -l < ${MICROBIOME}/accessions_tutorial_subset.txt)"
```

---

## Step 2: Download Reads from NCBI SRA

We use `fasterq-dump` (part of sra-tools) to download paired-end FASTQ
files for each accession. The `--split-files` flag produces separate `_1`
and `_2` files for forward and reverse reads.

> **Disk note:** Each sample is ~50–150 MB compressed. 16 samples ≈ 2–3 GB total.
> The download runs on a compute node via SLURM to avoid timeout issues on
> the login node.

```bash
cat > ${MICROBIOME}/scripts/02a_download_sra.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=sra_download
#SBATCH --output=logs/sra_download_%j.out
#SBATCH --error=logs/sra_download_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
CONTAINERS=${MICROBIOME}/containers
ACCESSIONS=${MICROBIOME}/accessions_tutorial_subset.txt
OUTDIR=${MICROBIOME}/00_raw_reads

echo "=== SRA Download: Edwards et al. 2015 subset ==="
echo "Started: $(date)"
echo "Accessions to download: $(wc -l < ${ACCESSIONS})"
echo ""

while IFS= read -r accession; do
  if [[ -z "${accession}" ]]; then
    continue
  fi

  echo "--- Downloading ${accession} ---"

  apptainer exec \
    --bind ${MICROBIOME}:/data \
    ${CONTAINERS}/sratools_3.2.1.sif \
    fasterq-dump \
      --outdir /data/00_raw_reads \
      --temp /tmp \
      --split-files \
      --threads 8 \
      --progress \
      "${accession}"

  echo "  Compressing ${accession}..."
  gzip -f ${OUTDIR}/${accession}_1.fastq
  gzip -f ${OUTDIR}/${accession}_2.fastq

  echo "  Done: ${accession}"
done < "${ACCESSIONS}"

echo ""
echo "=== Download complete ==="
echo "Files in ${OUTDIR}:"
ls -lh ${OUTDIR}/*.fastq.gz | awk '{print $5, $9}'
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/02a_download_sra.sh
squeue -u ${USER}
```

When the job finishes, verify the downloads:

```bash
ls ${MICROBIOME}/00_raw_reads/*.fastq.gz | wc -l
# Expected: 32  (16 samples × 2 reads each: _1.fastq.gz and _2.fastq.gz)
```

---

## Step 3: Quality Control with FastQC and MultiQC

Before any analysis, assess read quality. For amplicon data, pay attention to:

- **Per-base sequence quality:** Should be Q30+ for most positions; expect
  a drop at the 3′ end of R2 (reverse reads)
- **Sequence length distribution:** Should be tight around 250 bp for 2×250
  Illumina runs
- **Adapter content:** If primers were not computationally removed by the
  sequencing center, you may see adapter signal
- **Per-sequence GC content:** Amplicon data often shows a sharper GC peak
  than shotgun data because you are sequencing a specific genomic region

```bash
cat > ${MICROBIOME}/scripts/02b_fastqc.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=fastqc_multiqc
#SBATCH --output=logs/fastqc_multiqc_%j.out
#SBATCH --error=logs/fastqc_multiqc_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
CONTAINERS=${MICROBIOME}/containers

echo "=== FastQC: per-file quality assessment ==="
echo "Started: $(date)"

apptainer exec \
  --bind ${MICROBIOME}:/data \
  ${CONTAINERS}/fastqc_0.12.1.sif \
  fastqc \
    --outdir /data/01_qc \
    --threads 8 \
    /data/00_raw_reads/*.fastq.gz

echo ""
echo "=== MultiQC: aggregate QC report ==="

apptainer exec \
  --bind ${MICROBIOME}:/data \
  ${CONTAINERS}/multiqc_1.25.1.sif \
  multiqc \
    --outdir /data/01_qc/multiqc \
    --filename multiqc_report.html \
    /data/01_qc/

echo ""
echo "=== QC complete ==="
echo "Reports in ${MICROBIOME}/01_qc/"
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/02b_fastqc.sh
squeue -u ${USER}
```

Transfer the MultiQC report to your local computer:

```bash
# Run this on your LOCAL machine (not on OSC)
scp <your_username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/Microbiome/01_qc/multiqc/multiqc_report.html ~/Desktop/
```

Open `multiqc_report.html` in your browser and note:

1. Mean quality scores per position for R1 and R2
2. Whether R2 quality drops appreciably before position 200
3. Adapter content signal (guides whether trimming is needed)

> **For Edwards et al. data:** These reads were sequenced on Illumina HiSeq
> with 2×150 bp chemistry. Expect good R1 quality throughout, with some R2
> quality decline after ~130 bp. This will directly inform your DADA2
> truncation parameters in Module 3.

---

## Step 4: Build the QIIME2 Paired-End Manifest File

QIIME2 imports reads via a **manifest file** — a tab-separated table
that maps each sample ID to the absolute paths of its forward and reverse
reads. This is the bridge between your FASTQ files on disk and the
QIIME2 artifact system.

```bash
# Build the manifest from the metadata sample IDs and FASTQ filenames
# The SRR accessions serve as sample IDs here; they must match the metadata.tsv sample-id column

echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" \
  > ${MICROBIOME}/02_qiime2/import/pe-manifest.tsv

while IFS= read -r accession; do
  if [[ -z "${accession}" ]]; then continue; fi
  echo -e "${accession}\t${MICROBIOME}/00_raw_reads/${accession}_1.fastq.gz\t${MICROBIOME}/00_raw_reads/${accession}_2.fastq.gz" \
    >> ${MICROBIOME}/02_qiime2/import/pe-manifest.tsv
done < ${MICROBIOME}/accessions_tutorial_subset.txt

# Verify the manifest
echo "Manifest rows (including header):"
wc -l ${MICROBIOME}/02_qiime2/import/pe-manifest.tsv

echo ""
echo "First 3 rows:"
head -3 ${MICROBIOME}/02_qiime2/import/pe-manifest.tsv
```

> **Why absolute paths?** QIIME2 stores the manifest inside the `.qza`
> artifact. If you used relative paths, the artifact would be non-portable.
> Always use absolute paths in manifest files.

---

## Step 5: Import Reads into QIIME2

This step wraps your FASTQ files in a QIIME2 artifact with semantic type
`SampleData[PairedEndSequencesWithQuality]`. All subsequent QIIME2
commands operate on `.qza` artifacts — never on raw FASTQ files directly.

```bash
cat > ${MICROBIOME}/scripts/02c_import.sh << 'EOF'
#!/bin/bash
#SBATCH --account=PAS3260
#SBATCH --job-name=qiime2_import
#SBATCH --output=logs/qiime2_import_%j.out
#SBATCH --error=logs/qiime2_import_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G

set -euo pipefail

user_name=Jonathan    # ← replace with your OSC username
MICROBIOME=/fs/scratch/PAS3260/${user_name}/Microbiome
SHARED_Q2=/fs/scratch/PAS3260/Team_Project/Containers/QIIME2
Q2_CONTAINER=${SHARED_Q2}/qiime2_amplicon_2024.10.sif

echo "=== QIIME2 Import: PairedEndFastqManifestPhred33V2 ==="
echo "Started: $(date)"

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path /data/02_qiime2/import/pe-manifest.tsv \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path /data/02_qiime2/import/demux-paired-end.qza

echo ""
echo "=== Generate demux summary visualization ==="

apptainer exec \
  --bind ${MICROBIOME}:/data \
  --env MPLCONFIGDIR=/tmp \
  ${Q2_CONTAINER} \
  qiime demux summarize \
    --i-data /data/02_qiime2/import/demux-paired-end.qza \
    --o-visualization /data/03_visualizations/demux-summary.qzv

echo ""
echo "Import complete."
echo "Artifact:      ${MICROBIOME}/02_qiime2/import/demux-paired-end.qza"
echo "Visualization: ${MICROBIOME}/03_visualizations/demux-summary.qzv"
echo "Finished: $(date)"
EOF

cd ${MICROBIOME}
sbatch ${MICROBIOME}/scripts/02c_import.sh
squeue -u ${USER}
```

---

## Step 6: Inspect the Demux Summary Visualization

Transfer the `.qzv` file to your local computer and open it in
[view.qiime2.org](https://view.qiime2.org):

```bash
# Run on your LOCAL machine
scp <your_username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/Microbiome/03_visualizations/demux-summary.qzv ~/Desktop/
```

Then drag `demux-summary.qzv` into [view.qiime2.org](https://view.qiime2.org).

**What to look for:**

The interactive quality plot shows per-position median quality and interquartile
ranges for R1 (forward) and R2 (reverse) reads across all samples. This is the
key figure for choosing your DADA2 truncation lengths in Module 3.

**Guiding questions:**

1. At what position does R1 quality drop below Q25? This is your R1 truncation candidate.
2. At what position does R2 quality drop below Q25? This is your R2 truncation candidate.
3. Are any samples dramatically lower quality than others? Note the sample IDs.
4. How many reads per sample (check the "Per-sample sequence counts" tab)?
   Are any samples substantially under-sequenced?

> **Guideline for Edwards et al. V4 data (2×150 bp):**
> Typical truncation values are `--p-trunc-len-f 150` and `--p-trunc-len-r 130`.
> Your exact values should be determined by the interactive quality plot.
> The V4 amplicon is ~253 bp; truncating to 150 + 130 = 280 bp total gives
> a ~27 bp overlap, sufficient for merging.

---

## Expected Files After Module 2

```
${MICROBIOME}/
├── 00_raw_reads/
│   ├── SRR*.fastq.gz       ← 32 files (16 samples × R1 + R2)
├── 01_qc/
│   ├── *.html              ← FastQC reports (one per FASTQ file)
│   ├── *.zip
│   └── multiqc/
│       └── multiqc_report.html
├── 02_qiime2/
│   └── import/
│       ├── pe-manifest.tsv
│       └── demux-paired-end.qza   ← QIIME2 paired-end sequences artifact
└── 03_visualizations/
    └── demux-summary.qzv          ← Demux quality summary (open in view.qiime2.org)
```

---

## Checkpoint: Before You Proceed

* `ls ${MICROBIOME}/00_raw_reads/*.fastq.gz | wc -l` returns 32
* `ls ${MICROBIOME}/01_qc/multiqc/multiqc_report.html` exists
* `ls ${MICROBIOME}/02_qiime2/import/demux-paired-end.qza` exists
* `ls ${MICROBIOME}/03_visualizations/demux-summary.qzv` exists
* You have opened `demux-summary.qzv` and identified R1 and R2 truncation lengths
* You have noted how many reads per sample and confirmed no sample has fewer than ~5,000 reads
* You understand why absolute paths are required in the manifest file

---

*Previous: [Module 1 → Setup](01_setup.md) | Next: [Module 3 → DADA2 Denoising](03_dada2_denoising.md)*
