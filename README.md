# I2GDS2025-Final-Rui
Biological Question: How does single-cell DNA methylation reveal cell-to-cell heterogeneity and stable epigenetic states that are obscured in bulk measurements?

# Part1: Linux
My Linux workflow is the same as that used by Group 3, this time it includes the bamtools and Picard components that were not described in detail in the group assignment.

Bamtools is a toolkit for processing BAM files, primarily used for basic operations such as viewing, filtering, summarizing, and format conversion of alignment results. Because our data are paired-end sequencing results, after alignment we obtain two BAM files for each cell; therefore, bamtools is used to merge these two BAM files.

```
#!/bin/bash
#SBATCH --job-name=bamtools_merge
#SBATCH --account=chipseq
#SBATCH --partition=normal_q
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=05_bamtools/bamtools_all.out
#SBATCH --error=05_bamtools/bamtools_all.err
#SBATCH --array=1-1000


INDIR="/projects/lu_lab/Rui/qiang_NC_data/11.05/04_bismark_align"
OUTDIR="/projects/lu_lab/Rui/qiang_NC_data/11.05/05_bamtools"
ID_FILE="/projects/lu_lab/Rui/qiang_NC_data/11.05/03_idemp/cell_ids_1.txt"

mkdir -p "${OUTDIR}"

ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ID_FILE}")

R1_BAM="${INDIR}/SRR19391267_1_extracted_val_1.fq.gz_${ID}_bismark_bt2.bam"
R2_BAM="${INDIR}/SRR19391267_2_extracted_val_2.fq.gz_${ID}_bismark_bt2.bam"
MERGED_BAM="${OUTDIR}/${ID}.merged.bam"

if [[ -z "${ID}" ]]; then
  echo "[ERR] Empty cell ID at task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

if [[ ! -f "${R1_BAM}" || ! -f "${R2_BAM}" ]]; then
  echo "[WARN] Missing BAM(s) for cell ${ID}, skipping." >&2
  exit 0
fi

bamtools merge \
  -in "${R1_BAM}" \
  -in "${R2_BAM}" \
  -out "${MERGED_BAM}"

echo "[$(date)] Finished merge for cell ${ID}" >> 05_bamtools/bamtools_all.out
```

# bamtools output
<img width="1849" height="953" alt="14ea2666551debbd4563dc55b1c28efe" src="https://github.com/user-attachments/assets/74240611-bf12-49fe-96c3-bd0ec56f0329" />
This directory contains the merged BAM files produced by bamtools. For each single cell, paired-end alignment BAM files (R1 and R2) were merged into a single BAM file, resulting in one merged BAM per cell.
The variation in file size reflects differences in sequencing depth across individual cells.


# Potential mistakes for bamtools
* Because I need to process 5,000 single cells, with each cell treated as an independent task, and the university Linux cluster allows a maximum of 1,001 concurrent jobs, the cell ID list had to be split into five separate files, each containing 1,000 cell IDs. If the number of parallel tasks I submit exceeds the limit of 1,001, the job will fail to be submitted successfully.
* Input checking before bamtools is necessary because the dataset contains 5,000 single cells, and some samples are inevitably expected to fail; adding this check ensures that only valid BAM files are merged correctly.


Picard is a duplicate-marking tool that identifies and remove PCR duplicates.

```
#!/bin/bash
#SBATCH --job-name=picard_dedup
#SBATCH --account=chipseq
#SBATCH --partition=normal_q
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=06_picard/picard_all.out
#SBATCH --error=06_picard/picard_all.err
#SBATCH --array=1-1000

INDIR="/projects/lu_lab/Rui/qiang_NC_data/11.05/05_bamtools"
OUTDIR="/projects/lu_lab/Rui/qiang_NC_data/11.05/06_picard"
ID_FILE="/projects/lu_lab/Rui/qiang_NC_data/11.05/03_idemp/cell_ids_1.txt"

mkdir -p "${OUTDIR}"

ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ID_FILE}")

if [[ -z "${ID}" ]]; then
  echo "[ERR] Empty cell ID at task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

BAM_FILE="${INDIR}/${ID}.merged.bam"

if [[ ! -f "${BAM_FILE}" ]]; then
  echo "[WARN] Missing merged BAM for cell ${ID}, skipping." >&2
  exit 0
fi

SORTED="${OUTDIR}/${ID}.sorted.bam"
DEDUP="${OUTDIR}/${ID}.dedup.bam"
METRICS="${OUTDIR}/${ID}.dedup.metrics.txt"

samtools sort -@4 -o "${SORTED}" "${BAM_FILE}"

picard MarkDuplicates \
    I="${SORTED}" \
    O="${DEDUP}" \
    M="${METRICS}" \
    REMOVE_DUPLICATES=false \
    VALIDATION_STRINGENCY=SILENT

samtools index "${DEDUP}"

echo "[$(date)] Finished cell ${ID}" >> 06_picard/picard_all.out
```


# Picard output
<img width="1787" height="903" alt="40bde09c848d228f754ba0b24ece5dcd" src="https://github.com/user-attachments/assets/caa7e540-9634-46ce-9c4e-c44bc1e8b2c0" /> <img width="2541" height="320" alt="6c580577b3046e64c8b01355744590b7" src="https://github.com/user-attachments/assets/abe89465-dacb-4dc8-8e53-fb363ce4dc47" />


