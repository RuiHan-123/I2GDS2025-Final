# I2GDS2025-Final-Rui
Biological Question: How does single-cell DNA methylation reveal cell-to-cell heterogeneity and stable epigenetic states that are obscured in bulk measurements?

# Part1: Linux
My Linux workflow is the same as that used by Group 3, this time it includes the bamtools and Picard components that were not described in detail in the group assignment.

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

INDIR="/projects/lu_lab/Rui/qiang_NC_data/11.05/05_bamtools"
OUTDIR="/projects/lu_lab/Rui/qiang_NC_data/11.05/06_picard"
ID_FILE="/projects/lu_lab/Rui/qiang_NC_data/11.05/03_idemp/cell_ids_1.txt"

mkdir -p "${OUTDIR}"

ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ID_FILE}")

BAM_FILE="${INDIR}/${ID}.merged.bam"


SORTED="${OUTDIR}/${ID}.sorted.bam"
DEDUP="${OUTDIR}/${ID}.dedup.bam"
METRICS="${OUTDIR}/${ID}.dedup.metrics.txt"

echo "[$(date)] Start dedup for cell ${ID}" >> 06_picard/picard_all.out

samtools sort -@4 -o "${SORTED}" "${BAM_FILE}"

picard MarkDuplicates \
    I="${SORTED}" \
    O="${DEDUP}" \
    M="${METRICS}" \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT

samtools index "${DEDUP}"

echo "[$(date)] Finished cell ${ID}" >> 06_picard/picard_all.out
```
