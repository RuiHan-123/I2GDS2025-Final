# I2GDS2025-Final-Rui
Biological Question: How does single-cell DNA methylation reveal cell-to-cell heterogeneity and stable epigenetic states that are obscured in bulk measurements?

# Software versions
BAM-level processing was performed using samtools (v1.3.1) and bamtools (v2.5.1). PCR duplicate marking was conducted with Picard (v2.27.5). Data visualization and downstream analyses were performed in R (v4.4.2).


# Part1: Linux
My Linux workflow is the same as that used by Group 3, this time it includes the bamtools and Picard components that were not described in detail in the group assignment.

# Bamtool
Bamtools is a toolkit for processing BAM files, primarily used for basic operations such as viewing, filtering, summarizing, and format conversion of alignment results. Because our data are paired-end sequencing results, after alignment we obtain two BAM files for each cell; therefore, bamtools is used to merge these two BAM files.
First, I specify the job settings for the cluster, including the job name, account, runtime, memory, and the number of array tasks. Then, I define the input directory, output directory, and the file containing cell IDs. For each cell ID, the script selects the corresponding R1 and R2 BAM files and checks whether the cell and input files exist. Finally, bamtools is used to merge the paired-end BAM files for each cell.

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


# Picard
Picard is a duplicate-marking tool that identifies  PCR duplicates.
First, the script sets the cluster job parameters (job name, account, runtime, memory, and array size) and defines the input and output directories along with the cell ID file.
Each array task reads one cell ID and checks whether the corresponding merged BAM file exists.The BAM file is then sorted using samtools, which is required by Picard.
Finally, Picard MarkDuplicates is run to identify PCR duplicates and generate a duplication metrics file.

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
<img width="1787" height="903" alt="40bde09c848d228f754ba0b24ece5dcd" src="https://github.com/user-attachments/assets/caa7e540-9634-46ce-9c4e-c44bc1e8b2c0" />
This directory shows the output files generated by Picard MarkDuplicates for each single cell, including the deduplicated BAM file, its index file, and a duplication metrics report.

<img width="2541" height="320" alt="6c580577b3046e64c8b01355744590b7" src="https://github.com/user-attachments/assets/abe89465-dacb-4dc8-8e53-fb363ce4dc47" />
This is an example of a Picard metrics file showing that 717,235 reads were examined, of which 391,323 did not align to the reference genome, and 54.6% of the reads were identified as PCR duplicates; these values are reasonable for single-cell sequencing data.

# Potential mistakes for Picard
* The output and error logs for each cell should be merged into a single combined file, because running 5,000 cells would otherwise generate 5,000 .out files and 5,000 .err files. This can overwhelm the directory and cause severe lag or even force the file browser to crash when opening the folder.
* After merging paired-end BAM files, samtools sort is used to generate coordinate-sorted BAM files, which are required as input for Picard MarkDuplicates. Because Picard MarkDuplicates requires the input BAM file to be coordinate-sorted, the BAM file must be sorted in advance; otherwise, Picard cannot run properly.
* The -@4 option indicates that samtools uses four threads, which must match the number of CPUs requested with #SBATCH --cpus-per-task=4. If they do not match, it may cause CPU oversubscription, leading to inefficient performance or the job being terminated by the scheduler.


# Part2: R Visualization
In R, I want to create a figure where different cell types are classified into different clusters.
Here, UMAP is used to visualize single-cell methylation profiles in 2D, and k-means (k = 3) groups cells into three clusters based on their positions in the UMAP space. Then, the ggplot is used to visulaize the result.

# UMAP 
```
message("Running UMAP ...")
um_res <- umap(pc_scores, n_neighbors = 30, min_dist = 0.3)

umap_df <- data.frame(
  UMAP1 = um_res$layout[,1],
  UMAP2 = um_res$layout[,2],
  cell = rownames(pc_scores)
)

message("Running k-means clustering (k=3) ...")
set.seed(123)

km3 <- kmeans(umap_df[, c("UMAP1","UMAP2")], centers = 3, nstart = 50)
umap_df$cluster3 <- factor(km3$cluster)

p <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cluster3)) +
      geom_point(size = 1.5, alpha = 0.9) +
      theme_classic() +
      ggtitle("Drop-BS Single Cell Clustering (K-means = 3)") +
      scale_color_brewer(palette = "Set1")

print(p)

ggsave("dropbs_kmeans3_umap.png", p, width = 7, height = 6, dpi = 300)

message("Output saved: dropbs_umap.png")
```

# R output
<img width="1118" height="722" alt="d21a18d89d7886e4afec971d8b867f30" src="https://github.com/user-attachments/assets/10eb513e-576c-43b0-986a-d3fd5564232c" />
This UMAP visualization shows three well-separated clusters of single cells based on bin-level DNA methylation features. Each point represents one cell, and colors indicate cluster assignments derived from graph-based clustering. The clear separation between clusters suggests that global methylation patterns are sufficient to distinguish distinct cell populations.

# Potential mistakes for R
* UMAP projects high-dimensional features into two dimensions for visualization. If parameters such as n_neighbors and min_dist are not carefully tuned, clusters may appear poorly separated or misleading.
* Both UMAP and k-means involve randomness. Failing to set a random seed can lead to slightly different embeddings and cluster assignments across runs, making results difficult to reproduce.
