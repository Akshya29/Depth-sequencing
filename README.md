# Depth-sequencing
This script is used to do Peak Saturation processing and analysis. Peak Saturation consists of iteratively sub-sampling ChIP-seq bam files and repeatedly calling peaks to detect at what read depth the peak call hits a point of diminishing returns, i.e. it has reached saturation.
It requires macs2 and samtools.
```
bash run_subsample_macs2_pipeline.sh \
      -i /path/to/bam_dir \
      -o /path/to/output \
      -g 2.1e9 \
      -m samtools/1.18-gcc-12.3.0
```
This creates a peak_counts.tsv file which can be used to plot the saturation curve:
```
df <- read.table("peak_counts.tsv", header=TRUE, sep="\t")

ggplot(df, aes(x=depth, y=peaks, color=genotype)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x="reads (%)", y="Number of peaks")
```
