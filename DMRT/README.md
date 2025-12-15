# Workflow to reconstruct DMRT sequences

### 1. Select flnc IsoSeq transcripts

Full-length non-chimeric IsoSeq reads were mapped to the genome with minimap2.


### 2. Extract reads mapping to the 

```bash
#!/bin/bash

#PBS -N samtools
#PBS -l walltime=04:00:00
#PBS -l nodes=1:ppn=8

module load SAMtools/1.18-GCC-12.3.0

cd /kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/

# Read sample name from file based on array task ID
SAMPLE=$(sed -n "${PBS_ARRAYID}p" ./samples/samples.txt)

# Make directories to store reads
#mkdir ./DMRT_denovo/reads/"$SAMPLE"

# Define files
BED=/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/DMRT_denovo/bed/dmrt_conserved.bed
BAM_IN="/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_reseq/raw"
BAM_OUT="/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/DMRT_denovo/bam"
FASTQ_OUT="/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/DMRT_denovo/fastq"

# Extract reads
samtools view -b -L "$BED" "$BAM_IN/${SAMPLE}.bam" | samtools sort -n > "$BAM_OUT/${SAMPLE}.dmrt.bam"
samtools index "$BAM_OUT/${SAMPLE}.dmrt.bam"
samtools view -b "$BAM_OUT/${SAMPLE}.dmrt.bam" | samtools fastq -1 "$FASTQ_OUT/${SAMPLE}_R1.dmrt.fastq" -2 "$FASTQ_OUT/${SAMPLE}_R2.dmrt.fastq" -0 /dev/null -s /dev/null -n
```





