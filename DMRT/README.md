# Workflow to reconstruct DMRT sequences

### 1. Selection of dmrt transcripts

Unique clustered isoform transcripts (`./isoseq/OV210_03.flnc.clustered.hq.fasta`) mapping to the dmrt genes were manually selected in JBrowse. Four different isoforms were identified (iso1 - iso4) at the dmrt cluster at scaffold_11 and one at the dmrt_G cluster. Alignment of the isoforms in MEGA revealed that the dmrt_G isoform corresponded to isoform 2 at the dmrt (scaffold_11) cluster. 
Isoforms were manually curated in MEGA and stored in `./DMRT/transcripts_isoseq`:

- `dmrt_isoseq_full.fasta` (original isoseq full length dmrt transcripts, including UTR)
- `dmrt_isoseq_coding.fasta` (coding sequence of dmrt transcripts)
- `dmrt_isoseq_coding.aa` (translated coding sequence of dmrt transcripts)

### 2. Mapping of dmrt transcripts


```bash
#!/bin/bash
#PBS -N minimap2_isoseq
#PBS -l walltime=08:00:00
#PBS -l nodes=1:ppn=8

module load minimap2
module load SAMtools

cd /kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/isoseq

GENOME=/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/fasta/Ogib_2.0.fasta
READS=/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/isoseq/OV210_03.flnc.fasta
BAM_OUT=/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/minimap/IsoSeq_vs_Ogib_2.0/OV210_03.flnc.minimap2.bam


minimap2 -ax splice:hq --secondary=no -uf "$GENOME" "$READS" | samtools view -b | samtools sort -o "$BAM_OUT"
samtools index "$BAM_OUT"
````
Mappings were visualized in JBROWSE and transcript reads mapping to DMRT at scaffold_11 (dmrt) and scaffold_39 (dmrt_G) representing a full length unique haplotype or isoform were selected manually and stored in `./DMRT_snps/transcripts_flnc/transcripts_flnc_dmrt.fasta`.

Transcript reads were aligned (MEGA) and conserved transcript regions showing unambiguous alignment were selected and stored in `./DMRT_snps/transcripts_flnc/transcripts_flnc_dmrt_conserved.fasta`  



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





