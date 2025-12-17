# *dmrt* analysis

### 1. Selection of *dmrt* transcripts

Unique clustered isoform transcripts (`./isoseq/OV210_03.flnc.clustered.hq.fasta`) mapping to the *dmrt* genes were manually selected in JBrowse. Four different isoforms were identified (iso1 - iso4) at the *dmrt* cluster at scaffold_11 and one at the *dmrt_G* cluster. Alignment of the isoforms in MEGA revealed that the *dmrt_G* isoform corresponded to isoform 2 at the dmrt (scaffold_11) cluster. 
Isoforms were manually curated in MEGA and stored in `./DMRT/transcripts_isoseq`:

- `dmrt_isoseq_full.fasta` (original isoseq full length *dmrt* transcripts, including UTR)
- `dmrt_isoseq_coding.fasta` (coding sequence of *dmrt* transcripts)
- `dmrt_isoseq_coding.aa` (translated coding sequence of *dmrt* transcripts)

Visual inspection in MEGA showed that the coding sequences of isoform2 at scaffold_11 and scaffolds_39 aligned perfectly.  

### 2. Mapping of *dmrt* transcripts

Protein sequences of the *dmrt* isoseq transcripts were mapped to the genome with **miniprot**. Mapping were performed to both the entire genome - a reduced version of the genome was made that only includes the *dmrt* containg scaffolds i.e. scaffold_11 and scaffold_39 - and the G-locus (scaffold_39) only to identify the paralogous regions of the *dmrt* isoforms on the G-locus:

```bash
miniprot --gff ../genome/Ogib_2.0.reduced.fasta ../transcripts_isoseq/dmrt_isoseq_coding.aa > dmrt_isoseq_coding.gff
miniprot --gff ../genome/scaffold_39.fasta ../transcripts_isoseq/dmrt_isoseq_coding.aa > dmrt_isoseq_coding_scaf39.gff
```

### 3. Phylogenetic analysis of *dmrt* transcripts

We assessed the phylogenetic relationship between the two *dmrt* paralogs in relation to their sequence at the outgroup species. Due to the clear alignment of the different exons of isoform2, the analysis was restricted to this isoform only.  

#### 3.1. Preparation of BED files

Three BED files were generated for downstream analyses were generated:

- `dmrt_iso2_coding.bed` (BED file with the exons of isoform2 at both s11 and s39)
- `dmrt_iso2_coding_s11.bed` (BED file with the exons of isoform2 at s11)
- `dmrt_iso2_coding_s39.bed` (BED file with the exons of isoform2 at s39)

#### 3.2. SNP calling

Reconstruction of the *dmrt* isoform 2 sequences for all individuals, including outgroups, was performed via SNP calling using BCFtools. Because the *dmrt* gene has two paralogs, sequencing reads from outgroup species, presumed to have only a single *dmrt* copy, may map equally well to both copies, resulting in low mapping quality scores. To prevent these reads from being discarded, we used raw BAM files without filtering for low-quality mappings. SNPs were called with BCFtools and restricted to the exonic regions of interest. We further required that genotypes be called as heterozygous or homozygous for the alternative allele only if supported by at least two reads carrying the alternative allele; positions with only a single supporting read were called homozygous for the reference allele. This approach minimizes the risk of calling sequencing errors as heterozygous sites. SNP calling and filtering were performed using the following script:

```bash
#!/bin/bash
#PBS -N bcftools
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=8

cd /kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/DMRT

# Load modules
module load BCFtools

# Define directories and files
GENOME="/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/fasta/Ogib_2.0.fasta"
BAM="/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_reseq/raw"
VCF_RAW="./vcf/Ogib2_0.dmrt_iso2.raw.vcf.gz"
VCF_NOINDEL="./vcf/Ogib2_0.dmrt_iso2.noindel.vcf.gz"
BED="./bed/dmrt_iso2_coding.bed"

# SNP calling
bcftools mpileup -R "$BED" -a AD,DP,SP -Ou -f "$GENOME" "$BAM"/D1086_G.bam "$BAM"/D1090_G.bam "$BAM"/D1091_T.bam "$BAM"/D601_T.bam "$BAM"/H002_T.bam "$BAM"/H007_G.bam "$BAM"/Ofusc.bam "$BAM"/Oretu.bam "$BAM"/Otril.bam "$BAM"/OV001_G.bam "$BAM"/OV002_T.bam "$BAM"/OV208_T.bam "$BAM"/OV213_G.bam "$BAM"/PO002_T.bam "$BAM"/PO004_G.bam "$BAM"/PU001_G.bam "$BAM"/PU002_T.bam "$BAM"/SE001_G.bam "$BAM"/SE003_T.bam "$BAM"/W791_G.bam "$BAM"/W815_T.bam "$BAM"/W816_G.bam "$BAM"/W818_T.bam | bcftools call -m -f GQ,GP -Oz -o "$VCF_RAW"

bcftools view -V indels,mnps "$VCF_RAW" -Ou | bcftools +setGT -Ou  -- -t q -n 0 -i 'FMT/AD[*:1]<2' | bcftools +setGT -Oz -o "$VCF_NOINDEL"  -- -t q -n . -i 'FMT/DP<1'
```






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





