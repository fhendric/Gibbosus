# *dmrt* analysis

Note: analysis conducted within folder `~/Gibbosus/DMRT`

## 1. Relationship among all *O. gibbosus* *dmrt* genes

All putative *dmrt* genes of *O. gibbosus* were identified from the BRAKER-predicted gene set based on the presence of the term “Doublesex” in the functional annotation. Coding sequences were reconstructed using Iso-Seq transcripts when available; otherwise, transcripts predicted by StringTie were used.  

### 1.1 Selection and translation of *dmrt* isoseq transcripts

For *dmrt* genes where Iso-Seq transcripts were available (scaffold_11 and scaffold_39), unique clustered isoform transcripts (`./isoseq/OV210_03.flnc.clustered.hq.fasta`) mapping to the *dmrt* genes were manually selected in JBrowse. Four different isoforms were identified (iso1 - iso4) at the *dmrt* cluster at scaffold_11 and one at the *dmrt_G* cluster at scaffold_39. Alignment of the isoforms in MEGA revealed that the *dmrt_G* isoform corresponds to isoform 2 at the dmrt (scaffold_11) cluster. 
Isoforms were manually curated in MEGA and stored in `./DMRT/transcripts_isoseq`:

- `dmrt_isoseq_full.fasta` (original isoseq full length *dmrt* transcripts, including UTR)
- `dmrt_isoseq_coding.fasta` (coding sequence of *dmrt* transcripts)
- `dmrt_isoseq_coding.aa` (translated coding sequence of *dmrt* transcripts)

Visual inspection in MEGA showed that the coding sequences of isoform2 at scaffold_11 and scaffolds_39 aligned perfectly.  

### 1.2 Selection and translation of *dmrt* stringtie transcripts

Sequences of the *dmrt* genes that were not supported by Iso-Seq reads were retrieved from the StringTie predictions and manually selected in JBrowse. Transcript sequences were then translated using TransDecoder and stored in `./DMRT/transcripts_stringtie`::

```bash
module load TransDecoder
TransDecoder.LongOrfs -t dmrt_stringtie_full.fasta
TransDecoder.Predict  -t dmrt_stringtie_full.fasta --no_refine_starts
```
Resulting files `dmrt_stringtie_full.fasta.transdecoder.cds` and `dmrt_stringtie_full.fasta.transdecoder.pep` were then renamed, resulting in the following three StringTie transcript files:

- `dmrt_stringtie_full.fasta` (original stringtie full length *dmrt* transcripts, including UTR)
- `dmrt_stringtie_coding.fasta` (renamed `dmrt_stringtie_full.fasta.transdecoder.cds`file with coding sequence of *dmrt* transcripts)
- `dmrt_stringtie_coding.aa` (renamed `dmrt_stringtie_full.fasta.transdecoder.pep` file with translated coding sequence of *dmrt* transcripts)


## 2. Phylogenetic relationship of *dmrt* paralogs in relationship to the outgroup species

### 2.1 Selection and translation of of *dmrt* isoseq transcripts

=> See point 1.1.

### 2.2 Mapping of *dmrt* transcripts

Protein sequences of the *dmrt* isoseq transcripts were mapped to the genome with **miniprot**. Mapping were performed to both the entire genome - a reduced version of the genome was made that only includes the *dmrt* containg scaffolds i.e. scaffold_11 and scaffold_39 - and the G-locus (scaffold_39) only to identify the paralogous regions of the *dmrt* isoforms on the G-locus:

```bash
miniprot --gff ../genome/Ogib_2.0.reduced.fasta ../transcripts_isoseq/dmrt_isoseq_coding.aa > dmrt_isoseq_coding.gff
miniprot --gff ../genome/scaffold_39.fasta ../transcripts_isoseq/dmrt_isoseq_coding.aa > dmrt_isoseq_coding_scaf39.gff
```

### 2.3. Phylogenetic analysis of *dmrt* transcripts

We assessed the phylogenetic relationship between the two *dmrt* paralogs in relation to their sequence at the outgroup species. Due to the clear alignment of the different exons of isoform2, the analysis was restricted to this isoform only.  

#### 2.3.1. Preparation of BED files

Three BED files were generated from the `dmrt_isoseq_coding.gff` file for downstream analyses:

- `dmrt_iso2_coding.bed` (BED file with the exons of isoform2 at both s11 and s39)
- `dmrt_iso2_coding_s11.bed` (BED file with the exons of isoform2 at s11)
- `dmrt_iso2_coding_s39.bed` (BED file with the exons of isoform2 at s39)

#### 2.3.2. SNP calling

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
Output files are in./DMRT/vcf:

`Ogib2_0.dmrt_iso2.raw.vcf.gz`
`Ogib2_0.dmrt_iso2.noindel.vcf.gz`

#### 2.3.3. Reconstructing individual *dmrt* Sequences

Individual *dmrt* sequences were reconstructed based on the genotypes in `Ogib2_0.dmrt_iso2.noindel.vcf.gz` with `bcftools consensus` with the following script: 

```bash
#!/bin/bash
#PBS -N consensus
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=8

cd /kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/DMRT

# Load modules
module load BCFtools
module load SAMtools
module load BEDTools

# Read sample name from file based on array task ID
SAMPLE=$(sed -n "${PBS_ARRAYID}p" ./samples/samples.txt)

# Define directories and files
GENOME="/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/DMRT_snps/genome/Ogib_2.0.reduced.fasta"
VCF_NOINDEL="./vcf/Ogib2_0.dmrt_iso2.noindel.vcf.gz"
BED="./bed/dmrt_iso2_coding.bed"
CONSENSUS="./consensus"

# Generate consensus sequence
bcftools consensus --fasta-ref "$GENOME" --missing N --samples "$SAMPLE"  -o "$CONSENSUS/${SAMPLE}.fa" "$VCF_NOINDEL"

# Index consensus FASTA
samtools faidx "$CONSENSUS/${SAMPLE}.fa"

# Extract regions from fasta
bedtools getfasta -fi "$CONSENSUS/${SAMPLE}.fa" -bed ./bed/dmrt_iso2_coding_s11.bed -fo "$CONSENSUS/${SAMPLE}.scaffold_11.fa"
bedtools getfasta -fi "$CONSENSUS/${SAMPLE}.fa" -bed ./bed/dmrt_iso2_coding_s39.bed -fo "$CONSENSUS/${SAMPLE}.scaffold_39.fa"
```
Fasta files of scaffold_11 and scaffold_39 were then placed in the folders `./DMRT/consensus/scaffold_11` and `./DMRT/consensus/scaffold_39` respectively. A multifasta file that concatenates all sequences was generated with the script:

```bash
#!/bin/bash

# Folder containing your fasta files
input_folder="."   # change to your folder
output_file="dmrt_iso2_scaf39.fasta"

# Empty the output file if it exists
> "$output_file"

# Loop over all fasta files in the folder
for file in "$input_folder"/*.fasta "$input_folder"/*.fa; do
    # Get the filename without path and extension
    filename=$(basename "$file")
    individual_name="${filename%.*}"

    # Concatenate all sequences in the file into a single line, skipping any lines starting with ">"
    sequence=$(grep -v "^>" "$file" | tr -d "\n")

    # Write to the output file with the individual's name as header
    echo ">$individual_name" >> "$output_file"
    echo "$sequence" >> "$output_file"
done

echo "Concatenated multi-fasta written to $output_file"
```
#### 2.3.4. Obtaining the *dmrt* Sequences for outgroup species *Hylyphantes graminicola*

Analysis performed in the `./DMRT/hgram` folder. Map the *dmrt* coding region to the *H. graminicola* reference genome with miniprot:

```bash
cd ./DMRT/hgram
miniprot --gff ../../Hgram_genome/IOZCAS_Hgram_genomeAssembly_1.0.fa ../transcripts_isoseq/dmrt_isoseq_coding.aa > dmrt_iso2_hgram.gff
````
Generate a BED file from the `dmrt_iso2_hgram.gff` (`dmrt_iso2_hgram.bed`), which is then used to extract the sequence from the *H.graminicola* genome.

```bash
cd ./DMRT/hgram
bedtools getfasta -fi ../../Hgram_genome/IOZCAS_Hgram_genomeAssembly_1.0.fa -bed dmrt_iso2_hgram.bed -fo dmrt_iso2_hgram.fasta
```

#### 2.3.5. Generate aligned multifasta

DMRT_iso2 sequences of *Oedothorax* and *Hylyphantes* were combined and aligned with **MUSCLE** and stored in `./DMRT/multifasta/dmrt_iso2_multifasta.fasta`

#### 2.3.6. Phylogenetic analysis

Phylogenetic trees of the DMRT sequences were performed with IQtree and raxml.



```bash
module load RAxML-NG/1.2.0-GCC-12.3.0
raxml-ng --msa ../multifasta/dmrt_iso2_multifasta.fasta --model GTR+G --outgroup Hgram --bs-trees 1000 --threads 4 --prefix dmrt_iso2 --all
```

