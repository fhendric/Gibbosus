# *dmrt* analysis

Note: analysis conducted within folder `~/Gibbosus/DMRT`

## 1. Identfication of all putative *O. gibbosus* *dmrt* genes

All putative *dmrt* genes of *O. gibbosus* were identified from the BRAKER-predicted gene set (i) based on the presence of the term “Doublesex” in the functional annotation and (ii) by aligning the Dmrt sequences of *Parasteatoda tepidariorum*. If available, coding sequences were reconstructed using Iso-Seq transcripts, StringTie transcripts and aligned *P. tepidariorum* sequences.  

### 1.1 Selection and translation of *dmrt* isoseq transcripts

For *dmrt* genes where Iso-Seq transcripts were available (scaffold_11 and scaffold_39), unique clustered isoform transcripts (`./isoseq/OV210_03.flnc.clustered.hq.fasta`) mapping to the *dmrt* genes were manually selected in JBrowse. Multiple isoforms were identified at the *dmrt* cluster at scaffold_11 and the hunch-specific *dmrt* cluster at scaffold_39. 
Selected IsoSeq transcripts are stored in `./DMRT/transcripts/transcripts_isoseq`:

- `dmrt_isoseq_full.fasta` (original isoseq full length *dmrt* transcripts, including UTR)
- `dmrt_isoseq_coding.fasta` (coding sequence of *dmrt* transcripts)
- `dmrt_isoseq_coding.aa` (translated coding sequence of *dmrt* transcripts)

### 1.2 Selection and translation of *dmrt* stringtie transcripts

Sequences of the *dmrt* genes that were not supported by Iso-Seq reads were retrieved from the StringTie predictions and manually selected in JBrowse. Transcript sequences were then translated using TransDecoder and stored in `./DMRT/transcripts/transcripts_stringtie`:

```bash
module load TransDecoder
TransDecoder.LongOrfs -t dmrt_stringtie_full.fasta
TransDecoder.Predict  -t dmrt_stringtie_full.fasta --no_refine_starts
```
Resulting files `dmrt_stringtie_full.fasta.transdecoder.cds` and `dmrt_stringtie_full.fasta.transdecoder.pep` were then renamed, resulting in the following three StringTie transcript files:

- `dmrt_stringtie_full.fasta` (original stringtie full length *dmrt* transcripts, including UTR)
- `dmrt_stringtie_coding.fasta` (renamed `dmrt_stringtie_full.fasta.transdecoder.cds`file with coding sequence of *dmrt* transcripts)
- `dmrt_stringtie_coding.aa` (renamed `dmrt_stringtie_full.fasta.transdecoder.pep` file with translated coding sequence of *dmrt* transcripts)

### 1.3 Selection of *dmrt* genes without transcript evidence

Sequences of the *dmrt* genes that were not supported by Iso-Seq of StringTie predictions were obtained from the Braker predicted coding sequences
  
### 1.4 Compilation of all *dmrt* genes

Sequences of all *dmrt* genes were combined and stored in ´./DMRT/transcripts/transcripts_all/` and renamed based on homology or alignment with *P. tepidariorum* sequences resulting in the following 14 Dmrt genes:

- Dmrt11E
- Dmrt93B
- Dmrt99B
- DmrtL1A-B
- DmrtL2A-G
- DmrtL2B_G and DmrtL2D_G

The final set af all *dmrt* genes are compiled in the following files:
- `dmrt_all_full.fasta` (full length sequences of all *dmrt* genes, including UTR)
- `dmrt_all_coding.fasta` (coding sequences of all *dmrt* genes)
- `dmrt_all_coding.aa` (peptide sequences of all *dmrt* genes)

We also aligned the entire set of *dmrt* genes to the genome to identify the exact location of the coding sequences on the Ogib.2.0 genome:

```bash
miniprot ../fasta/Ogib_2.0.fasta ./transcripts/transcripts_all/dmrt_all_coding.aa --gtf --trans > ./miniprot/dmrt_vs_genome/dmrt_all_coding.gtf
```

### 1.5 Alignment of *O.gibbosus* and *P. tepidariorum* *dmrt* transcripts

Peptide sequences of the *O. gibbosus* and *P. tepidariorum* *dmrt* genes were combined and aligned with COBALT. The alignment was stored under ´./DMRT/transcripts/transcripts_all/dmrt_all_coding.cobalt.fasta`

### 1.6 Identification of putative *dmrtL*2 paralogs in the hunch-specific sequence

We identified the presence of all DmrtL2 genes within the hunch-specific sequence by aligning the *dmrtL*2 coding sequences to the hunch-specific scaffold harboring the duplicated *dmrtL*2 cluster (scaffold_39) with miniprot:

```bash
miniprot scaffold_39.fasta ./transcripts/transcripts_all/dmrt_all_coding.aa --gtf --trans > ./miniprot/dmrt_vs_scaffold_39/dmrt_all_coding.s_39.gtf
```



## 2. Phylogenetic relationship of *dmrt* paralogs in relationship to the outgroup species

### 2.1. Phylogenetic analysis of *dmrt* transcripts

We assessed the phylogenetic relationship between the two *dmrt* paralogs in relation to their sequence at the outgroup species. 

#### 2.3.1. Preparation of BED files

A BED files was generated that includes all exons of the *dmrtL*2 clusters on scaffold_11 and scaffold_39. Location of the exons was inferred from the miniprot mappings specified in 1.4 (´./miniprot/dmrt_vs_genome/dmrt_all_coding.gtf´) and 1.6 (`./miniprot/dmrt_vs_scaffold_39/dmrt_all_coding.s_39.gtf`).

```bash
./bed/dmrtL2_CDS.bed
```
#### 2.3.2. SNP calling

Reconstruction of the *dmrt* exons for all individuals, including outgroups, was performed via SNP calling using BCFtools. Because the *dmrt* gene has two paralogs, sequencing reads from outgroup species, presumed to have only a single *dmrt* copy, may map equally well to both copies, resulting in low mapping quality scores. To prevent these reads from being discarded, we used raw BAM files without filtering for low-quality mappings. SNPs were called with BCFtools and restricted to the exonic regions of interest. We further required that genotypes be called as heterozygous or homozygous for the alternative allele only if supported by at least two reads carrying the alternative allele; positions with only a single supporting read were called homozygous for the reference allele. This approach minimizes the risk of calling sequencing errors as heterozygous sites. SNP calling and filtering were performed using the following script:

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
VCF_RAW="./vcf/Ogib2_0.dmrtL2.raw.vcf.gz"
VCF_NOINDEL="./vcf/Ogib2_0.dmrtL2.noindel.vcf.gz"
BED="./bed/dmrtL2_CDS.bed"

# SNP calling
bcftools mpileup -R "$BED" -a AD,DP,SP -Ou -f "$GENOME" "$BAM"/D1086_G.bam "$BAM"/D1090_G.bam "$BAM"/D1091_T.bam "$BAM"/D601_T.bam "$BAM"/H002_T.bam "$BAM"/H007_G.bam "$BAM"/Ofusc.bam "$BAM"/Oretu.bam "$BAM"/Otril.bam "$BAM"/OV001_G.bam "$BAM"/OV002_T.bam "$BAM"/OV208_T.bam "$BAM"/OV213_G.bam "$BAM"/PO002_T.bam "$BAM"/PO004_G.bam "$BAM"/PU001_G.bam "$BAM"/PU002_T.bam "$BAM"/SE001_G.bam "$BAM"/SE003_T.bam "$BAM"/W791_G.bam "$BAM"/W815_T.bam "$BAM"/W816_G.bam "$BAM"/W818_T.bam | bcftools call -m -f GQ,GP -Oz -o "$VCF_RAW"

bcftools view -V indels,mnps "$VCF_RAW" -Ou | bcftools +setGT -Ou  -- -t q -n 0 -i 'FMT/AD[*:1]<2' | bcftools +setGT -Oz -o "$VCF_NOINDEL"  -- -t q -n . -i 'FMT/DP<1'
```
Output files are in./DMRT/vcf:

`Ogib2_0.dmrtL2.raw.vcf.gz`
`Ogib2_0.dmrtL2.noindel.vcf.gz`

#### 2.3.3. Reconstructing individual *dmrt* Sequences

Individual *dmrt* sequences were reconstructed based on the genotypes in `Ogib2_0.dmrtL2.noindel.vcf.gz` with `bcftools consensus` with the following script: 

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
GENOME="/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/DMRT/genome/Ogib_2.0.reduced.fasta"
VCF_NOINDEL="./vcf/Ogib2_0.dmrtL2.noindel.vcf.gz"
BED="./bed/dmrtL2_CDS.bed"
CONSENSUS="./consensus"

# Generate consensus sequence
bcftools consensus --fasta-ref "$GENOME" --missing N --samples "$SAMPLE"  -o "$CONSENSUS/${SAMPLE}.fa" "$VCF_NOINDEL"

# Index consensus FASTA
samtools faidx "$CONSENSUS/${SAMPLE}.fa"

# Extract regions from fasta
bedtools getfasta -fi "$CONSENSUS/${SAMPLE}.fa" -bed ./bed/dmrtL2_CDS_s11.bed -fo "$CONSENSUS/${SAMPLE}.scaffold_11.fa"
bedtools getfasta -fi "$CONSENSUS/${SAMPLE}.fa" -bed ./bed/dmrtL2_CDS_s39.bed -fo "$CONSENSUS/${SAMPLE}.scaffold_39.fa"
```
Fasta files of scaffold_11 and scaffold_39 were then placed in the folders `./DMRT/consensus/scaffold_11` and `./DMRT/consensus/scaffold_39` respectively. A multifasta file that concatenates all sequences was generated with the script:

```bash
#!/bin/bash

CONSENSUS="./consensus"

# Folder containing fasta files
input_folder="$CONSENSUS"  
output_file="dmrtL2_scaf39.fasta"

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

