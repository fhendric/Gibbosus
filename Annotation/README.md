# Annotation

Genome annotation using:
(i) *abinitio* gene prediction using BRAKER/AUGUSTUS
(ii) mapping proteins from the related species *H. graminicola* with **miniprot**
(iii) transcript assembly using **Stringtie** 
(iv) mapping of PacBio HiFi mRNA reads. 

The different types of evidence were then compined with **EVidenceModeler (EVM)** to obtain a final set of gene predictions. 

## Braker/Augustus

Braker was run using the following script:
```bash
#!/bin/bash
 
 
#PBS -N braker_mRNA
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=9
 
 
module load Perl/5.32.0-GCCcore-10.2.0
module load BRAKER/2.1.6-foss-2020b
module load AUGUSTUS/3.4.0-foss-2020b
module load DIAMOND/2.0.7-GCC-10.2.0
module load CDBtools/0.99-GCC-10.2.0
 
 
cd /kyukon/scratch/gent/vo/000/gvo00032/PacBio_gibo_wtdbg2/braker_mRNA
 
export AUGUSTUS_CONFIG_PATH=/kyukon/scratch/gent/vo/000/gvo00032/PacBio_gibo_wtdbg2/AUGUSTUS-config/3.3.2/config
 
braker.pl --cores=9 --species=ogibo_mRNA --genome=/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/fasta/Ogib_2.0.red_masked.fa --bam=/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV200_03_aft.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV200_17_aft.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV204_25_aft.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV207_01_afg.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV207_02_aft.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV207_09_amt.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV207_12_afg.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV207_13_afg.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV207_15_amg.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV207_17_aft.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV207_19_amg.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV207_22_amg.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV207_26_amt.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV207_28_amg.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV209_01_amt.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV209_07_afg.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV209_13_amt.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV209_16_amg.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV209_24_amt.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV210_03_amg.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV210_05_amt.bam,/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/bam_mRNA/OV210_07_afg.bam --softmasking
```

## Protein mappings

## Stringtie transcripts (short-read mRNA)

## IsoSeq transcripts (PacBio HiFi mRNA reads)

### 1. Used files:

Raw ccs read files (Macrogen):

```
OV210_03.hifi_reads.bam
OV210_03.hifi_reads.bam.pbi
OV210_03_HiFi.fastq.gz
```

### 2. Primer removal

Reads were demultiplexed by Macrogen (multiplex barcode removed), but still contain the cDNA primers and polyA tail. These were removed with **lima** (performed at HPC KULeuven as **lima** could not be installed on the HPC).

Used primer sequences are:

```
>primer_5p
GCAATGAAGTCGCAGGGTTGGG
>primer_3p
AAGCAGTGGTATCAACGCAGAGTAC
```

Resulting files are:

```
OV210_03.hifi_reads.demux.bam
OV210_03.hifi_reads.demux.bam.pbi
```


### 3. Clean and polyA-tail removal

Reads were then cleaned and polyA tail removed with isoseq refine (`run_isoseq.sh`):

```bash
module load Isoseq
isoseq refine --require-polya OV210_03.hifi_reads.demux.bam primers.fasta OV210_03.flnc.bam
```

Resulting files are:

```
OV210_03.flnc.bam
OV210_03.flnc.bam.pbi
OV210_03.flnc.consensusreadset.xml
OV210_03.flnc.filter_summary.report.json
OV210_03.flnc.report.csv
````

### 4. Clustering transcripts

Transcripts were subsequently clustered with `isoseq cluster`:

```
isoseq cluster OV210_03.flnc.bam OV210_03.flnc.clustered.bam
```

with resulting files:

```
OV210_03.flnc.clustered.bam
OV210_03.flnc.clustered.bam.pbi
```

and all remaining `OV210_03.flnc.clustered* files`.


### 5. Mapping of clustered transcripts

Clustered transcripts were mapped to the genome with minimap2.

```
minimap2 -t 9 -ax splice:hq -uf /kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/fasta/Ogib_2.0.fasta /kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/isoseq/OV210_03.flnc.clustered.hq.fasta.gz | samtools sort -@10 -o Ogib_2.0.isoseq.clustered.hq.bam
```
with resulting files:

```
Ogib_2.0.isoseq.clustered.hq.bam
Ogib_2.0.isoseq.clustered.hq.bam.bai
```


### 6. Mapping reads with pbmm2

Install pbmm2

```bash
module load Miniconda3
conda create --prefix /user/gent/404/vsc40419/scratch_vo/Gibbosus/isoseq/pbmm2  pbmm2
```

Index genome

```bash
../isoseq/pbmm2/bin/pbmm2 index Ogib_2.0.fasta Ogib_2.0.mmi
```

Map flnc reads to genome

```bash
./pbmm2/bin/pbmm2 align ../fasta/Ogib_2.0.mmi OV210_03.flnc.bam Ogib_2.0.flnc.pbmm2.bam --preset HIFI --sort
```

Map clustered reads to genome

```bash
./pbmm2/bin/pbmm2 align ../fasta/Ogib_2.0.mmi OV210_03.flnc.clustered.hq.bam Ogib_2.0.flnc.clustered.hq.pbmm2.bam --preset ISOSEQ --sort
```

### 7. Collapse reads and generate GFF

```bash
isoseq collapse --do-not-collapse-extra-5exons Ogib_2.0.flnc.clustered.hq.pbmm2.bam OV210_03.flnc.bam OV210_03.flnc.collapsed.hq.gff
```

### 8. Add read counts
Make table to convert collapsed versus clustered names and add read counts

```bash
grep '>' OV210_03.flnc.collapsed.hq.fasta | sed 's/|P.*|/ /g' | sed 's/>//g' > OV210_03.flnc.clustered2collapsed_conversion.hq.txt
less OV210_03.flnc.collapsed.hq.flnc_count.txt | sed 's/^.*,//g' > abundances.txt
paste  OV210_03.flnc.clustered2collapsed_conversion.hq.txt abundances.txt > brol.txt
mv brol.txt OV210_03.flnc.clustered2collapsed_conversion.hq.txt
```

### 9. Final GFF

```bash
OV210_03.flnc.collapsed.hq.gff
```

## EVidenceModeler (EVM)

Before running EVM, GTF/GFF files were first parsed to EVM compatible GFF3 format. The following directory structure was created to store the EVM-compatible GFF3 files:

```
evm_inputs/
├── abinitio/
├── transcripts/
└── proteins/
```

### Parse Braker GTF file

GTF file `augustus.hints.mRNA.gtf` produced by **BRAKER** was parsed using the `braker_GTF_to_EVM_GFF3.pl` tool available in the EvmUtils in EVM using the following script:
```bash
module load EVidenceModeler/2.1.0-foss-2024a
$EVM_HOME/EvmUtils/misc/braker_GTF_to_EVM_GFF3.pl augustus.hints.mRNA.gtf > ./evm_inputs/abinitio/Ogibo.braker.evm.gff3
```
### Parse miniprot Protein GFF file

GFF2 file `Ogib_2.0_HgramProt.gff` produced by **miniprot** was parsed using the `miniprot_GFF_2_EVM_GFF3.py` tool available in the EvmUtils in EVM using the following script:
```bash
module load EVidenceModeler/2.1.0-foss-2024a
python $EVM_HOME/EvmUtils/misc/miniprot_GFF_2_EVM_GFF3.py Ogib_2.0_HgramProt.gff > ./evm_inputs/proteins/Ogibo.HgramProt.evm.gff3
```

### Parse StringTie transcripts file

The merged **StringTie** output `stringtie_merged.gtf`contains both the transcripts assembled using mRNAseq data and the non-overlapping *abinitio* gene predictions from Braker/Augustus. First remove the Braker/Augustus records as they are already included as evidence the Braker GTF file:

`grep -v 'AUGUSTUS' stringtie_merged.gtf > stringtie.gtf`




