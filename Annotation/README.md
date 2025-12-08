# Annotation

Genome annotation using (i) *abinitio* gene prediction using BRAKER/AUGUSTUS, (ii) proteins from *H. graminicola*, (iii) transcript assembly using Stringtie and (iv) mapping of PacBio HiFi mRNA reads. The different types of evidence were then compined with EVidenceModeler to obtain a final set of gene predictions. 

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

## EVidenceModeler (EVM)

Before running EVM, GTF/GFF files were first parsed to GFF3 format that is compatible with EVM. 

### Parse Braker GTF file

Parsing the Braker GTF file `augustus.hints.mRNA.gtf` was done
### Parse Protein GTF file
