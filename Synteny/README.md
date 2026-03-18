# Synteny analysis between Ogib_2.0, Ogib_1.0 and Hgram genomes

## 1. Map sequences with minimap

Genomic sequences of Ogib_2.0 genome were mapped against Ogib_1.0 and Hgram genomes using following script:

```bash
#!/bin/bash
#PBS -N minimap2
#PBS -l walltime=04:00:00
#PBS -l nodes=1:ppn=8

module load minimap2


GENOME=/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/fasta/Ogib_2.0.fasta
#QUERY=/kyukon/scratch/gent/vo/000/gvo00032/PacBio_gibo_wtdbg2/fasta/Ogib_1.0.fasta
QUERY=/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/Hgram_genome/IOZCAS_Hgram_genomeAssembly_1.0.fa
#OUT=/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/minimap/Ogib_1.0_vs_Ogib_2.0/Ogib_1.0_vs_Ogib_2.0.synteny.paf
OUT=/kyukon/scratch/gent/vo/000/gvo00032/Gibbosus/minimap/Ogib_2.0_vs_Hgram/Ogib_2.0_vs_Hgram.synteny.paf

minimap2 -x asm10 -N 0 $GENOME $QUERY > $OUT
```

## 2. Filter paf

Unplaced contigs and alignment blocks < 20kb were filter from the paf files using the following commands.

For the ```Ogib_1.0_vs_Ogib_2.0.synteny.paf``` :

```less Ogib_1.0_vs_Ogib_2.0.synteny.paf | grep 'CM0332' | sed 's/scaffold_//'g | awk '$6<14' | awk '$11>20000' > Ogib_1.0_vs_Ogib_2.0.synteny.filt.paf```
