# IsoSeq transcripts (PacBio HiFi mRNA reads)

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

Reads were then cleaned and polyA tail removed with isoseq refine (`run_isoseq.sh`) to generate **FLNC (Full-Length Non-Chimeric)** reads:

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
paste  OV210_03.flnc.clustered2collapsed_conversion.hq.txt abundances.txt > OV210_03.flnc.clustered2collapsed_conversion.hq.txt
```

### 9. Final GFF

```bash
OV210_03.flnc.collapsed.hq.gff
```

