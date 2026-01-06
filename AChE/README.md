# AChE analysis

## 1. Analysis of all AChE transcripts

### 1.1 Selection of AChE transcripts

In folder `./AChE/all`. We first selected all isoseq transcripts identified as "Cholinesterase" and "Acetylcholinesterase" using the command:
`grep 'holinesterase' ./blast_isoseq/Ogib_2.0_IsoSeq.swissprot.blastout.txt`
and generated a list of all respective isoseq genes, transcripts and isoforms. For transcripts, only the most abundant transcript for each gene was selected. 

`ache.isoseq.all.genes.list`       => List of all AChE genes (isoseq)

`ache.isoseq.all.isoforms.list`    => List of all AChE isoforms (isoseq)

`ache.isoseq.all.transcripts.list` => List of most abundant transcript of each AChE gene (isoseq)

Nucleotide sequences of all AChE transcripts were selected:
```bash
grep -w -f ache.isoseq.all.transcripts.list ../../isoseq/OV210_03.flnc.collapsed.hq.fasta > ache.isoseq.all.transcripts.fasta
```

### 1.2 Translation of AChE transcripts

Identification of the candidate coding sequences and the predicted protein sequences was performed using TransDecoder. 

```bash
module load TransDecoder/5.7.1-GCC-12.3.0
curl -L https://cpanmin.us | perl - App::cpanminus
cpanm install DB_File
cpanm install URI::Escape
TransDecoder.LongOrfs -t ache.isoseq.all.transcripts.fasta
TransDecoder.Predict -t ache.isoseq.all.transcripts.fasta
```
Output files are `ache.isoseq.all.transcripts.fasta.transdecoder.cds` and `ache.isoseq.all.transcripts.fasta.transdecoder.pep`

### 1.3 COBALT alignment of AChE transcripts

Peptide sequences were aligned using COBALT (https://www.ncbi.nlm.nih.gov/tools/cobalt/re_cobalt.cgi) and stored under `./all/cobalt/ache.isoseq.all.transcripts.fasta.transdecoder.cobalt.fasta`. Fasta headers were simplified using the command:

```bash
cat ache.isoseq.all.transcripts.fasta.transdecoder.cobalt.fasta | sed 's/^>[^|]*|/>/' | sed 's/GENE.*$//' > ache.isoseq.all.transcripts.fasta.transdecoder.cobalt.fasta
```
One sequence was clearly partial and removed: PB.7749.4:784184-799177(-)|transcript/17778.p1

Poorly aligned parts of the alignment were removed with trimAl:

```bash
module load trimAl/1.5.0-GCCcore-13.3.0
trimal -in ache.isoseq.all.transcripts.fasta.transdecoder.cobalt.fasta -out ache.isoseq.all.transcripts.fasta.transdecoder.cobalt.trimal.fasta -strict
```









