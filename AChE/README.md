# AChE analysis

### Selection of AChE transcripts

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

### Alignment of AChE transcripts





