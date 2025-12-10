# EVidenceModeler (EVM)

Before running EVM, GTF/GFF files were first parsed to EVM compatible GFF3 format. The following directory structure was created to store the EVM-compatible GFF3 files:

```
evm_inputs/
├── abinitio/
├── transcripts/
└── proteins/
```

To check if a gff3 file is compatible for EVM, you can use the `/gff3_gene_prediction_file_validator.pl` tool available in the EvmUtils in EVM: 

```bash
module load EVidenceModeler/2.1.0-foss-2024a
$EVM_HOME/EvmUtils/gff3_gene_prediction_file_validator.pl input.gff3
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

$EVM_HOME/EvmUtils/misc/taco_gtf_to_alignment_gff3.pl OV210_03.flnc.collapsed.hq.gff > test

