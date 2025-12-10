# Annotation

The Ogib_2.0 genome was annotated using the following gene prediction tools:

- [BRAKER/AUGUSTUS](./Braker/): *abinitio* gene prediction using BRAKER/AUGUSTUS trained with mRNAseq data
- [Protein mappings](./Proteins/): mapping protein sequences from the related species *H. graminicola* with **miniprot**
- [StringTie](./StringTie/): transcript assembly of short-read mRNAseq using **Stringtie**
- [IsoSeq](./IsoSeq/): mapping of PacBio HiFi reads of full-length mRNA (IsoSeq)

The different types of evidence were then compined with **EVidenceModeler (EVM)** to obtain a final set of gene predictions. 
- [EvidenceModeler](./EVidenceModeller/): mapping of PacBio HiFi reads of full-length mRNA (IsoSeq)

