# Conjugation mediates large-scale chromosomal transfer in Streptomyces driving diversification of antibiotic biosynthetic gene clusters
Caroline Choufa<sup>1</sup>, Pauline Gascht<sup>1</sup>, Hugo Leblond<sup>1</sup>, Anthony Gauthier<sup>1</sup>, Michiel Vos<sup>2</sup>, Cyril Bontemps<sup>1*</sup>, Pierre Leblond<sup>1*</sup>

<sup>1</sup>Université de Lorraine, INRAe, DynAMic, F-54000, Nancy, France.

<sup>2</sup>European Centre for Environment and Human Health, University of Exeter Medical School, Environment and Sustainability Institute, Penryn Campus, TR10 9FE, United Kingdom.

<sup>*</sup>co-corresponding authors


## Abstract
Streptomyces are ubiquitous soil dwelling bacteria of special importance as a source of metabolites used in human and veterinary medicine, agronomy and industry. Conjugation is the main mechanism of Streptomyces Horizontal Gene Transfer, and this process has long been known to be accompanied by mobilization of chromosomal DNA. However, the magnitude of DNA transfer, or the localization of acquired DNA across the chromosome, has remained undetermined. We here show that conjugative crossings in sympatric strains of Streptomyces result in the large-scale, genome-wide distributed replacement of up to one third of the recipient chromosome, a phenomenon for which we propose the name ‘Streptomyces Chromosomal Transfer’ (SCT). Such chromosome blending resulted in the acquisition, loss and hybridization of antibiotic Biosynthetic Gene Clusters leading to a novel metabolic arsenal in conjugant offspring. Harnessing conjugation-mediated BGC diversification holds great promise in the search for new Streptomyces functional diversity.

## Introduction
To decipher how- and to what extent the transfer of conjugative elements promotes the transfer of chromosomal DNA in Streptomyces, we performed mating experiments with closely related strains isolated, selected recombinants and sequenced them at high coverage (MiSeq sequencing, Illumina, CA, USA). The sequence for each recombinant strain was deduced from the alignment of sequencing reads on the parental chromosomes and SNP analysis was used to identify the parental origin of each recombinant region. The genome data are available on the NCBI database for the parental strains and raw sequencing data for the recombinant strains were deposited at SRA (NCBI, Bioproject PRJNA912173). The alignments were performed using the Geneious prime platform reference algorithm (Invitrogen Corp., version 2022.0.1). The consensus of the recombinant sequence was aligned with the chromosomal sequences of the parents using the progressive Mauve algorithm. The Single Nucleotide Polymorphism (SNP) call table was exported with MAUVE tool ‘export SNPs’. Custom scripts (Python V.3.9) were used to call SNPs of both parental strains along the alignment (Github). The presence of two or more successive SNPs was used to assign a genomic region to one of the two parents (to buffer against possible effects of (single) point mutations).


## Cloning the Repository
The repository contains submodules, thus please check it out with 
```shell
# SSH
git clone git@github.com:leblond14u/Large-scale-chromosomal-DNA-transfer.git
```
or
```shell
# HTTPS
git clone https://github.com/leblond14u/Large-scale-chromosomal-DNA-transfer.git
```

## Overview
The codebase has **2 main** components:

A numpy based python backend file **SNP.py** containing the following functions :

- **SNP_parent()** : Take a SNP list from MAUVE and label each SNP as P1 or P2 (parent 1 or parent 2 of the crossing).
- **SNP_chaine()** : Define some SNP blocs depending on the number of successive SNP to switch from P1 to P2 as well as from P2 to P1.
- **compteur_snp()** : Count the number of SNP in each blocs of SNP.
- **chaine_borne()** : Convert each labelled SNP blocs to regions of coordinates (i.e. P1 regions, start, end).
- **calcul_des_sommes_rapport_recombinant()** : Compute the sum of the regions size relative to the receiver.
- **calcul_des_sommes_rapport_donneur()** : Compute the sum of the regions size relative to the donor.
- **creer_excel()** : Creates a Excel file containing the computed informations.


A frontend python file **Main.py** :

- The frontend file is used as a user interface. Two MAUVE files sorted by P1 and P2 are referenced. The single nucleotide polymorphism (SNP) call table exported from MAUVE's "SNP export" tool is sorted according to donor and recipient genome positions to give the excel files 'name_triéP1' and 'name_triéP2' respectively. These files are then used as input files. The script then assigns the donor parent and the recipient parent according to the sum of the SNP region sizes. And then uses the **SNP.py** functions to create an output excel file with the right informations.


