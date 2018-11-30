# SNPore
Forensic Nanopore SNP sequencing

## Background

DNA analysis has become the cornerstone of contemporary forensic science. Altough most forensic DNA testing still use PCR and capillary electrophoresis (CE)-based analysis methods, a transition to sequencing techniques is at hand.
Next generation targetted sequencing allows forensic scientists worldwide to harness the full potential of their DNA sample. 

## Aim

The aim of the paper (Forensic SNP Analysis using Nanopore MinION™ Sequencing) associating this repository was to asses the capability of Oxford Nanopore Technologies’ (ONT) handheld sequencer MinION™ to created a forensic SNP based profile. 

## Motivation

The cloud based base calling software Metrichor being used by ONT offers a fully integrated, bespoke, real-time analysis solution. However, as ONT sequencing was primarerly designed to analyze long reads (>1kb) the base calling software lacks the capability to analyze short (<100bp) sequences. 
To bypass this hiatus, the short PCR amplicons were be pooled and ligated randomly to artifically create longer sequencable fragments.
This repository contains the code to identify, split and allocate the concatenated reads into to subreads and the subsequent SNP detection. 

## Method

Part 1: Splicing reads into subreads based on identified primers
        A fuzzy regex is used allowing up to 3 sequencing errors to occur at a random position in the primer sequence. 

Part 2: Variant identification 
        The subreads are mapped against a reference database of all loci consisting of the SNP and 25 nucleotides of flanking region on         either side.
        The aligment data is used to generate a table of all nucleotide variants and create a SNP profile
        
## Links

Oxford Nanopore Technologies
https://nanoporetech.com/

Metrichor
https://metrichor.com/

European Nucleotide Archieve
https://www.ebi.ac.uk

## Contributors

Senne Cornelis,
Yannick Gansemans,
Sander Willems,
Christophe Van Neste,
Filip Van Nieuwerburgh 
