# SNPore
Forensic Nanopore SNP sequencing

## Background

DNA analysis has become the cornerstone of contemporary forensic science. Altough most forensic DNA testing stil utilizing PCR and capillary electrophoresis (CE)-based analysis methods a transition to sequencing techniques is at hand.
Next generation targetted sequencing allows forensic scientists worldwide to harness the full potential of their DNA sample. 
Targetted sequencing allows casework and database efforts to be directed toward the genomic regions that best answer forensic questions

##Aim

The aim of the paper (Forensic SNP Analysis using Nanopore MinION™ Sequencing) associating this repository was to asses the capability of Oxford Nanopore Technologies’ (ONT) handheld sequencer MinION™ to created a forensic SNP based profile. 

## Motivation

The cloud based base calling software Metrichor being used default by ONT offers fully integrated, bespoke, real-time analysis solutions However, as ONT sequencing was primarerly designed to analyze long (>1kb) reads the base calling software lacks the possibility to analyze short (<100bp) fragments. 
To bypass this hiatus, the PCR amplicons were be pooled and ligated randomly to artifically create longer fragments. 
This repository contains the code to identify, split and allocate the concatenated reads into to subreads and the subsequent SNP detection. 

## Installation

Provide code examples and explanations of how to get the project.

## Contributors

Senne Cornelis,
Yannick Gansemans,
Sander Willems,
Filip Van Nieuwerburgh 
