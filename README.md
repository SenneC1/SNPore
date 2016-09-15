# SNPore
Forensic Nanopore SNP sequencing

## Synopsis

The aim of the paper (Forensic SNP Analysis using Nanopore MinION™ Sequencing) associating this repository was to asses the capability of Oxford Nanopore Technologies’ (ONT) handheld sequencer MinION™ to created a forensic SNP based profile. A 

## Motivation

The cloud based base calling software Metrichor being used default by ONT offers fully integrated, bespoke, real-time analysis solutions based on scalable, real time nanopore sensing technologies. As ONT sequencing was primarerly designed to analyze long (>1kb) reads the base calling lacks the possibility to analyze short (<100bp) fragments. 
To bypass this minimum fragment length requirement, the PCR amplicons have to be pooled and ligated randomly to create longer fragments. This repository contains the code to split the concatenated reads into to subreads and the subsequent SNP detection. 

## Installation

Provide code examples and explanations of how to get the project.

## Contributors

Senne Cornelis,
Yannick Gansemans,
Sander Willems,
Filip Van Nieuwerburgh 
