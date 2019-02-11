# BCR_Evaluator
BCR evaluator is used to evaluate bisulfite conversion ratio in WGSBS experiment

Whole genome shotgun bisulfite sequencing (WGSBS) also known as BS-seq has been widely used to measure the methylation of whole genome at single-base resolution. One of the key steps in the assay is converting unmethylated cytosines into thymines (BS conversion). Incomplete conversion of unmethylated cytosines can introduces false positive methylation call. Developing a quick method to evaluate bisulfite conversion ratio (BCR) is benefit for both quality control and data analysis of WGSBS. Here BCReval is a small python script to estimate the unconverted rate (UCR) by using telomeric repetitive DNA as native spike-in control.

## Usage information
Usage: pythob BCReval.py -n <rep_number> -i <input_files> -o <output_file>

It has been tested in python 2.7
