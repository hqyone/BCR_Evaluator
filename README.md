# BCR_Evaluator
BCR evaluator is used to evaluate bisulfite conversion ratio in WGSBS experiment

Whole genome shotgun bisulfite sequencing (WGSBS) also known as BS-seq has been widely used to measure the methylation of whole genome at single-base resolution. One of the key steps in the assay is converting unmethylated cytosines into thymines (BS conversion). Incomplete conversion of unmethylated cytosines can introduces false positive methylation call. Developing a quick method to evaluate bisulfite conversion ratio (BCR) is benefit for both quality control and data analysis of WGSBS. Here BCReval is a small python script to estimate the unconverted rate (UCR) by using telomeric repetitive DNA as native spike-in control.

## Usage information
Usage: pythob BCReval.py -n <rep_number> -i <input_files> -o <output_file>

The input_file should be a FASTQ file or a list of FASTQ files which are seperated by commas. 
It has been tested in python 2.7

## Output format
The output is a string including 15 fields which are seperated by tabs.

Example:
/home/server3/data2/zmq/ENCFF156UKB.fastq.gz	+	8	925473550	125151	0	2022759	90386	0,273921715,20071116,304864,21496,9547,9123,9530,9139,8617,8387,7941,7707,7244,6903,6392,5836,5323,4954,4450,4221,4119,3765,3524,16106,6095,0,0,0,0	83503771,542711711,4739860,42520,2172,544,250,180,176,108,111,107,111,110,96,76,105,102,109,163,163,210,194,200,1713,574,0,0,0,0	89547,290,300,244,0,4,1,0	304,76,2761,91,1383,4504,5226,1918028	0.0032	0.0034	0.0028

In which:
1. File: home/server3/data2/zmq/ENCFF329YKL.fastq.gz
2. Type : +             //from FASTQ_1(+) or FASTQ_2(-)
3. Min_unit_count : 8   //The minimal number of telomeric rep units
4. Total_reads: 925473550  //Total number of reads in the FASTQ file
5. Telomeric reads: 125151 . // Telomeric reads numbers
6. Total_n3_reads:0     //The total number of reads which contain N3 units that has three methylated cytosines 
7. total_unit_num:2022759 . //The number of telomeric units in the files (All reads)
8. c_strand_unit_num: 90386 . //The number of telomeric units of C-strand derived telomeric reads
9. g_len_series:0,273921715,20071116,....6095,0,0,0,0 . // The number of reads containing certain number of G-telomeric derived units (TTAGGG). For example here 273921715 is the number of reads containing one G-unit.
10. c_len_series:83503771,542711711,4739860,....574,0,0,0,0 . // The number of reads containing certain number of C-telomeric derived units (CCCTAA).
11. a_count_series:89547,290,300,244,0,4,1,0  //The numbers of all possible units (A series, from FASTQ_1)
12. b_count_series:304,76,2761,91,1383,4504,5226,1918028  //The numbers of all possible units (B series, from FASTQ_2), you can omit it here because this data is from FASTQ1 file.
13. rc1: 0.0032 . //The unconverted/methylated ratio of the first C
14. rc2: 0.0034 . //The unconverted/methylated ratio of the second C
15. rc3: 0.0028 . //The unconverted/methylated ratio of the third C
