#!/usr/bin/python
# BCR evaluator is used to evaluate bisulfite conversion ratio in WGSBS experiment
# Author : Quanyuan He Ph.D
# Email : hqyone@hunnu.edu.com
# Insititution: School of Medicine, Hunan Normal University
# Licensed under the MIT License
# Last update: 2019/02/11
# version: 1.3

import sys, getopt
import re , os
import gzip

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

# All possible telomeric block
tel_aseq_ls0 = ['TTTTAA','CTTTAA','TCTTAA','TTCTAA','CCTTAA','TCCTAA','CTCTAA','CCCTAA'] 
tel_aseq_ls1 = ['TTAAAA','TTAAAG','TTAAGA','TTAGAA','TTAAGG','TTAGGA','TTAGAG','TTAGGG'] 

# For test only
testTeloSeq = "CCCTAAGGTCCTAACTCTAACCTTAATTCTAAGGTTTTAATCTTAATTTTAA"

class ResultReport:
    def __init__(self, file):
        self.file= file
        self.type = "" #+/-
        self.min_unit_count= 0      # The mininal broke number
        self.total_reads = 0        # The total number of reads in the fastq file
        self.total_telo_reads = 0   # The total number of telomeric reads
        self.total_n3_reads = 0     # The number of c strand orignal N3 reads, which have at least one N3 block with three potential methylated sites
        self.total_unit_num = 0     # The total number of telomeric brock/units
        self.c_strand_unit_num = 0  # The number of brocks from C strands
        self.g_len_series = ""      # The histogram of brock number in g strand original sequences
        self.c_len_series = ""      # The histogram of brock number in c strand original sequences
        self.a_count_series = ""    # The statistics of all possible telomeric blocks from FASTQ_1
        self.b_count_series = ""    # The statistics of all possible telomeric blocks from FASTQ_2
        self.rc1 = float(0)         # The methylation ratio of site1
        self.rc2 = float(0)         # The methylation ratio of site2
        self.rc3 = float(0)         # The methylation ratio of site3

    def getTitle(self):
        return "file\t"+"type\t"\
            +"min_unit_count\ttotal_reads\t"\
            +"total_telo_reads\ttotal_n3_reads\ttotal_unit_num\t"\
            +"c_strand_unit_num\t"\
            +"g_len_series\t"\
            +"c_len_series\t"+",".join(tel_aseq_ls0)\
            +"\t"+",".join(tel_aseq_ls1)+"\t"+"rc1\t"+"rc2\t"+"rc3"
    
    def getValues(self):
        return self.file+"\t"\
                + str(self.type)+"\t"\
                + str(self.min_unit_count)+"\t"\
                + str(self.total_reads)+"\t"\
                + str(self.total_telo_reads)+"\t" \
                + str(self.total_n3_reads) + "\t" \
                + str(self.total_unit_num) + "\t" \
                + str(self.c_strand_unit_num) + "\t" \
                + self.g_len_series+"\t"\
                + self.c_len_series+"\t"\
                + self.a_count_series+"\t"\
                + self.b_count_series+"\t"\
                + str(round(self.rc1, 4))+"\t"\
                + str(round(self.rc2, 4))+"\t"\
                + str(round(self.rc3, 4))

# Patterns of C strand
p0 = '(([C,T]{3}TAA)+)'
# Patterns of G strand
p1 = '((TTA[G,A]{3})+)'

def findLongestTelomereMatch(seq, pattern):
    best_len = 0
    best_match_str = ""
    matchs = re.findall(pattern, seq)
    if matchs:
        for m in matchs:
            if len(m[0])>best_len:
                best_len = len(m[0])
                best_match_str = m[0]
    return best_match_str

def findBestTelomereMatch(seq):
    bs1= findLongestTelomereMatch(seq, p0)
    bs2 = findLongestTelomereMatch(seq, p1)
    if len(bs1)>=len(bs2):
        return [0,bs1]   # C strand pattern  If the reads can't match both of strands, just count as 0.
    else:
        return [1,bs2]   # G strand pattern

#a = findLongestTelomereMatch(testTeloSeq, p0)

def getTeloRepCount(min_unit_count, fastq_file):
    len_dic0 = [0]*30
    len_dic1 = [0]*30
    report = ResultReport(fastq_file)

    python_version = sys.version_info[0]

    report.min_unit_count = min_unit_count
    result_str =""

    # tel_aseq_ls0 = ['TTTTAA','CTTTAA','TCTTAA','TTCTAA','CCTTAA','TCCTAA','CTCTAA','CCCTAA'] 
    tel_seq_dic0 = {
        'CCCTAA': 0,
        'TCCTAA': 0,
        'CTCTAA': 0,
        'CCTTAA': 0,
        'TTCTAA': 0,
        'CTTTAA': 0,
        'TCTTAA': 0,
        'TTTTAA': 0
    }
    site1_ls0=[tel_aseq_ls0[1], tel_aseq_ls0[4], tel_aseq_ls0[6], tel_aseq_ls0[7]]
    site2_ls0=[tel_aseq_ls0[2], tel_aseq_ls0[4], tel_aseq_ls0[5], tel_aseq_ls0[7]]
    site3_ls0=[tel_aseq_ls0[3], tel_aseq_ls0[5], tel_aseq_ls0[6], tel_aseq_ls0[7]]
    

    # tel_aseq_ls1 = ['TTAAAA','TTAAAG','TTAAGA','TTAGAA','TTAAGG','TTAGGA','TTAGAG','TTAGGG'] 
    tel_seq_dic1 = {
        'TTAGGG': 0,
        'TTAGGA': 0,
        'TTAGAG': 0,
        'TTAAGG': 0,
        'TTAGAA': 0,
        'TTAAAG': 0,
        'TTAAGA': 0,
        'TTAAAA': 0
    }
    site1_ls1=[tel_aseq_ls1[1], tel_aseq_ls1[4], tel_aseq_ls1[6], tel_aseq_ls1[7]]
    site2_ls1=[tel_aseq_ls1[2], tel_aseq_ls1[4], tel_aseq_ls1[5], tel_aseq_ls1[7]]
    site3_ls1=[tel_aseq_ls1[3], tel_aseq_ls1[5], tel_aseq_ls1[6], tel_aseq_ls1[7]]

    unit_len = 6  # The lengtho of "TTAGGG"
    total_read_number = 0
    total_telo_number = 0
    c_n3_read_number = 0
    g_n3_read_number = 0

    #FASTQ = open(fastq_file,"r")
    # Go through fastq_file to classify sequences
    FASTQ=gzip.open(fastq_file, "rb")
    for line in FASTQ:
        if line.startswith("@".encode('utf-8')):
            total_read_number+=1
            n = 1
            if python_version>=3:
                line = FASTQ.readline().strip()
            else:
                line = FASTQ.next().strip()
            [type, bs] = findBestTelomereMatch(line)
            if type == 0:  # Seqeuence containing more (CCCTAA) blocks
                len_dic0[int(len(bs)/unit_len)]+=1
                if len(bs)/unit_len >= min_unit_count:
                    total_telo_number+=1
                    if "CCCTAA" in line:
                        c_n3_read_number+=1
                    for i in range(0, len(bs),unit_len):
                        if bs[i:i+unit_len] in tel_seq_dic0:
                            tel_seq_dic0[bs[i:i+unit_len]]+=1
                        else:
                            print (bs[i:i+unit_len])
            elif type == 1:  # Seqeuence containing more (TTAGGG) blocks
                len_dic1[int(len(bs)/unit_len)]+=1
                if len(bs)/unit_len >= min_unit_count:
                    total_telo_number+=1
                    if "TTAGGG" in line:
                        g_n3_read_number += 1
                    for i in range(0, len(bs), unit_len):
                        if bs[i:i+unit_len] in tel_seq_dic1:
                            tel_seq_dic1[bs[i:i+unit_len]]+=1
                        else:
                            print (bs[i:i+unit_len])

    report.c_len_series = ",".join(str(x) for x in len_dic0)
    report.g_len_series = ",".join(str(x) for x in len_dic1)
    report.total_reads = total_read_number
    report.total_telo_reads = total_telo_number

    #result_str += "len_dic0:\t"+"\t".join(str(x) for x in len_dic0)+"\n"
    #result_str += "len_dic1:\t"+"\t".join(str(x) for x in len_dic1)+"\n"

    temp_ls = []
    total_unit_num = 0
    for i in tel_aseq_ls0:
        total_unit_num += tel_seq_dic0[i]
        temp_ls.append(str(tel_seq_dic0[i]))

    report.a_count_series = ",".join(temp_ls)
    #result_str += "\t".join(tel_aseq_ls0)+"\n"
    #result_str += "\t".join(temp_ls)+"\n"

    temp_ls = []
    for i in tel_aseq_ls1:
        total_unit_num+=tel_seq_dic1[i]
        temp_ls.append(str(tel_seq_dic1[i]))

    report.b_count_series = ",".join(temp_ls)
    #result_str += "\t".join(tel_aseq_ls1)+"\n"
    #result_str += "\t".join(temp_ls)+"\n"
    
    
    #Summary for each sites
    #transformed_unit_num = total_unit_num-tel_seq_dic0['CCCTAA']-tel_seq_dic1['TTAGGG']
    c_strand_unit_num = 0
    c_strand_completed_trans_unit_num = 0  #All C are unmethylated
    site1=0
    site2=0
    site3=0

    #Decide what the strand represented dy the fastq file
    if (tel_seq_dic0[tel_aseq_ls0[7]]<tel_seq_dic1[tel_aseq_ls1[7]]):
        # if the number of 'CCCTAA' blocks < number of 'TTAGGG' blocks
        # The file is from FASTQ_1
        # So only count transformation in TTTTAA related reads
        report.type = "+"
        for key in tel_seq_dic0:
            c_strand_unit_num += tel_seq_dic0[key]
        c_strand_completed_trans_unit_num =  tel_seq_dic0['TTTTAA'] 
        for i in site1_ls0:
            if i in tel_seq_dic0:
                site1 += tel_seq_dic0[i]
        for i in site2_ls0:
            if i in tel_seq_dic0:
                site2 += tel_seq_dic0[i]
        for i in site3_ls0:
            if i in tel_seq_dic0:
                site3 += tel_seq_dic0[i]
    else:
        # The file is from FASTQ_2
        # So only count transformation in TTAAAA related reads
        report.type = "-"
        for key in tel_seq_dic1:
            c_strand_unit_num+=tel_seq_dic1[key]
        c_strand_completed_trans_unit_num =  tel_seq_dic1['TTAAAA']
        for i in site1_ls1:
            if i in tel_seq_dic1:
                site1 += tel_seq_dic1[i]
        for i in site2_ls1:
            if i in tel_seq_dic1:
                site2 += tel_seq_dic1[i]
        for i in site3_ls1:
            if i in tel_seq_dic1:
                site3 += tel_seq_dic1[i]
    if report.type=="+":
        report.total_n3_reads = c_n3_read_number
    else:
        report.total_n3_reads = g_n3_read_number
    report.total_unit_num = total_unit_num
    report.c_strand_unit_num = c_strand_unit_num
    report.rc1 = round(float(site1) / c_strand_unit_num, 4)
    report.rc2 = round(float(site2) / c_strand_unit_num, 4)
    report.rc3 = round(float(site3) / c_strand_unit_num, 4)

    """
    complete_transform_unit_num = tel_seq_dic0['TTTTAA']+tel_seq_dic1['TTAAAA']
    result_str += "Total_Unit_Number:\t"+str(total_unit_num)+"\n"
    result_str += "C_Stand_Unit_Number:\t"+str(c_strand_unit_num)+"\n"
    result_str += "C_Stand_Complete_Transform_Unit_num:\t"+str(c_strand_completed_trans_unit_num)+"\n"
    result_str += str(site1)+"\t"+str(site2)+"\t"+str(site3)+"\n"
    result_str += str(round(float(site1)/c_strand_unit_num,4))+"\t"+str(round(float(site2)/c_strand_unit_num,4))+"\t"+str(round(float(site3)/c_strand_unit_num,4))+"\n"
    """
    #print(report.getTitle())
    #print(report.getValues())
    return report
                
"""
# Inline testing part
findLongestTelomereMatch(testTeloSeq, p1)
n = 6
#rep_unit = "TTGCAA"
#fastq_file = "/home/server3/test.fastq.gz"
fastq_file="/home/server3/data2/zmq/ENCFF567DAI.fastq.gz"
result = getTeloRepCount(n,fastq_file)
print(result.getValues())

# python BCReval.py -n 8 -i /home/server3/data/zmq/ENCFF156UKB.fastq.gz -o ./FF156UKB_8.txt
"""

help_text="BCR evaluator is used to evaluate bisulfite conversion ratio in WGSBS experiment\n"\
          +"Author : Quanyuan He Ph.D\n"\
          +"Email : hqyone@hunnu.edu.com\n"\
          +"Last update: 2018/12/16\n"\
          +"version: 1.2\n"\
          +"\n"\
          +"Usage: BCReval.py -n <rep_number> -i <inputfiles> -o <outputfile>\n"\
          +"-n\t\tmin_rep_number\tmininal number of repeat units (default:6)','\n"\
          +"-i\t--ifile\tinnputfiles\tfastq files seperate by ','\n"\
          +"-o\t--ofile\toutputfile\toutput will be printed in the file or at console if it is null\n"


def main(argv):
   min_rep_number = 6
   inputfiles = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hn:i:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print ('BCReval.py -n <rep_number> -i <inputfiles> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print (help_text)
         sys.exit()
      elif opt in ("-n"):
          try:
              min_rep_number=int(arg)
          except:
              print (arg)
              print("min_rep_number should be a integer!")
              sys.exit()
      elif opt in ("-i", "--ifile"):
          inputfiles = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   fastq_files = inputfiles.split(",")
   OUT = open(outputfile, "w")
   if OUT:
       if outputfile!="" and  os.path.isfile(outputfile):
           a = ResultReport("")
           OUT.write(a.getTitle())
       for f in fastq_files:
           if os.path.isfile(f):
               rep = getTeloRepCount(min_rep_number, f)
               out_str = rep.getValues()
               if OUT:
                   OUT.write(out_str+"\n")
               else :
                   print (out_str+"\n")
           else :
               if f=="":
                   print ("Where are FASTQ files?")
               else:
                   print (f+" is not a file, omit.\n" )
   if OUT:
       OUT.close()

if __name__ == "__main__":
   main(sys.argv[1:])
