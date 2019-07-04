#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 18:37:42 2019

@author: veredas
"""

#!/usr/bin/env/python

import sys
import pandas

## First of all the functions needed are imported.

# 'read_fasta' takes a string which is the name of the FASTA file which contains the reads obtained from the sequencing.
# The function returns a dictionary: the keys are the head and the values the sequence of each read. 
def read_fasta (fasta_file):
    multifasta={}
    for lines in open(fasta_file, "r"):
        lines = lines.rstrip()
        if lines.startswith('>'):
            head = lines[1:]
            multifasta[head] = ''
        else:
            multifasta[head] += lines
    return multifasta

def read_fasta_gen (fasta_file):
    fasta = open(fasta_file, "r")
    line = fasta.readline()
    head = ""
    sequence = ""
    while line:
        line = line.rstrip("\n")
        if ">" in line:
            head = line
        else:
            sequence = sequence+line
        line = fasta.readline()
    return head, sequence


# 'read_bed' takes a string which is the name of the BED file which contains the chromosome, init and end position of each read in the chromosome, 
    #the ID of the read, the quality and the strand (+: forward, -:reverse) of the read.
# The function returns a list by row of the BED file.
def read_bed (bed_file):
    content = []
    with open(bed_file)as f:
        for line in f:
            content.append(line.strip().split())
    return content


# 'pos_reads' takes 5 arguments:
    # 1) A string (the sequence) of the gene studied.
    # 2) The lists with the rows of the bed file.
    # 3) The dictionary of the fasta file with the information of the reads
    # 4) The gene's first position in the choromosome (as a number).
    # 5) The gene's end position in the choromosome (as a number).
# The function returns a dictionary: the keys are the position of the gene and the values the bases and the number of each base in
    #that position of the gene. 
def pos_reads(gene, bed_file, fasta_file, pos_init, pos_end):
    dicc = {}
    value=[]
    diferents=[]
    for j in range(len(gene)):
        diferents = [] 
        for i in bed_file[0:len(bed_file)-1]:
            if i[5] == '+':
                init = int(i[1]) - pos_init
                end = int(i[2]) - pos_init - 1     
                long = list(range(init, end))       
                ID = i[3]        
                if j in long and end-init == len(fasta_file[ID])-1:
                    pos_read = j - init -1 
                    base_read = fasta_file[ID][pos_read]
                    value = [ID, pos_read, base_read]
                    diferents.append(value) 
            
            elif i[5] == '-':
                init = pos_end - int(i[1]) -1
                end = pos_end - int(i[2])     
                long = list(range(end, init))  
                ID = i[3]        
                if j in long and init-end == len(fasta_file[ID])-1:
                    pos_read = j - end-1
                    base_read = fasta_file[ID][pos_read]
                    value = [ID, pos_read, base_read]
                    diferents.append(value)                
        dicc[j] = diferents
        
    cuant={}
    for i in dicc:
        alelos = []
        count_A = 0
        count_C = 0
        count_G = 0
        count_T = 0
        for j in range(len(dicc[i])):
            if len(dicc[i][j]) > 0:
                base = dicc[i][j][2]
                if base == "A":
                    count_A += 1
                elif base == "T":
                    count_T += 1
                elif base == "C":
                    count_C += 1
                elif base == "G":
                    count_G += 1
            alelos = ["A:", count_A, "T:", count_T, "C:", count_C, "G:", count_G]
        cuant[i] = alelos
    return cuant


# 'pos_reads_pathogenic' takes 3 arguments:
    # 1) A string with the name of the downloaded or filtered CSV file  from Ensembl, which contains the positions in the
    #chromosome of the pathogenic variants in the studied gene
    # 2) A dictionary which keys are the positions of the gene and the values the number of each base which is located in that position in each read. 
# The function returns a dictionary which keys are the positions in the gene which are pathogenic (the values are the bases 
    #at those positions and the number of them in each position, like in the 'pos_reads' fucntion). 
def pos_reads_pathogenic (csv_pat_Ensembl, pos_init, pos_reads):
    pos_pat = []
    for i in csv_pat_Ensembl["Location"]:
        p = str(i)
        p = p.replace('1:', '')   
        pos = int(p) -pos_init
        pos_pat.append(pos)
    
    
    pat = {}
    for i in pos_reads:
        if i in pos_pat:
            pat[i] = pos_reads[i]
    return pat

# 'alleles_quantify' takes 2 arguments:
    # 1) A list with the positions in the gene which are compiled like variants 
    # 2) A dictionary which keys are the positions in the gene which are pathogenic and the values the number of each base which is located in that position in each read.
#The function returns a list with a sublist with the position of the pathogenic variant and the percentaje of reads for each allele.    
def alleles_quantify(pos_mut, reads_pat):
    mut_pat=[]
    reads_cuant=[]
    for i in pos_mut:
        j = i - 173872942
        if j in reads_pat.keys():
            mut = j+173872942
            mut_pat.append(mut)
            alleles = reads_pat[j]
            total = alleles[1] + alleles[3] + alleles[5] + alleles[7]
            for i in range(1, 8, 2):
                if alleles[i] != 0:
                    a_cuant = (alleles[i]*100)/total
                    allele = alleles[i-1]
                    lisst = [j+173872942, allele, a_cuant]
                    reads_cuant.append(lisst)
    return reads_cuant


# if the first argument is '-h', the help is shown
if sys.argv[1] == "-h":       
    print("Usage: SCRIPT_MOSAICISMO.py [reads_file.fa] [reads_file.bed] [gene_sequence.fa] [pathogenic_variants_Ensembl.csv] [gene start position in chromosome] [gene end position in chromosome] [gene_variants_positions.csv] [name of the CSV output file]")
    exit ()

# First argument: reads file name in fasta format
fasta_file = sys.argv[1]
# Second argument: reads file name in bed format
bed_file = sys.argv[2]
# Third argument: gene file name in fasta format
gene = sys.argv[3]
# Fourth argument: information of pathogenic variants in CSV format (downloaded from Ensembl)
pat_variants = sys.argv[4]
# Fifth: number of the chromosome position where the gene starts
init = int(sys.argv[5])
# Sixth: number of the chromosome position where the gene ends
end = int(sys.argv[6])
# Seventh: variants positions file name in CSV format
pos_variants_gene = sys.argv[7]
# Eighth: name of the output file 
output = sys.argv[8]
extension = output+".csv" #the output file will be a csv file

#Reads bed file
bed = read_bed(bed_file)
#Reads fasta file
fasta = read_fasta(fasta_file)
#Gene sequence fasta file
fasta37 = read_fasta_gen(gene)
gen = fasta37[1]
#Pathogenic variants from Ensembl
variantes = pandas.read_csv(pat_variants)
#Gene variants
var = pandas.read_csv(pos_variants_gene)
pos_var = []
for i in var.Pos:
    pos_var.append(i)
    

# First, the nucleotide by position in each read is calculated.
reads_pac = pos_reads(gen, bed, fasta, init, end)

# Second, the nucleotide by pathogenic position in each read is calculated.
reads_pat_pac = pos_reads_pathogenic(variantes, init, reads_pac)

#Third, the reads are quantifies by the percentage of the nucleotides in pathogenic positions.
quant_pac = alleles_quantify(pos_var, reads_pat_pac)

# The CSV output file is written
with open(extension, 'w') as f:

#Column names
    f.write("%s,%s,%s\n"%("Position", "Allele", "Percentage of reads")) 
    for i in quant_pac:
    
#column data
        f.write("%s,%s,%s\n"%(i[0], i[1], i[2])) 

