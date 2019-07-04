#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 12:20:40 2019

@author: veredas
"""
import pandas
import sys

## First, the functions which are needed are imported.

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

# 'hotspots' takes 3 arguments:
    # 1) A string with the sequence of the gene which is studied.
    # 2) A string with the pattern which is serached in the gene sequence.
    # 3) A list with the positions of the nucleotides located in the exons of the gene. 
# The function returns two lists, the first one is the position of the nucleotides of the pattern which are located in 
    #the exons, and the second one is the positions in the introns. 
def hotspots (gene, sequence, exons, pos_init):
    len_seq = len(sequence)
    pos = []
    in_exons = []
    in_introns = []
    for i in range(len(gene) - len_seq + 1):
        if gene[i:i+len_seq] == sequence:
            for j in range(len_seq):
                pos.append(i+j+pos_init)
    for h in pos:
        if h in exons:
            in_exons.append(h)
        else:
            in_introns.append(h)
    return in_exons, in_introns

def hotspots_read (read, sequence, exons, pos_init):
    len_seq = len(sequence)
    pos = []
    in_exons = []
    in_introns = []
    for r in read:
        seq = read[r]
        for i in range(len(seq) - len_seq + 1):
            if seq[i:i+len_seq] == sequence:
                for j in range(len_seq):
                    pos.append(i+j+pos_init)
        for h in pos:
            if h in exons:
                in_exons.append(h)
            else:
                in_introns.append(h)
    return in_exons, in_introns


# 'forward_reverse' takes a list of positions in a chromosome. 
# The function returns a list with the positions in a chromosome but in the forward string 
def forward_reverse (positions, pos_init, pos_end):
    chain = []
    for p in positions:
        pos = (pos_end+1) - p
        pos = pos + (pos_init+1)
        chain.append(pos)
    return chain

# 'pathogenic' takes 2 arguments:
    # 1) A list with the positions in the chromosome of the nucleotides wich are in a hotspot. 
    # 2) A string with the name of the CSV downloaded or filetered from Ensemble which contains the information
    #of the pathogenic variants in the gene. 
# The function returns a list with the position of the nucleotides which are in a pathogenic location.
def pathogenic(position, variability_pat):
    csv_pat = pandas.read_csv(variability_pat)
    loc = csv_pat["Location"]
    locat=[]
    for i in loc:
        l = i.replace("1:", "")
        l = int(l)
        locat.append(l)
        
    is_pat=[]
    for p in position:
        if p in locat:
            p = str(p)
            #pos = "1:"+p
            is_pat.append(p)
    return is_pat


# if the first argument is '-h', the help is shown
if sys.argv[1] == "-h":       
    print("Usage: SCRIPT_HOTSPOTS.py [gen/reads.fa] [hotspot pattern] [init-end_exons.csv] [pathogenic_variants.csv] [gene start position in chromosome] [gene end position in chromosome] [output name of CSV file]")
    exit ()

# First argument: reads or gene file name in fasta format
fasta_file = sys.argv[1]
# Second argument: sequence of the hotspot
pattern = sys.argv[2]
# Third argument: position of exons file name in csv format
exons_input = sys.argv[3]
exons_f = pandas.read_csv(exons_input)
# Third argument: information of ensembl pathogenic variants in csv format
variability = sys.argv[4]
# Forth argument: gene start position in chromosome
pos_init = int(sys.argv[5])
# Forth argument: gene end position in chromosome
pos_end = int(sys.argv[6])
# Seventh: name of the output file 
output = sys.argv[7]
extension = output+".csv" #the output file will be a csv file



# First, the position of exons is calculated in the reverse strand. 
start=[]
end=[]
exons = []
for i in exons_f.start:
    start.append(i)
for i in exons_f.end:
    end.append(i)
for j in range(len(start)):
    exon = list(range(start[j], end[j]))
    exons = exons+exon

exons_r = forward_reverse(exons, pos_init, pos_end)


# Second, hotspot are searched in the sequence.
fasta = read_fasta(fasta_file)
if len(fasta) == 1: #there is only a sequence (probably it is a gene)
    fasta = read_fasta_gen(fasta_file)
    gen = fasta[1] #gen in the sequence, without head
    hotspot_pos = hotspots(gen, pattern, exons_r, pos_init)
else: 
    multifasta = read_fasta_gen(fasta_file)
    hotspot_pos= hotspots_read(multifasta, pattern, exons_r, pos_init)

    
#Third, the positions are calculated in the forward strand of the gene.
hotspot_ex = forward_reverse(hotspot_pos[0], pos_init, pos_end)
hotspot_in = forward_reverse(hotspot_pos[1], pos_init, pos_end)


#Forth, the pathogenic hotspots are calculated.
hotspot_ex_pat = pathogenic(hotspot_ex, variability)
hotspot_in_pat = pathogenic(hotspot_in, variability)

result = {"Exons": [len(hotspot_ex), hotspot_ex], 
          "Introns" :  [len(hotspot_in), hotspot_in], 
          "Exons pathogenic" :  [len(hotspot_ex_pat), hotspot_ex_pat],
          "Introns pathogenic" :  [len(hotspot_in_pat), hotspot_in_pat]}

# The CSV output file is written
with open(extension, 'w') as f:

#Column names
    f.write("%s,%s,%s\n"%("Type of pathogenic hot-spot", "Number of hot-spot", "Locations")) 
    for i in result:  
        for j in result[i][1]:
#column data
            f.write("%s,%s,%s\n"%(i, result[i][0], j)) 


