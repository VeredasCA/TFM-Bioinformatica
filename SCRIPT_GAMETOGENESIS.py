#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 17:37:41 2019

@author: veredas
"""

import numpy
import sys

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
    for j in range(len(gene)):
        diferents = [] 
        for i in bed_file:
            if i[5] == '+':
                init = int(i[1]) - pos_init
                end = int(i[2]) - pos_init - 1    
                long = list(range(init, end))       
                ID = i[3]        
                if j in long and end-init == len(fasta_file[ID])-1:
                    pos_read = j - init -1 
                    base_read = fasta_file[ID][pos_read]
                    value = [ID, pos_read, base_read, i[5]]
                    diferents.append(value) 
            
            elif i[5] == '-':
                init = pos_end - int(i[1]) -1
                end = pos_end - int(i[2])     
                long = list(range(end, init))  
                ID = i[3]        
                if j in long and init-end == len(fasta_file[ID])-1:
                    pos_read = j - end-1
                    base_read = fasta_file[ID][pos_read]
                    value = [ID, pos_read, base_read, i[5]]
                    diferents.append(value)         
        #print("Count bases by reads")
        dicc[j] = diferents
    
    cuant={}

    for i in dicc:
        alelos = []
        count_A = 0
        count_C = 0
        count_G = 0
        count_T = 0
        count_r = 0
        count_f = 0
        for j in range(len(dicc[i])): 
            if len(dicc[i][j]) > 0:
                base = dicc[i][j][2]
                strand = dicc[i][j][3]
                if base == "A":
                    count_A += 1
                elif base == "T":
                    count_T += 1
                elif base == "C":
                    count_C += 1
                elif base == "G":
                    count_G += 1
                if strand == "+":
                    count_f += 1
                elif strand == "-":
                    count_r += 1
            alelos = ["A:", count_A, "T:", count_T, "C:", count_C, "G:", count_G, "-", count_r, "+", count_f]
            
            total = alelos[1] + alelos[3] + alelos[5] + alelos[7]
            
            if count_r > n_lec or count_f > n_lec:
                if (alelos[1]*100)/total > 79 or (alelos[3]*100)/total > 79 or (alelos[5]*100)/total > 79 or (alelos[7]*100)/total > 79:
                        continue
                else:                             
                    cuant[i] = alelos
            
    return cuant

# if the first argument is '-h', the help is shown
if sys.argv[1] == "-h":       
    print("Usage: SCRIPT_GAMETOGENESIS.py [reads_file_patient.fa] [reads_file_patient.bed] [gene_sequence.fa] [gene start position in chromosome] [gene end position in chromosome] [gene mutation position in chromosome] [minimun number of reads with the pathogenic variation] [reads_file_progenitor.fa] [reads_file_progenitor.bed] [name of the CSV output file]")
    exit ()

# First argument: reads file name in fasta format
fasta_file = read_fasta(sys.argv[1])
# Second argument: reads file name in bed format
bed_file = read_bed(sys.argv[2])
lec_pat = len(bed_file)
# Third argument: gene file name in fasta format
gene_file = read_fasta_gen(sys.argv[3])
gene=gene_file[1]
# Forth: number of the chromosome position where the gene starts
pos_init = int(sys.argv[4])
# Fifth: number of the chromosome position where the gene ends
pos_end = int(sys.argv[5])
# Sixth: number of chromosome position where the mutation is located
pos_mut = int(sys.argv[6])
# Seventh: cut-off of lectures
n_lec = int(sys.argv[7])
# Eigth argument: parental reads file name in fasta format
fasta_progenitor = read_fasta(sys.argv[8])
# Ninth argument: parental reads file name in bed format
bed_progenitor = read_bed(sys.argv[9])
lec_prog = len(bed_progenitor)
# Tenth: name of the output file 
output = sys.argv[10]
extension = output+".csv" #the output file will be a csv file



#Forward strand of the gene
gen_f=''
reverse = gene[::-1]
for i in reverse:
    if i == 'A':
        gen_f = gen_f+'T'
    elif i == 'T':
        gen_f = gen_f+'A'
    elif i == 'C':
        gen_f = gen_f+'G'
    elif i == 'G':
        gen_f = gen_f+'C'


#new bed: only contains the reads which contains the location of the pathogenic mutation
newbed=[]
for i in bed_file:
    init = int(i[1]) 
    end = int(i[2]) - 1  
    ID = i[3]
    if pos_mut > init and pos_mut < end and end-init == len(fasta_file[ID])-1:
        newbed.append(i)

newbed_progenitor=[]
for i in bed_progenitor:
    init = int(i[1]) 
    end = int(i[2]) - 1  
    ID = i[3]
    if pos_mut > init and pos_mut < end and end-init == len(fasta_progenitor[ID])-1:
        newbed_progenitor.append(i)


#Extract variants from the patient and the progenitor
print("Extracting variants positions in the patient...")
patient = pos_reads(gene, newbed, fasta_file, pos_init, pos_end)
print("Variants in reads with mutation in patient:", patient)

print("Extracting variants positions in the progenitor...")
progenitor = pos_reads(gene, newbed_progenitor, fasta_progenitor, pos_init, pos_end)
print("Variants in reads with mutation in progenitor:", progenitor)
 


#csv with variants
with open(output+"_DIC"+".csv", 'w') as a:
    
    #Column names
    a.write("%s,%s\n"%("variant position in gene", "alleles"))
    for i in patient:
        #column data
        a.write("%s,%s\n"%(i, patient[i]))
    a.write("%s,%s\n"%("progenitor variants", ""))
    for j in progenitor:
        #column data
        a.write("%s,%s\n"%(j, progenitor[j]))

# shared variants
variants_r_dic = {}
variants_f_dic = {}

for i in patient:
    if i in progenitor.keys():
        if patient[i][9] > n_lec:
            variants_r_dic[i] = patient[i][9]
        elif patient[i][11] > n_lec:
            variants_f_dic[i] = patient[i][11]
            
variants_r = list(variants_r_dic.keys())
variants_f = list(variants_f_dic.keys())
print("Shared variants in reverse strand:", variants_r, "\nShared variants in forward strand:", variants_f)


#csv with shared variants
with open(output+"_variants_shared"+".csv", 'w') as a:
    
    #Column names
    a.write("%s,%s,%s\n"%("variant position in gene", "reverse", "forward"))
    for i in variants_r:
        #column data
        a.write("%s,%s,%s\n"%(i, "rev", "-"))
    for j in variants_f:
        #column data
        a.write("%s,%s,%s\n"%(j, "-", "forward"))


## Find variants in normal or mutated allele
patogenica_f = gen_f[(pos_mut-pos_init)]
patogenica_r = gene[(pos_end-pos_mut)-1]
pos_variante_patogenica = pos_mut 
result_mut = []
result_nor = []
result_mut_f = []
result_nor_f = []
dic_mut = {}
dic_nor = {}
dic_mut_f = {}
dic_nor_f = {}
for i in newbed:
    if i[5] == "+":
        ID = i[3]
        posin = int(i[1])
        posend = int(i[2])
        seq = fasta_file[ID]
        variante_pato = seq[pos_variante_patogenica-posin-1]
        if variante_pato != patogenica_f:
            for j in variants_f_dic:
                pos_variante_global = pos_init + j
                pos_variante_read = pos_variante_global-posin-1
                if pos_variante_read < len(seq) and pos_variante_global > posin-1  and pos_variante_global < posend+1:
                    variante = seq[pos_variante_read]
                    pos_var_gen = gen_f[j]
                    if pos_var_gen != variante:
                        result_mut_f.append(j+pos_init)                        
                        print ("MUTATED ALLELE: Shared variant found (ref-variant in + strand):, {}-{} {} ({} chromosome position)\n".format(pos_var_gen,variante,j,pos_variante_global))
                        dic_mut_f[pos_variante_global] = [pos_var_gen, variante, len(newbed), len(newbed_progenitor)]
        else:
            for j in variants_f_dic:
                pos_variante_global = pos_init + j
                pos_variante_read = pos_variante_global-posin-1
                if pos_variante_read < len(seq) and pos_variante_global > posin-1  and pos_variante_global < posend+1:
                    variante = seq[pos_variante_read]
                    pos_var_gen = gen_f[j]
                    if pos_var_gen != variante:
                        result_nor_f.append(j+pos_init)
                        print ("NORMAL ALLELE: Shared variant found (ref-variant in + strand):, {}-{} {} ({} chromosome position)\n".format(pos_var_gen,variante,j,pos_variante_global))
                        dic_nor_f[pos_variante_global] = [pos_var_gen, variante, len(newbed), len(newbed_progenitor)]
    
    else:
        ID = i[3]
        posin = int(i[1])
        posend = int(i[2])
        seq = fasta_file[ID]
        variante_pato = seq[posend-pos_variante_patogenica-1]
        if variante_pato != patogenica_r:
            for j in variants_r_dic:
                pos_variante_global = pos_end - (j-1)
                pos_variante_read = posend-pos_variante_global
                if pos_variante_read < len(seq) and pos_variante_global > posin-1  and pos_variante_global < posend+1: 
                    variante = seq[pos_variante_read]
                    pos_var_gen = gene[j-1]
                    if pos_var_gen != variante:
                        result_mut.append(pos_end-(j-1))
                        print ("MUTATED ALLELE: Shared variant found (ref-variant in - strand):, {}-{} {} ({} chromosome position)\n".format(pos_var_gen,variante,j-1,pos_variante_global))
                        dic_mut[pos_variante_global] = [pos_var_gen, variante, len(newbed), len(newbed_progenitor)]
                        
        else:  
            for j in variants_r_dic:
                pos_variante_global = pos_end - (j-1)
                pos_variante_read = posend-pos_variante_global
                if pos_variante_read < len(seq) and pos_variante_global > posin-1 and pos_variante_global < posend+1: 
                    variante = seq[pos_variante_read]
                    pos_var_gen = gene[j-1]
                    if pos_var_gen != variante:
                        result_nor.append(pos_end-(j-1))
                        print ("NORMAL ALLELE: Shared variant found (ref-variant in - strand):, {}-{} {} ({} chromosome position)\n".format(pos_var_gen,variante,j-1,pos_variante_global))
                        dic_nor[pos_variante_global] = [pos_var_gen, variante, len(newbed), len(newbed_progenitor)]

final_f = numpy.unique(result_mut_f)    
print("There are {} shared variants (+) in mutated allele between patient and progenitor, in positions {}".format(len(final_f), final_f))
final_n_f = numpy.unique(result_nor_f)      
print("There are {} shared variants (+) in normal allele between patient and progenitor, in positions {}".format(len(final_n_f), final_n_f))

final = numpy.unique(result_mut)    
print("There are {} shared variants (-) in mutated allele between patient and progenitor, in positions {}".format(len(final), final))
final_n = numpy.unique(result_nor)      
print("There are {} shared variants (-) in normal allele between patient and progenitor, in positions {}".format(len(final_n), final_n))

 
#The CSV output file is written
with open(extension, 'w') as f:

#Column names
    f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%("Strand","Position_mutation", "Position_shared_variant_mutated_allele", "Position_shared_variant_normal_allele", "Allele _ref", "Allele_variant", "Reads_with_variant_in_mutated_allele", "Reads_with_variant_in_normal_allele", "Reads_with_pathogenic_patient", "Reads_with_pathogenic_progenitor", "Total_patient_reads", "Total_progenitor_reads", "Percentage")) 
#column data
    for a in dic_mut_f:     
        contador = result_mut_f.count(a)          
        f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%("forward",pos_mut, a, "-", dic_mut_f[a][0], dic_mut_f[a][1], contador, "-", dic_mut_f[a][2], dic_mut_f[a][3], lec_pat, lec_prog )) 

    for b in dic_nor_f:
        contador = result_nor_f.count(b)
        f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%("forward","-", "-", b, dic_nor_f[b][0], dic_nor_f[b][1], "-", contador, dic_nor_f[b][2], dic_nor_f[b][3], lec_pat, lec_prog)) 

    for i in dic_mut:
        contador = result_mut.count(i)
        f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%("reverse",pos_mut, i, "-", dic_mut[i][0], dic_mut[i][1], contador, "-", dic_mut[i][2], dic_mut[i][3], lec_pat, lec_prog )) 

    for j in dic_nor:
        contador = result_nor.count(j)
        f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%("reverse","-", "-", j, dic_nor[j][0], dic_nor[j][1], "-", contador, dic_nor[j][2], dic_nor[j][3], lec_pat, lec_prog)) 
