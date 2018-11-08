#ORFDiscovery

import numpy as np
import matplotlib as plt
import pandas as ps
import collections

from DNA_dictionaries import *



#Opening FASTA file

def open_fasta(filename):
    
    with open(filename) as f:
        
        description = f.readline()
        
        dna_seq_raw = f.read()
        
        dna_seq = dna_seq_raw.replace('\n','')
        
    return {'info': description, 'data': dna_seq}
    
    
    
####################
#Sense strand function 

def get_sense_strand(dna_seq):
    
    template_seq = ''
    
    final_template_seq = ''
    
    length_dna = len(dna_seq)
    
    for i in range(0,length_dna):
        
        letter = dna_seq[i]
        
        template_letter = DNA_comp_seq[letter]
        
        template_seq = template_seq + template_letter
        
    for i in range(0,len(template_seq)):
        
        base = template_seq[i]
        
        template_rna = transcription[base]
        
        final_template_seq = final_template_seq + template_rna.upper()
        
    return final_template_seq
    
############
#Antisense strand function

def get_antisense_strand(dna_seq):

    comp_seq = ""

    final_RNA_sequence = ""

    length_dna=len(dna_seq)

    for i in range(0,length_dna):

        letter=dna_seq[i]

        comp_letter=DNA_comp_seq[letter]

        comp_seq=comp_seq+comp_letter

        reversed_seq=comp_seq[::-1]

    length_rna = len(reversed_seq)

    for i in range(0,length_rna):

        rna_base = reversed_seq[i]

        rna_seq = rna_bases[rna_base]

        final_RNA_sequence = final_RNA_sequence + rna_seq

    return final_RNA_sequence.upper()

############

def get_codon(sequence,position):
    
    codon = ""
    
    if position <= len(sequence)-3:
        
        codon = sequence[position:position+3]
        
    return codon
    
##########

    
def rna_to_codon_sFrame(rna_seq, frame_number):
    
    #get the codon sequence from a Specific Frame number

    codon_list = [ ]

    start_pos = frame_number-1

    pos = start_pos

    while pos<len(rna_seq):

        codon = get_codon(rna_seq,pos)

        if codon!= "":

            codon_list.append(codon)

        pos = pos+3

    return codon_list
    
###########

def rna_to_codon_aFrame(rna_seq):
    
    #get codon sequences for All three frames
    
    codon_list_first = []
    codon_list_second = []
    codon_list_third = []
    
    pos = [0,1,2]
    
    
    while pos[0]<len(rna_seq):
        
        codon = get_codon(rna_seq, pos[0])
        
        if codon != "":
            
            codon_list_first.append(codon)
            
        pos[0] = pos[0] + 3
        
    while pos[1]<len(rna_seq):
        
        codon = get_codon(rna_seq, pos[1])
        
        if codon != "":
            
            codon_list_second.append(codon)
            
        pos[1] = pos[1] + 3
        
    while pos[2] < len(rna_seq):
        
        codon = get_codon(rna_seq, pos[2])
        
        if codon != "":
            
            codon_list_third.append(codon)
            
        pos[2] = pos[2] + 3
        
    return codon_list_first, codon_list_second, codon_list_third
    
###############



def codon_sequence_to_polypeptide(codon_list, frames = True):
    
    #This frames argument is very important here:
    #If one has used the rna_to_codon_sFrame() function, thereby returning
    #a codon list from one frame, then this code will return an amino acid list (a polypeptide sequence)
    #of that frame only and FRAMES FALSE
    
    #However, if one has used the rna_to_codon_aFrame() function, thereby returning 
    #a codon list for all three frames, then this code will return an amino acid list of
    #all three frames. These three amino acid sequences can still individually be found
    #by using the subscripting method (i.e. [ 0 ] at the end of the function). 
    #THEREFORE FRAMES SHOULD BE TRUE
    
    if frames == True:
        
        polypeptide_first = ""
        polypeptide_second = ""
        polypeptide_third = ""
        
        if codon_list:
            
            for i in range(0,3):
                
                for j in range(0, len(codon_list[i])):
                    
                    codon_of_interest = codon_list[i][j]
                    
                    if codon_of_interest == None:
                        
                        break
                        
                    aa_of_interest = dict_codon_to_aa[codon_of_interest]
                    
                    if i == 0:
                        
                        polypeptide_first = polypeptide_first + aa_of_interest
                    
                    elif i == 1: 
                        
                        polypeptide_second = polypeptide_second + aa_of_interest
                        
                    elif i == 2:
                        
                        polypeptide_third = polypeptide_third + aa_of_interest
    
        return polypeptide_first, polypeptide_second, polypeptide_third
                        
    
    elif frames != True:
        
        polypeptide = "" 
        
        for k in range(0,len(codon_list)):
            
            codon_of_interest = codon_list[k]
            
            if codon_of_interest == None:
                
                break
                
            aa_of_interest = dict_codon_to_aa[codon_of_interest]
            
            polypeptide = polypeptide + aa_of_interest
            
        return polypeptide

#############################



#DNA_sequence to polypeptide sequence

def DNA_sequence_to_polypeptide(dna_seq, PRINT = True):
    
    #This code takes the dna sequence from the fasta file and returns the amino acid
    #sequence for all three frames, for both sense and antisense strands
    
    #The print function is just to see what the frames actually are and is not necessary
    
    #The sense amino acid sequences can be subscripted from the output: 
    #sequence_to_polypeptide(sequence, PRINT= False)[0]
    
    #From this we can then access all three frames by further subscripting: 
    #sequence_to_polypeptide(sequence, PRINT= False)[0][0] 
    #^ this for instance, is the first frame amino acid sequence of the sense strand
    
    #The antisense strand can be found by:
    
    #sequence_to_polypeptide(sequence, PRINT= False)[1]
    #and the second frame amino acid sequence of the antisense strand is:
    #sequence_to_polypeptide(sequence, PRINT= False)[1][1]
    
    
    #get RNA strands
    
    sense_strand = get_sense_strand(dna_seq)
    
    antisense_strand = get_antisense_strand(dna_seq)
    
    #sense codon sequences
    
    sense_codons = rna_to_codon_aFrame(sense_strand)
    
    #antisense codon sequences
    
    antisense_codons = rna_to_codon_aFrame(antisense_strand)
    
    #sense codon polypeptides
    
    sense_polypeps = codon_sequence_to_polypeptide(sense_codons, frames = True)
    
    #antisense codon polypeptides
    
    antisense_polypeps = codon_sequence_to_polypeptide(antisense_codons, frames = True)
    
    if PRINT == True:
        
        print("Sense polypeptides from sequence are: \n Frame1: {} \n Frame2: {} \n Frame3: {}".format(sense_polypeps[0],
                                                                                                       sense_polypeps[1],
                                                                                                       sense_polypeps[2]))
    
    
        print("Antisense polypeptides from sequence are: \n Frame1: {} \n Frame2: {} \n Frame3: {}".format(antisense_polypeps[0],
                                                                                                        antisense_polypeps[1],
                                                                                                        antisense_polypeps[2]))
    
    sense_polypep_list = []
    sense_polypep_list.append(sense_polypeps[0])
    sense_polypep_list.append(sense_polypeps[1])
    sense_polypep_list.append(sense_polypeps[2])
    
    
    antisense_polypep_list = []
    antisense_polypep_list.append(antisense_polypeps[0])
    antisense_polypep_list.append(antisense_polypeps[1])
    antisense_polypep_list.append(antisense_polypeps[2])
    
    
    
    return sense_polypep_list, antisense_polypep_list
    
####################


def ORF_presence(polypep):
    
    #Will find the ORFs for a given polypeptide
    
    # can return more than one orf seq
    orf_seq_list=[]
    orf_seq_lengths = []

    for i in range(len(polypep)):

        current_pep1 = polypep[i]
        
        # check if we have start if so read to the next end pos

        if  current_pep1=='M': 

            for j in range(i,len(polypep)):

                current_pep2 = polypep[j]
                
                # check if we have end if so we have an orf

                if current_pep2=='*': 

                    orfseq=polypep[i:j+1]
    
                    orf_seq_list.append(orfseq)
                    
                    orf_seq_lengths.append(len(orfseq))

                    break    
    
    return orf_seq_list
    

#################################

def ORF_multiple_sense(polypep_list):
    
    #to be used after using the DNA_sequence_to_polypeptide() function
    #first returned object is sense, second is antisense
    
    ORF_sense_list = []
    ORF_antisense_list = []
    
    sense_polypeps = polypep_list[0]
    
    antisense_polypeps = polypep_list[1]
    
    for m in range(0, len(sense_polypeps)):
        
        sORF = ORF_presence(sense_polypeps[m])
        
        ORF_sense_list.append(sORF)
        
    for n in range(0, len(antisense_polypeps)):
        
        asORF = ORF_presence(antisense_polypeps[n])
        
        ORF_antisense_list.append(asORF)
        
        
    ORF_sense_list = [x for y in ORF_sense_list for x in y]
    ORF_antisense_list = [x for y in ORF_antisense_list for x in y]
    
    return ORF_sense_list, ORF_antisense_list
    
    
    
#####################################




def ORF_lengths(ORFS, counter = True):
    
    ORF_lengths = []
    
    for i in ORFS:
        
        ORF_lengths.append(len(i))
        
    if counter != True:
        
        return(ORF_lengths)
    
    else:
        
        lengthCount = collections.Counter(ORF_lengths)
        length, count = zip(*lengthCount.items())
        
        return length, count
        

        
#################        



#Thresholds definition

##### Create a class to enable different threshold lengths

class thresholds():
    
    a = 10
    b = 25
    c = 50
    d = 100
    e = 200
    f = 300
    
    def threshold_print(self):
        print("The amino acid (AA) thresholds used are: {}, {}, {}, {}, {}, {} AAs".format(self.a,self.b,self.c,self.d,self.e,self.f))
        
        
        
        
######################################
mythresh = thresholds()

######################################


def measuringORFs(ORF_lengths, dataframe = True):
    
    ORF_MORF_less_a = []
    ORF_MORF_above_a = []
    ORF_MORF_b = []
    ORF_MORF_c = []
    ORF_MORF_d = []
    ORF_MORF_e = []
    ORF_MORF_f = []
    
    for i in range(0,len(ORF_lengths)):
        
        if ORF_lengths[i] < mythresh.a:
            
            ORF_MORF_less_a.append(ORF_lengths[i])
        
        if ORF_lengths[i] >= mythresh.a:
            
            ORF_MORF_above_a.append(ORF_lengths[i])
            
        if ORF_lengths[i] >= mythresh.b:
            
            ORF_MORF_b.append(ORF_lengths[i])
            
        if ORF_lengths[i] >= mythresh.c:
            
            ORF_MORF_c.append(ORF_lengths[i])
            
        if ORF_lengths[i] >= mythresh.d:
            
            ORF_MORF_d.append(ORF_lengths[i])
            
        if ORF_lengths[i] >= mythresh.e:
            
            ORF_MORF_e.append(ORF_lengths[i])
            
        if ORF_lengths[i] >= mythresh.f:
            
            ORF_MORF_f.append(ORF_lengths[i])
            
    if dataframe == True: 
        
        MORF_list = [ORF_MORF_less_a, ORF_MORF_above_a, ORF_MORF_b, 
                         ORF_MORF_c, ORF_MORF_d, ORF_MORF_e, ORF_MORF_f]
        
        MORF_DF = ps.DataFrame(index = ['ORF < {}'.format(mythresh.a),'ORF > {}'.format(mythresh.a), 
                                    'ORF > {}'.format(mythresh.b),'ORF > {}'.format(mythresh.c),
                                    'ORF > {}'.format(mythresh.d),'ORF > {}'.format(mythresh.e), 
                                    'ORF > {}'.format(mythresh.f)], 
                           columns = ['Frequency'])
        
        for i in range(0, len(MORF_DF.index)):
        
            MORF_DF.iloc[i].Frequency = len(MORF_list[i])
            
        return(MORF_DF)
        
    else:
        
        return ORF_MORF_less_a, ORF_MORF_above_a, ORF_MORF_b, ORF_MORF_c, ORF_MORF_d, ORF_MORF_e, ORF_MORF_f
        
        