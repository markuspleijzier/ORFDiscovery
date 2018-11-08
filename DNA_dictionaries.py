#Dictionaries 

#complementary nucleotide dictionary
DNA_comp_seq = {"a":"t", "t":"a","g":"c","c":"g"}

#Transcription dictionary
transcription = {'a':'u','t':'a','g':'c','c':'g'}

#RNA base dictionary
rna_bases = {"a":"a", "t":"u", "g":"g", "c":"c"}

#Translation dictionary- RNA sequences to amino acids

dict_codon_to_aa = {

'UUU': 'F', 'UCU': 'S','UAU': 'Y','UGU': 'C',

'UUC': 'F', 'UCC': 'S','UAC': 'Y','UGC': 'C',

'UUA': 'L', 'UCA': 'S','UAA': '*','UGA': '*',

'UUG': 'L', 'UCG': 'S','UAG': '*','UGG': 'W',

'CUU': 'L', 'CCU': 'P','CAU': 'H','CGU': 'R',

'CUC': 'L', 'CCC': 'P','CAC': 'H','CGC': 'R',

'CUA': 'L', 'CCA': 'P','CAA': 'Q','CGA': 'R',

'CUG': 'L', 'CCG': 'P','CAG': 'Q','CGG': 'R',

'AUU': 'I', 'ACU': 'T','AAU': 'N','AGU': 'S',

'AUC': 'I', 'ACC': 'T','AAC': 'N','AGC': 'S',

'AUA': 'I', 'ACA': 'T','AAA': 'K','AGA': 'R',

'AUG': 'M', 'ACG': 'T','AAG': 'K','AGG': 'R',

'GUU': 'V', 'GCU': 'A','GAU': 'D','GGU': 'G',

'GUC': 'V', 'GCC': 'A','GAC': 'D','GGC': 'G',

'GUA': 'V', 'GCA': 'A','GAA': 'E','GGA': 'G',

'GUG': 'V', 'GCG': 'A','GAG': 'E','GGG': 'G',

}