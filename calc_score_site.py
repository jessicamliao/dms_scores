import csv
import pandas as pd
import sys
from Bio import SeqIO
bases = "tcag"
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
def read_fasta(path):
    with open(path, 'r') as ffile:
        for record in SeqIO.parse(ffile, 'fasta'):
            return str(record.seq)

def translate(seq):
    seq = seq.lower().replace('\n', '').replace(' ', '')
    aa_seq = ''
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        #print(codon)
        aa = codon_table.get(codon, '*')
        if aa != '*':
            aa_seq += aa
        elif aa == '*':
            break

    return aa_seq

site_df = pd.read_csv('escape.csv')

ref_spike = 'ref_spike.fasta'
seq_path = sys.argv[1]
ab = sys.argv[2]
ab_df = site_df[site_df['antibody']==ab]

seq = read_fasta(seq_path)
aa_ref_spike = translate(read_fasta(ref_spike))
aa_seq = translate(seq)
escape_score = 0
current_site = ''
#print(mut_df['site'])
#print(len(aa_seq))
for i in range(0, len(aa_seq)):
    if aa_seq[i:i+1] != aa_ref_spike[i:i+1]:
        for idx, row in ab_df.iterrows():
            if ab_df.loc[idx, 'site'] == i+1:
                #print(ab_df.loc[idx, 'escape'])
                escape_score += ab_df.loc[idx, 'escape']
print(escape_score)
