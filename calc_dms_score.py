from Bio import SeqIO
import pandas as pd
import csv
import sys
import os
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

def get_escape(df_path):
    escape_score = 0
    prev_loc = ''
    df = pd.read_csv(df_path)
    for column in df.columns:
        if df[column].dtype == 'int64':
            df[column] = df[column].astype(str)
    for i in range(0, len(aa_seq)):
        for index, row in df.iterrows():
            if str(i+1) == df.loc[index, 'site']:
                if df.loc[index, 'site'] != prev_loc:

                    if aa_seq[i:i+1] == df.loc[index, 'wildtype']:
                        prev_loc = df.loc[index, 'site']
                    if aa_seq[i:i+1] == df.loc[index, 'mutant']:
                        escape_score += df.loc[index, 'escape_median']
                        #print(df.loc[index, 'site'] + df.loc[index, 'mutant'] + ' ' + str(df.loc[index, 'escape_median']))
                        prev_loc = df.loc[index, 'site']
                    #if aa_seq[i:i+1] != df.loc[index, 'wildtype']:
                    #    print(df.loc[index, 'site'] + df.loc[index, 'wildtype'] + aa_seq[i:i+1])
    return escape_score
def get_func_score(csv_path):
    func_score = 0
    df = pd.read_csv(csv_path)
    for i in range(0, len(aa_seq)):

        for index, row in df.iterrows():
            if str(i+1) == df.loc[index, 'sequential_site']:
                if aa_seq[i:i+1] == df.loc[index, 'mutant']:
                        escape_score += df.loc[index, 'effect']
    return func_score
seq_path = sys.argv[1]
lineage = sys.argv[2]
seq = read_fasta(seq_path)

aa_seq = translate(seq)
escape_score = 0
if lineage == 'BA.1':
    escape_dir = 'omicron_dms_data/escape_data'
    for file in os.listdir(escape_dir):
        file = os.path.join(escape_dir, file)
        if os.path.isfile(file):
            escape_score += get_escape(file)
            #print(file + ' ' + str(escape_score))
    #escape_score = get_escape('omicron_dms_data/escape_data/CC67.105_avg.csv') + get_escape('omicron_dms_data/escape_data/CC9.104_avg.csv') + get_escape('omicron_dms_data/escape_data/LyCoV-1404_avg.csv') + get_escape('omicron_dms_data/escape_data/NTD_5-7_avg.csv')
if lineage == 'B.1.617.2':
    escape_dir = 'delta_dms_data/escape_data'
    
    for file in os.listdir(escape_dir):
        file = os.path.join(escape_dir, file)
        if os.path.isfile(file):
            escape_score += get_escape(file)
            #print(file + ' ' + str(escape_score))
    
    #escape_score = get_escape('delta_dms_data/escape_data/267C_avg.csv') + get_escape('delta_dms_data/escape_data/279C_avg.csv') + get_escape('delta_dms_data/escape_data/REGN10933_avg.csv')

print(escape_score)
