import escapecalculator
from Bio import SeqIO
import sys
import pandas as pd
def read_fasta(path):
    with open(path, 'r') as ffile:
        for record in SeqIO.parse(ffile, 'fasta'):
            return str(record.seq)

def get_ab_sources(lineage):
    ab_source_df = pd.read_csv('/home/jeliao/dms/antibody_sources.csv')
    antibodies = []
    added_abs = []
    #ab_lineage_df = ab_source_df['BA.1' in ab_source_df['source']]
    for index, row in ab_source_df.iterrows():
        if lineage in ab_source_df.loc[index, 'source'] and ab_source_df.loc[index, 'antibody'] not in antibodies:
            #antibodies.append(ab_source_df.loc[index, 'antibody'])
            if ab_source_df.loc[index, 'source'] not in added_abs:
                added_abs.append(ab_source_df.loc[index, 'source'])
    return added_abs
def translate(seq):
    bases = "tcag"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
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

def get_mut_sites(aa_seq):
    mut_sites = []
    aa_ref_spike = translate(read_fasta('/home/jeliao/dms/ref_spike.fasta'))
    #print(aa_ref_spike)
    for i in range(0, len(aa_seq)):
        if aa_seq[i:i+1] != aa_ref_spike[i:i+1]:
            mut_sites.append(i+1)
    return mut_sites

lineage = sys.argv[2]
ab_sources = get_ab_sources(lineage)
#calc = escapecalculator.EscapeCalculator(weight_by_neg_log_ic50 = False, study='any', virus='BA.1', sources = {"include_exclude": "include", "sources": ['BA.1 convalescents', 'BA.1 BTI', 'long-term BA.1 BTI', 'BA.1 BTI + BA.5/BF.7 infection', 'BA.1 + BA.5/BF.7 infection', 'long-term BA.1 convalescents', 'BA.1 convalescents reinfection']}) 
calc = escapecalculator.EscapeCalculator(weight_by_neg_log_ic50 = False, study='any', virus = lineage, sources = {"include_exclude": "include", "sources": ab_sources})
#mut_sites = [339, 371, 373, 375, 440, 446, 477, 478, 484, 493, 496, 498, 501, 505]
input_seq = sys.argv[1]
aa_seq = translate(read_fasta(input_seq))
all_mut_sites = get_mut_sites(aa_seq)
mut_sites = []
for site in all_mut_sites:
    if site > 318 and site < 542:
        mut_sites.append(site)
df = calc.escape_per_site(mut_sites).query("site in @mut_sites").round(3)
#print(df)
print('original ' + str(df['original_escape'].sum()))
print('retained ' + str(df['retained_escape'].sum()))
score = 0
#print(calc.escape_per_site([378])
"""
for site in mut_sites:
    single_site = [site]
    df = calc.escape_per_site([site]).query("site in @single_site").round(3)
    #print(df)
    score += df['retained_escape'].sum()
print(score)
"""
