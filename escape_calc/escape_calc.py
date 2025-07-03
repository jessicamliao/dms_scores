import csv
import escapecalculator
import pandas as pd
import tempfile
import sys
from Bio import SeqIO
import subprocess
from tempfile import NamedTemporaryFile
from pathlib import Path
import os
from io import StringIO
import glob
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

def run_cmd(cmd):
    #result = subprocess.run(cmd, shell=True,stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, capture_output=True)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(1)
    return result.stdout
def get_consensus(ref: str, seq: str, prefix: str):
    output_folder = 'consensus_fastas'
    os.makedirs(output_folder, exist_ok=True)
    spike_outputs = 'spikes'
    os.makedirs(spike_outputs, exist_ok = True)
    sam_file = f"{prefix}_aln.sam"
    vcf_file = f"{prefix}_variants.call.vcf.gz"
    consensus_file = os.path.join(output_folder, f"{prefix}_consensus.fa")
    sorted_bam = f"{prefix}_aln_sorted.bam"
    run_cmd(
            f"minimap2 -ax asm5 {ref} {seq} | samtools sort -o {sorted_bam}"
            )
    #bam_file = f"{prefix}_aln.bam"
    #run_cmd(
            #f"samtools view -S -b {sam_file} > {bam_file}"
            #)
    #sorted_bam = f"{prefix}_aln_sorted.bam"
    #run_cmd(
            #f"samtools sort {bam_file} -o {sorted_bam}"
            #)
    
    run_cmd(
            f"bcftools mpileup -Ou -f {ref} {sorted_bam} | bcftools call -mv -Oz -o {vcf_file}"
            )
    run_cmd(
            f"bcftools index {vcf_file}"
            )
    run_cmd(
            f"cat {ref} | bcftools consensus {vcf_file} > {consensus_file}"
            )
    
    """
    dna_seq = read_fasta(f"{consensus_file}")
    spike = dna_seq[21562:25384]
    aa_seq = translate(spike)
    return aa_seq
    """
    spike_bam = f"{prefix}_spike.sorted.bam"
    ref_spike = '/home/jeliao/dms/ref_spike.fasta'
    run_cmd(
            f"minimap2 -ax asm5 {consensus_file} {ref_spike} | samtools sort -o {spike_bam}"
            )
    bedtools_result = run_cmd(
            f"bedtools bamtobed -i {spike_bam}"
            )
    #print(bedtools_result)
    
    bedtools_data = pd.read_csv(StringIO(bedtools_result), sep="\t", header=None,
                       names=["chrom", "start", "end", "name", "score", "strand"])
    #print(bedtools_data)
    starts = bedtools_data["start"].tolist()[0] + 1
    ends = bedtools_data["end"].tolist()[0]
    region = f"{starts}-{ends}"
    #print(region)
    output_file = os.path.join(spike_outputs, f"{prefix}_spike.fasta")
    run_cmd(
            f"samtools faidx {consensus_file} NC_045512.2:{region} > {output_file}"
            )
    spike_seq = read_fasta(f"{output_file}")
    aa_spike_seq = translate(spike_seq)
    for file in glob.glob("*.bam"):
        os.remove(file)
    #os.remove(bam_file)
    os.remove(vcf_file)
    #os.remove(sorted_bam)
    for file in glob.glob("*.csi"):
        os.remove(file)
    return aa_spike_seq

def get_abs(lineage):
    ab_source_df = pd.read_csv('/home/jeliao/dms/antibody_sources.csv')
    antibodies = []
    added_abs = []
    #ab_lineage_df = ab_source_df['BA.1' in ab_source_df['source']]
    for index, row in ab_source_df.iterrows():
        if lineage in ab_source_df.loc[index, 'source'] and ab_source_df.loc[index, 'antibody'] not in antibodies:
            antibodies.append(ab_source_df.loc[index, 'antibody'])
            if ab_source_df.loc[index, 'source'] not in added_abs:
                added_abs.append(ab_source_df.loc[index, 'source'])
    return antibodies
    #return ab_lineage_df['antibody'].unique()

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

def get_mut_sites(aa_seq):
    all_mut_sites = []
    aa_ref_spike = translate(read_fasta('/home/jeliao/dms/ref_spike.fasta'))
    for i in range(0, len(aa_seq)):
        if aa_seq[i:i+1] != aa_ref_spike[i:i+1]:
            all_mut_sites.append(i+1)
    
    mut_sites = []
    for site in all_mut_sites:
        if site > 330 and site < 532:
            mut_sites.append(site)
    return mut_sites

def calc_score(aa_seq, seq_name):
    lineage_df = pd.read_csv('/home/jeliao/dms/100_seq_nextclade.tsv', delimiter = '\t')
    for index, row in lineage_df.iterrows():
        if lineage_df.loc[index, 'seqName'] == seq_name:
            sublineage = lineage_df.loc[index, 'partiallyAliased']
            #print(sublineage)
            ab_lineages = ['BA.5', 'BA.2']
            for ab_lineage in ab_lineages:
                if sublineage == 'XDA':
                    current_ab_lineage = 'BA.2'
                if 'XBB' in sublineage:
                    current_ab_lineage = 'BA.2'
                if ab_lineage in sublineage:
                    current_ab_lineage = ab_lineage
                    break
            mut_sites = get_mut_sites(aa_seq)
            #print(mut_sites)
            viruses = ['BA.1', 'BA.2.75', 'BA.2.86', 'BA.2', 'BA.5', 'BQ.1.1', 'D614G', 'HK.3.1', 'JN.1', 'KP.2', 'KP.3', 'SARS', 'XBB.1.5', 'XBB']
            lineage_df['clade_display'] = lineage_df['clade_display'].str.extract(r'\((.*?)\)')
            clade = lineage_df.loc[index, 'clade_display']
            for virus in viruses:
                if sublineage == 'XDA':
                    lineage = 'XBB'
                if virus in clade or virus in sublineage:
                    lineage = virus
                    break
            ab_sources = get_ab_sources(current_ab_lineage)
            calc = escapecalculator.EscapeCalculator(weight_by_neg_log_ic50 = False, study='any', virus = lineage, sources = {"include_exclude": "include", "sources": ab_sources})
            df = calc.escape_per_site(mut_sites).query("site in @mut_sites").round(3)
            escape_score = df['retained_escape'].sum()
            score_dict = {}
            return clade, escape_score
input_fasta = '/home/jeliao/dms/gisaid_hcov-19_2025_06_23_19.fasta'
output_file = 'output.csv'
with open(output_file, 'w') as out_file:
    writer = csv.writer(out_file)
    writer.writerow(['sublineage', 'score'])
    for record in SeqIO.parse(input_fasta, 'fasta'):
        name = record.id
        name_parts = name.split("|")
        gisaid = name_parts[-2]
        with tempfile.NamedTemporaryFile(mode='w+', suffix=".fasta", delete=False) as tmp:
            SeqIO.write(record, tmp, "fasta")
            temp_path = tmp.name
        aa_seq = get_consensus('/home/jeliao/Desktop/sequence.fasta', temp_path, gisaid)
        score_dict = calc_score(aa_seq, name)
        writer.writerows([calc_score(aa_seq, name)])
    #print(calc_score(aa_seq, name))
