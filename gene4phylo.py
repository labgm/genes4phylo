from Bio import SeqIO
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import argparse

parser = argparse.ArgumentParser(description="Generate separate directories for genes of interest for later phylogenetic analysis performance.")
parser.add_argument('--input', type=str, required=True, help='list of input genbank files.')
parser.add_argument('--genes', type=str, required=True, help='List of genes of interest.')
parser.add_argument('--output', type=str, required=True, help='Output directory.')

args = parser.parse_args()

genes = pd.read_csv('genes_of_interest.tsv',sep='\t')
genes_list = genes['Gene'].to_list()

with open(args.genes,'r') as genes:
    genes_list = genes.readlines()
    genes_list = [i.strip('\n') for i in genes_list]

with open(args.input,'r') as filelist:
    files = filelist.readlines()
    files = [i.strip('\n') for i in files]
    

if not os.path.exists(args.output):
    os.makedirs(args.output)

for file in files:
    filename = file.replace('.gb','')
    for record in SeqIO.parse(file, 'genbank'):
        for feature in record.features:
            if feature.type == 'CDS':
                if 'gene' in feature.qualifiers:
                    gene = feature.qualifiers['gene'][0]
                    if gene in genes_list:
                        if not os.path.exists(f'{args.output}/{gene}'):
                            os.makedirs(f'{args.output}/{gene}')
                            
                        start = feature.location.start
                        end = feature.location.end
                        cord = SeqFeature(FeatureLocation(start, end))
                        seq = cord.extract(record.seq)

                        ids = f'{gene}'
                        desc = f'Genome: {record.id}\tlocation: {feature.location}'

                        to_extract = SeqRecord(seq, id=ids, description=desc)
                            
                        output_file = f'{args.output}/{gene}/{gene}_{filename}.fasta'
                        with open(output_file, 'w') as outfile:
                            SeqIO.write(to_extract, outfile, 'fasta')
