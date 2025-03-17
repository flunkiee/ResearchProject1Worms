#the bit that calls all prereqs
import os
import re
import pandas as pd
from collections import defaultdict

def parse_gtf(file_path):
    #the bit that makes a big dictionary of gene ids
    gene_exons = defaultdict(int)
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith("#"):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != "exon":
                    continue

                attributes = fields[8]
                gene_id_match = re.search(r'gene_id "(.*?)"', attributes)
                if gene_id_match:
                    gene_id = gene_id_match.group(1)
                    gene_exons[gene_id] += 1
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
    
    return gene_exons

def calculate_averages(gene_exons):
    #the bit that calculates a bunch of averages
    total_exons = sum(gene_exons.values())
    total_genes = len(gene_exons)

    avg_exons = total_exons / total_genes
    avg_introns = (total_exons - total_genes) / total_genes

    return avg_exons, avg_introns

def main():
    #the bit that makes me wish i did this with a bit that goes through all the files in a dir and looks for gtf files instead
    gtf_files = [
    "A.australiensis.gtf",
    "C_formosanus.gtf",
    "C_fukii.gtf",
    "G_aquaticus.gtf",
    "G_montsenyensis.gtf",
    "G_parenensis.gtf",
    "N_munidae.gtf",
    "P_varius.gtf"
]

    results = []

    for gtf_file in gtf_files:
        if not os.path.exists(gtf_file):
            print(f"File not found: {gtf_file}")
            continue

        gene_exons = parse_gtf(gtf_file)
        avg_exons, avg_introns = calculate_averages(gene_exons)

        results.append({
            "File": gtf_file,
            "Average Exons per Gene": avg_exons,
            "Average Introns per Gene": avg_introns
        })

    #the bit that throws it into a dataframe
    results_df = pd.DataFrame(results)
    print(results_df)
#the bit that calls main
main()
