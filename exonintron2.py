
import pandas as pd

# the list of gtf files again
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

with open('intronexon2out.txt', 'w') as output_file:
    # loop over each file
    for gtf_file in gtf_files:
        try:
            #pd df
            df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

            #default gtf col names
            df.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
            exons = df[df['feature'] == 'exon'].copy()
            exons['length'] = exons['end'] - exons['start'] + 1 #1 based inclusive format so +1 fixes this
            exons['gene_id'] = exons['attribute'].str.extract('gene_id "([^"]+)"') #regex to grab gene id
            exonssorted = exons.sort_values(['gene_id', 'start'])
            introns = []
            for gene_id, group in exonssorted.groupby('gene_id'):
                for i in range(1, len(group)):
                    intron_length = group.iloc[i]['start'] - group.iloc[i-1]['end'] - 1 #anytime you see a random +1 or -1 assume its because of the format
                    if intron_length > 0:
                        introns.append(intron_length)

            #avg stat calcs
            mean_exon = exonssorted['length'].mean()
            median_exon = exonssorted['length'].median()
            mean_intron = pd.Series(introns).mean()
            median_intron = pd.Series(introns).median()

            #wall of outputs
            output_file.write(f"Results for {gtf_file}:\n")
            output_file.write(f"  Mean exon length: {mean_exon}\n")
            output_file.write(f"  Median exon length: {median_exon}\n")
            output_file.write(f"  Mean intron length: {mean_intron}\n")
            output_file.write(f"  Median intron length: {median_intron}\n")
            output_file.write("\n" + "-"*50 + "\n") #looks pretty cool. makes a bunch of dashes to split up outputs!

        except Exception as e:
            output_file.write(f"something broke when processing {gtf_file}: {e}\n") #this shouldnt run but it will
            output_file.write("\n" + "-"*50 + "\n")

