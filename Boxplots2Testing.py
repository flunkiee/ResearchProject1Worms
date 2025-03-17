#library stuff
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns #makes matplots look good.
import numpy as np
import re

def parse_gtf(file_path):
    #parse gtf create df 
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'] #default gtf columns
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=columns)
    return df

def calculate_lengths_and_stats(df):
    #calculate mean and median for introns and exons
    exon_lengths = []
    intron_lengths = []
    
    #grab geneid
    df['gene_id'] = df['attribute'].apply(
        lambda x: re.search('gene_id "([^"]+)"', x).group(1) if re.search('gene_id "([^"]+)"', x) else None #behold, regex for gene id again
    )

    #group by genid
    grouped = df.groupby('gene_id')
    
    for gene, group in grouped:
        #only exons
        exons = group[group['feature'] == 'exon']
        if exons.empty:
            continue
        
        #sort by start
        exons_sorted = exons.sort_values(by='start')
        
        #work out exon lengths
        exon_lengths.extend(exons_sorted['end'] - exons_sorted['start'] + 1) #as mentinoed in that other script, 1-based inclusive format, so this is done so its not offset by 1
        
        #introns calculated by the gaps between exons
        if len(exons_sorted) > 1:
            for i in range(1, len(exons_sorted)):
                intron_start = exons_sorted.iloc[i-1]['end'] + 1 #same thing as above
                intron_end = exons_sorted.iloc[i]['start'] - 1#same thing as above
                if intron_end >= intron_start:
                    length = intron_end - intron_start + 1#you get the idea
                    if length > 0:
                        intron_lengths.append(length)
    
    mean_exon_length = np.mean(exon_lengths) if exon_lengths else 0
    median_exon_length = np.median(exon_lengths) if exon_lengths else 0
    mean_intron_length = np.mean(intron_lengths) if intron_lengths else 0
    median_intron_length = np.median(intron_lengths) if intron_lengths else 0 #youd be surprised how often this bit is relevant. stops the program from exploding when theres an empty list 
    
    return exon_lengths, intron_lengths, mean_exon_length, median_exon_length, mean_intron_length, median_intron_length

def process_gtf_files(directory):
    #process all gtf junk
    exon_data = []
    intron_data = []
    #NO MORE LISTING ALL THE FILES MANUALLY!!
    for file_name in os.listdir(directory):
        if file_name.endswith('.gtf'):
            species_name = file_name.replace('.gtf', '')  # get rid of extension
            file_path = os.path.join(directory, file_name)
            print(f"Processing file: {file_name}")
            df = parse_gtf(file_path)
            exon_lengths, intron_lengths, mean_exon, median_exon, mean_intron, median_intron = calculate_lengths_and_stats(df)
            
            # add exon datapoints
            for length in exon_lengths:
                exon_data.append({'Species': species_name, 'Exon Length': length})
            
            # same for introns
            for length in intron_lengths:
                intron_data.append({'Species': species_name, 'Intron Length': length})
    
    exon_df = pd.DataFrame(exon_data)
    intron_df = pd.DataFrame(intron_data)
    return exon_df, intron_df

def filter_outliers_by_percentile(df, column, percentile=98):
#outlier filtering
    def filter_group(group):
        threshold = np.percentile(group[column], percentile)
        return group[group[column] <= threshold]
    
    filtered_df = df.groupby('Species', as_index=False).apply(filter_group, include_groups=False)
    filtered_df.reset_index(drop=True, inplace=True) #keeps spitting a deprecation warning but i did the thing it suggested how irritating
    return filtered_df

def plot_boxplots(exon_df, intron_df, percentile=98):
#actually plot the stuff
    exon_filtered = filter_outliers_by_percentile(exon_df, 'Exon Length', percentile)
    intron_filtered = filter_outliers_by_percentile(intron_df, 'Intron Length', percentile)
    
    sns.set_theme(style="whitegrid")
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    
    #exons
    ax_exon = axs[0]
    sns.boxplot(
        x='Species', y='Exon Length', hue='Species', data=exon_filtered,
        ax=ax_exon, showfliers=False
    )
    ax_exon.set_title(f'Exon Lengths per Species')
    ax_exon.tick_params(axis='x', rotation=45)
    if ax_exon.get_legend() is not None:
        ax_exon.get_legend().remove()
    
    #introns
    ax_intron = axs[1]
    sns.boxplot(
        x='Species', y='Intron Length', hue='Species', data=intron_filtered,
        ax=ax_intron, showfliers=False
    )
    ax_intron.set_title(f'Intron Lengths per Species')
    ax_intron.tick_params(axis='x', rotation=45)
    if ax_intron.get_legend() is not None:
        ax_intron.get_legend().remove()
    plt.tight_layout()
    plt.show()


#set dir here
directory_path = r'put your directory here' # oops left this as my directory originally 
    
    #call functinos
exon_df, intron_df = process_gtf_files(directory_path)
    
    #plot the stuff (call function)
plot_boxplots(exon_df, intron_df, percentile=98)

#if you made it here, apologies. 
#the plots this thing spits out look ok though! i probably shouldve made this two different scripts though.