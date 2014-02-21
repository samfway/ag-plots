#!/usr/bin/env python

__author__ = "Sam Way"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Sam Way", "Rob Knight"] 
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Sam Way"
__email__ = "samfway@gmail.com"
__status__ = "Development"

import argparse
from numpy import asarray, array, zeros, argsort, cumsum, arange
from qiime.parse import parse_mapping_file_to_dict, parse_taxa_summary_table
import matplotlib.pyplot as plt
import brewer2mpl

def interface():
    args = argparse.ArgumentParser() 
    args.add_argument('-m', '--mapping-file', help='Mapping file', required=True)
    args.add_argument('-t', '--taxa-file', help='Taxa summary file', required=True)
    args.add_argument('-c', '--metadata-category', help='Metadata category', default='SIMPLE_MATTER')
    args.add_argument('-v', '--metadata-value', help='Specific metadata value', default=None)
    args.add_argument('-l', '--ylabel', help='Y-axis label', default='Phylum')
    args = args.parse_args()
    return args

def get_simple_taxa_summary(mapping_file, taxa_summary_file, metadata_category, \
    metadata_value, top_n_taxa=7, select_taxa=None): 
    """ Get a simplified taxonomy table. 

        Inputs:
        mapping_file - Input mapping file (sample ids match taxa file)
        taxa_summary_file - Taxonomy summary file 
        metadata_category - Used to select specific samples from the taxa file
        metadata_value - Value used to select specific samples from the taxa file
        top_n_taxa - If taxonomy groups aren't specified use the top N most abundant
        select_taxa - List of desired taxonomic groups

        Outputs:
        filtered_sample_ids - selected sample ids
        taxa_labels - taxonomic labels (including "Other")
        collapsed_taxa_table - simplified taxonomy table 
    """ 

    mapping_fp = open(mapping_file, 'rU')
    mapping_dict, comments = parse_mapping_file_to_dict(mapping_fp)
    taxa_fp = open(taxa_summary_file, 'rU')
    sample_ids, taxa_ids, taxa_table = parse_taxa_summary_table(taxa_fp)
    taxa_ids = [ taxa_id.split('__')[-1] for taxa_id in taxa_ids ] 

    selected_ids = [ key for key in mapping_dict.keys() if \
        mapping_dict[key][metadata_category] == metadata_value ] 
    sample_id_indices = [ i for i in xrange(len(sample_ids)) if sample_ids[i] in selected_ids ] 
    filtered_taxa_table = taxa_table[:, sample_id_indices]
    filtered_sample_ids = [ sample_ids[idx] for idx in sample_id_indices ]

    if select_taxa is None:
        # If select_taxa is None, take the top N most abundant 
        totals = filtered_taxa_table.sum(axis=1)
        taxa_indices = argsort(-totals)
        top_taxa = taxa_indices[:top_n_taxa]
        other_taxa = taxa_indices[top_n_taxa:]
        taxa_labels = [ taxa_ids[idx] for idx in top_taxa ]
    else:
        # List of taxa was supplied, use those 
        top_taxa = [ taxa_ids.index(x) for x in select_taxa ]
        other_taxa = [ t for t in xrange(len(taxa_ids)) if t not in top_taxa ] 
        taxa_labels = select_taxa  

    taxa_labels.append('Other')
    N = len(taxa_labels) # Number of classes/labels
    M = filtered_taxa_table.shape[1] # Number of samples after filtering

    # Sort samples by most_abundant_taxa
    sort_sample_indices = argsort(-filtered_taxa_table[top_taxa[0], :])
    filtered_taxa_table = filtered_taxa_table[:, sort_sample_indices]
    filtered_sample_ids = [ filtered_sample_ids[idx] for idx in sort_sample_indices ] 

    # Collapse "Others" rows into single row
    collapsed_taxa_table = zeros((N, filtered_taxa_table.shape[1]))
    collapsed_taxa_table[:-1, :] = filtered_taxa_table[top_taxa, :]
    collapsed_taxa_table[-1, :] = filtered_taxa_table[other_taxa, :].sum(axis=0) 
    total = collapsed_taxa_table.sum(axis=0)
    collapsed_taxa_table = collapsed_taxa_table / total

    return filtered_sample_ids, taxa_labels, collapsed_taxa_table 

def make_stacked_plot(filtered_sample_ids, taxa_labels, collapsed_taxa_table, \
    ylabel, colors, sample_ticks=None):
    """ Create a stacked plot. """ 
    N = len(taxa_labels)
    M = collapsed_taxa_table.shape[1]
    x = arange(M)
    cumulative = cumsum(collapsed_taxa_table, axis=0)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    # Get xticks, if any
    xticks = []
    xtick_labels = [] 
    if sample_ticks is not None:
        for sample_id, sample_label in sample_ticks:
            sample_index = filtered_sample_ids.index(sample_id)
            xticks.append(sample_index)
            xtick_labels.append(sample_label)

    ax1.fill_between(x, 1, 1-cumulative[0,:], color=colors[0])
    for i in xrange(1,N):
        ax1.fill_between(x, 1-cumulative[i-1,:], 1-cumulative[i,:], color=colors[i])
    yticks = arange(0.0, 1.25, .25)
    ytick_labels = [ str(int(ytick*100)) for ytick in yticks ] 
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xtick_labels, rotation=40, ha='right')
    ax1.xaxis.set_tick_params(width=1, length=10, pad=7, direction='out', top='off')
    ax1.set_yticks(yticks)
    ax1.set_yticklabels(ytick_labels)
    ax1.yaxis.set_tick_params(width=1, length=5, direction='out', right='off')
    plt.ylabel('%s Frequency (%%)' % (ylabel), fontsize=16)
    plt.ylim([0,1])
    plt.xlim([1,M])
    plt.subplots_adjust(bottom=0.25)
    plt.show()

def make_pie_chart(collapsed_taxa_table, colors):
    """ Create a simple pie chart """ 
    fractions = [ 100*x for x in collapsed_taxa_table.mean(axis=1) ] 
    wedges, texts = plt.pie(fractions, colors=colors)
    
    for w in wedges:
        w.set_linewidth(0)

    plt.axis('equal')
    plt.show()

if __name__=="__main__":
    args = interface() 

    special_labels = [ \
        ('000007120.1075677', 'Michael Pollan'), \
        ('000007108.1075657', 'Michael Pollan pre'), \
        ('000007118.1075682', 'Michael Pollan post'), \
        ('000007077.1075690', 'Jeff Leach'), \
        ('000007127.1075706', 'Shannon Ford') ]

    select_taxa = ['Firmicutes', 'Bacteroidetes', 'Proteobacteria', \
        'Verrucomicrobia', 'Actinobacteria', 'Tenericutes', 'Cyanobacteria']

    filtered_sample_ids, taxa_labels, collapsed_taxa_table = \
        get_simple_taxa_summary(args.mapping_file, args.taxa_file, \
        args.metadata_category, args.metadata_value, select_taxa=select_taxa)

    colors = brewer2mpl.get_map('Spectral', 'Diverging', len(taxa_labels)).mpl_colors

    make_stacked_plot(filtered_sample_ids, taxa_labels, collapsed_taxa_table, \
        args.ylabel, colors, sample_ticks=special_labels)

    make_pie_chart(collapsed_taxa_table, colors)
        
