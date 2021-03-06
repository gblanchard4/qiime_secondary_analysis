#!/usr/bin/env python

__author__ = "Gene Blanchard"
__email__ = "me@geneblanchard.com"

# Imports
import argparse
import sys
import subprocess
import datetime
import csv
from itertools import izip

'''
Python secondary analysis script

This script takes the basic primary analysis files (biom, map, tre, params) and runs the SOP commands
'''


# Function to make sure qiime is loaded in the PATH
def check_qiime(log):
    try:
        subprocess.Popen("print_qiime_config.py", stdout=log)

    except OSError:
        print "Qiime is not loaded into your path"
        sys.exit()


# Convert the map file to a dictionary, removing default categories
def map_to_dictionary(mapfile):
    dictionary = dict()
    with open(mapfile, "rb") as map_handle:
        mapcsv = izip(*csv.reader(map_handle, delimiter='\t'))
        for row in mapcsv:
            header = row[0]
            unique_values = list(set(row[1:]))
            unique_values.insert(0, header)
            dictionary[unique_values[0]] = unique_values[1:]
    try:
        del dictionary['#SampleID']
        del dictionary['Description']
        del dictionary['BarcodeSequence']
        del dictionary['LinkerPrimerSequence']
    except KeyError:
        pass
    return dictionary


# Save the biom table stats to a file
def create_biom_summary(biom_file, batch):
    command = "print_biom_table_summary.py -i {} -o table_summary.txt\n".format(biom_file)
    batch.write(command)


def summarize_taxa(biom, mapfile, params, categories, batch):
    basic_taxa_command = "summarize_taxa_through_plots.py -s -i {} -m {} -p {} -o taxa_summary/taxa_individual\n".format(biom, mapfile, params)
    batch.write(basic_taxa_command)
    # Heatmap PDF/PNG
    for level in ["L2", "L3", "L4", "L5", "L6", "L7"]:
        taxa_level = "otus_sorted_{}.biom".format(level)
        heatmap_pdf_command = "make_otu_heatmap.py -i taxa_summary/taxa_individual/{} -o taxa_summary/taxa_individual/{}.pdf\n".format(taxa_level, level)
        heatmap_png_command = "make_otu_heatmap.py -i taxa_summary/taxa_individual/{} -o taxa_summary/taxa_individual/{}.png -g png\n".format(taxa_level, level)
        batch.write(heatmap_pdf_command)
        batch.write(heatmap_png_command)
    for category in categories:
        category_taxa_command = "summarize_taxa_through_plots.py -i {} -m {} -p {} -o taxa_summary/taxa_{} -c {} -s\n".format(biom, mapfile, params, category, category)
        batch.write(category_taxa_command)
        # for level in ["L2", "L3", "L4", "L5", "L6", "L7"]:
        #    taxa_level = "{}_otu_table_sorted_{}.biom".format(category, level)
        #    heatmap_pdf_command = "make_otu_heatmap.py -i taxa_summary/taxa_{}/{} -o taxa_summary/taxa_{}/{}.pdf -m {} -c {}\n".format(category, taxa_level, category, level, mapfile, category)
        #    heatmap_png_command = "make_otu_heatmap.py -i taxa_summary/taxa_{}/{} -o taxa_summary/taxa_{}/{}.png -m {} -c {} -g png\n".format(category, taxa_level, category, level, mapfile, category)
        #    batch.write(heatmap_pdf_command)
        #    batch.write(heatmap_png_command)


def alpha_diversity(biom, mapfile, params, tre, batch):
    alpha_command = "alpha_rarefaction.py -i {} -m {} -p {} -t {} -a -O 24 -o alpha_diversity\n".format(biom, mapfile, params, tre)
    batch.write(alpha_command)


def compare_alpha_diversity(mapfile, categories, batch):
    alpha_diversity_path = 'alpha_diversity/alpha_div_collated/'
    metrics = ['chao1', 'observed_species', 'PD_whole_tree', 'shannon', 'simpson']
    for metric in metrics:
        for category in categories:
            metric_file = alpha_diversity_path + metric
            output = alpha_diversity_path + metric + '_' + category
            compare_alpha_diversity_command = "compare_alpha_diversity.py -m {} -n 9999 -c {} -i {}.txt -o {} \n".format(mapfile, category, metric_file, output)
            batch.write(compare_alpha_diversity_command)


def beta_diversity(biom, mapfile, tre, batch):
    beta_diversity_command = "beta_diversity_through_plots.py -i {} -m {} -t {} -o beta_diversity\n".format(biom, mapfile, tre)
    batch.write(beta_diversity_command)
    # Make the 2d plots
    make_2d_unweighted_command = "make_2d_plots.py -i beta_diversity/unweighted_unifrac_pc.txt -m {} -o beta_diversity/2d_unweighted_unifrac_plots\n".format(mapfile)
    make_2d_weighted_command = "make_2d_plots.py -i beta_diversity/weighted_unifrac_pc.txt -m {} -o beta_diversity/2d_weighted_unifrac_plots\n".format(mapfile)
    batch.write(make_2d_unweighted_command)
    batch.write(make_2d_weighted_command)


def compare_beta(mapfile, categories, batch):
    for category in categories:
        compare_beta_unweighted_command = "compare_categories.py --method anosim -i beta_diversity/unweighted_unifrac_dm.txt -m {} -c {} -n 9999 -o beta_diversity/ANOSIM_{}_unweighted\n".format(mapfile, category, category)
        compare_beta_weighted_command = "compare_categories.py --method anosim -i beta_diversity/weighted_unifrac_dm.txt -m {} -c {} -n 9999 -o beta_diversity/ANOSIM_{}_weighted\n".format(mapfile, category, category)
        batch.write(compare_beta_unweighted_command)
        batch.write(compare_beta_weighted_command)


def compute_core_microbiome(biom, mapfile, categories, batch):
    state_dictionary = map_to_dictionary(mapfile)
    for category in categories:
        for state in state_dictionary[category]:
            core_microbiome_command = "compute_core_microbiome.py -i {} --mapping_fp {} --num_fraction_for_core_steps 6 -o core_microbiome/core_microbiome_{}_{} --valid_states {}:{}\n".format(biom, mapfile, category, state, category, state)
            batch.write(core_microbiome_command)


def qiime17_otu_category_sig(biom, mapfile, categories, batch):
    tests = ['g_test', 'ANOVA']
    for category in categories:
        for test in tests:
            otu_cat_sig_command = "otu_category_significance.py -i {} -m {} -c {} -s {} -o taxa_summary/taxa_{}/{}.txt".format(biom, mapfile, category, test, category, test)
            batch.write(otu_cat_sig_command)


def get_time():
    now = datetime.datetime.now().ctime()
    return str(now) + '\n'


def hms_string(sec_elapsed):
    h = int(sec_elapsed / (60 * 60))
    m = int((sec_elapsed % (60 * 60)) / 60)
    s = sec_elapsed % 60.
    return "{}:{:>02}:{:>05.2f}".format(h, m, s)


def group_significance(biom, mapfile, categories, batch):
    tests = ['nonparametric_t_test', 'bootstrap_mann_whitney_u', 'ANOVA', 'kruskal_wallis', 'g_test', 'parametric_t_test', 'mann_whitney_u']
    for category in categories:
        for test in tests:
            otu_group_significance_command = "group_significance.py -i {} -m {} -c {} -s {} -o taxa_summary/{}_{}.txt\n".format(biom, mapfile, category, test, category, test)
            batch.write(otu_group_significance_command)


def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="This script runs through a standard QIIME secondary analysis pipeline. The required input files are the biom, map, tre, and params. ")

    # biom -b --biom
    parser.add_argument("-b", "--biom", dest="biom", required=True, help="The biom file")
    # map -m --map
    parser.add_argument("-m", "--map", dest="map", required=True, help="The mapping file")
    # parameters -p --params
    parser.add_argument("-p", "--params", dest="params", required=True, help="The parameters file")
    # tre -t --tre
    parser.add_argument("-t", "--tre", dest="tre", required=True, help="The tre file")
    # categories -c --categories
    parser.add_argument("-c", "--categories", dest="categories", help="The metadata categories to compute. Must be colon seperated")

    # Parse the arguments
    args = parser.parse_args()

    # Assign variables
    biom_file = args.biom
    mapping_file = args.map
    params_file = args.params
    tre_file = args.tre

    valid_params = """
summarize_taxa:level 2,3,4,5,6,7
plot_taxa_summary:labels Phylum,Class,Order,Family,Genus,Species
alpha_diversity:metrics shannon,simpson,PD_whole_tree,chao1,observed_species
multiple_rarefactions:min 100
multiple_rarefactions:max 18000
multiple_rarefactions:step 500
beta_diversity_through_plots:seqs_per_sample 18000
"""

    # Make sure a valid parameters file has been input
    required_values = ['summarize_taxa:level', 'plot_taxa_summary:labels', 'alpha_diversity:metrics', 'multiple_rarefactions:min', 'multiple_rarefactions:max', 'multiple_rarefactions:step', 'beta_diversity_through_plots:seqs_per_sample']
    found_values = []
    with open(params_file, 'r') as params:
        for line in params:
            value = line.split(' ')[0]
            found_values.append(value)
    for required in required_values:
        if required not in found_values:
            print required + " not found in parameters file.\n"
            print valid_params
            sys.exit()

    # Read the mapping file into a dictionary
    # The key is the column header (category)
    categories_dictionary = map_to_dictionary(mapping_file)
    # Check if the user has inputed their own categories
    # If no categories have been selected, use all categories in the categories_dictionary
    if args.categories is None:
        categories = categories_dictionary.keys()
    # User selected categories
    else:
        categories = []
        # Split the categories into list and make sure they are valid
        for category in args.categories.split(':'):
            if category in categories_dictionary.keys():
                categories.append(category)
            else:
                print "ERROR: %s not found in mapping file, ommiting\n" % category

    # Command lists
    with open('qiime19.sh', 'w') as qiime19:
        # Mkdir Taxa summary <<<<<<<<<<
        qiime19.write("mkdir taxa_summary/\n")
        # Taxa summary commands
        summarize_taxa(biom_file, mapping_file, params_file, categories, qiime19)
        # Alpha div
        alpha_diversity(biom_file, mapping_file, params_file, tre_file, qiime19)
        # Alpha div colated
        compare_alpha_diversity(mapping_file, categories, qiime19)
        # Beta diversity
        beta_diversity(biom_file, mapping_file, tre_file, qiime19)
        # Comapre categories
        compare_beta(mapping_file, categories, qiime19)
        # Core microbiome
        qiime19.write("mkdir core_microbiome/\n")
        compute_core_microbiome(biom_file, mapping_file, categories, qiime19)
        # Group Significance
        group_significance(biom_file, mapping_file, categories, qiime19)

if __name__ == "__main__":
    main()
