#!/usr/bin/env python

__author__ = "Gene Blanchard"
__email__ = "me@geneblanchard.com"

# Imports
import argparse
#from optparse import OptionParser
import sys
import subprocess
import os
import time
import datetime
import csv
from itertools import izip

'''
Python secondary analysis script

This script takes the basic primary analysis files (biom, map, tre, params) and runs the SOP commands

'''

# Validate that the parameters file has all the values needed
def valid_params(params_file):
	required_values =['summarize_taxa:level','plot_taxa_summary:labels','alpha_diversity:metrics','multiple_rarefactions:min','multiple_rarefactions:max','multiple_rarefactions:step','beta_diversity_through_plots:seqs_per_sample']
	found_values = []
	with open(params_file, 'r') as params:
		for line in params:
			value = line.split(' ')[0]
			found_values.append(value)
	for required in required_values:
		if required not in found_values:
			print required+" not found in parameters file.\n"
			print "\nI.E.:\nsummarize_taxa:level 2,3,4,5,6,7\nplot_taxa_summary:labels Phylum,Class,Order,Family,Genus,Species\nalpha_diversity:metrics shannon,simpson,PD_whole_tree,chao1,observed_species\nmultiple_rarefactions:min 100\nmultiple_rarefactions:max 18000\nmultiple_rarefactions:step 500\nbeta_diversity_through_plots:seqs_per_sample 18000\n"
			sys.exit()

# Function to make sure qiime is loaded in the PATH
def check_qiime(log):
	try:
		subprocess.Popen("print_qiime_config.py", stdout=log) 
		
	except OSError:
		print "Qiime is not loaded into your path"
		sys.exit()

# Convert the map file to a dictionary 
def map_to_dictionary(mapfile):
	dictionary =  dict()
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
def create_biom_summary(biom_file, log):
	command = "print_biom_table_summary.py -i %s -o table_summary.txt" % biom_file
	log.write(command+'\n')
	biom_proc = subprocess.Popen(command.split(' '))
	biom_proc.wait()

# Get the metadata columns 
def read_mapfile_metadata(map_file):
	with open(map_file, 'r') as mapping_file:
		categories = mapping_file.readline()
	metadata_categories = categories.rstrip('\n').split('\t')
	try:
		metadata_categories.remove('#SampleID')
		metadata_categories.remove('Description')
		metadata_categories.remove('BarcodeSequence')
		metadata_categories.remove('LinkerPrimerSequence')
	except ValueError:
		pass
	return metadata_categories

def qiime17_otu_category_sig(biom, mapfile, categories, log):
	tests = ['g_test', 'ANOVA']
	temp_stamp = str(int(time.time()))
	with open('otu_category_significance_batch_'+temp_stamp, 'w') as batch:
		batch.write('#!/bin/bash\n. sh /media/nfs_opt/qiime17/activate.sh\n')
		for category in categories:
			for test in tests:
				otu_cat_sig_command = "otu_category_significance.py -i %s -m %s -c %s -s %s -o taxa_summary/taxa_%s/%s.txt&" % (biom, mapfile, category, test, category, test)
				batch.write(otu_cat_sig_command+'\n')
				log.write(otu_cat_sig_command+'\n')
	os.chmod('otu_category_significance_batch_'+temp_stamp, 0755)
	otu_cat_sig_proc = subprocess.Popen(['./otu_category_significance_batch_'+temp_stamp])
	otu_cat_sig_proc.wait()

def summarize_taxa(biom, mapfile, params, categories, log):
	# Make Taxa Summary Directory
	os.mkdir('taxa_summary')
	basic_taxa_command = "summarize_taxa_through_plots.py -s -i %s -m %s -p %s -o taxa_summary/taxa_individual" % (biom, mapfile, params)
	log.write(basic_taxa_command+'\n')
	basic_taxa_proc = subprocess.Popen(basic_taxa_command.split(' '))
	basic_taxa_proc.wait()
	for category in categories:
		category_taxa_command = "summarize_taxa_through_plots.py -i %s -m %s -p %s -o taxa_summary/taxa_%s -c %s -s" % (biom, mapfile, params, category, category)
		log.write(category_taxa_command+'\n')
		category_taxa_proc = subprocess.Popen(category_taxa_command.split(' '))
		category_taxa_proc.wait()

def alpha_diversity(biom, mapfile, params, tre, log):
	alpha_command = "alpha_rarefaction.py -i %s -m %s -p %s -t %s -a -O 24 -o alpha_diversity" % (biom, mapfile, params, tre)
	log.write(alpha_command+'\n')
	alpha_proc = subprocess.Popen(alpha_command.split(' '))
	alpha_proc.wait()

def compare_alpha_diversity(cwd, mapfile, categories, log):
	alpha_diversity_path = 'alpha_diversity/alpha_div_collated/'
	metrics = ['chao1', 'observed_species', 'PD_whole_tree', 'shannon', 'simpson']
	temp_stamp = str(time.time())[:10]
	batchfile_name = alpha_diversity_path+'batch'+temp_stamp
	for metric in metrics:
		for category in categories:
			metric_file = alpha_diversity_path+metric
			output =  alpha_diversity_path+metric+'_'+category
			compare_alpha_diversity_command = "compare_alpha_diversity.py -m %s -n 9999 -c %s -i %s.txt -o %s &\n" % (mapfile, category, metric_file, output)
			with open(batchfile_name, 'a') as batchfile:
				batchfile.write(compare_alpha_diversity_command)
				log.write(compare_alpha_diversity_command)
	compare_alpha_diversity_proc = subprocess.Popen(['sh',batchfile_name])
	compare_alpha_diversity_proc.wait()

def beta_diversity(biom, mapfile, tre, log):
	beta_diversity_command = "beta_diversity_through_plots.py -i %s -m %s -t %s -o beta_diversity" % (biom, mapfile, tre)
	log.write(beta_diversity_command+'\n')
	beta_diversity_proc =  subprocess.Popen(beta_diversity_command.split(' '))
	beta_diversity_proc.wait()
	# Make the 2d plots
	make_2d_unweighted_command = "make_2d_plots.py -i beta_diversity/unweighted_unifrac_pc.txt -m %s -o beta_diversity/2d_unweighted_unifrac_plots" % mapfile
	make_2d_weighted_command = "make_2d_plots.py -i beta_diversity/weighted_unifrac_pc.txt -m %s -o beta_diversity/2d_weighted_unifrac_plots" % mapfile
	log.write(make_2d_unweighted_command+'\n')
	log.write(make_2d_weighted_command+'\n')
	make_2d_unweighted_proc = subprocess.Popen(make_2d_unweighted_command.split(' '), stderr=log)
	make_2d_unweighted_proc.wait()
	make_2d_weighted_proc = subprocess.Popen(make_2d_weighted_command.split(' '), stderr=log)
	make_2d_weighted_proc.wait()

def compare_beta(mapfile, categories, log):
	for category in categories:
		compare_beta_unweighted_command = "compare_categories.py --method anosim -i beta_diversity/unweighted_unifrac_dm.txt -m %s -c %s -n 9999 -o beta_diversity/ANOSIM_%s_unweighted" % (mapfile, category, category)
		compare_beta_weighted_command = "compare_categories.py --method anosim -i beta_diversity/weighted_unifrac_dm.txt -m %s -c %s -n 9999 -o beta_diversity/ANOSIM_%s_weighted" % (mapfile, category, category)
		log.write(compare_beta_unweighted_command+'\n')
		log.write(compare_beta_weighted_command+'\n')
		compare_beta_unweighted_proc = subprocess.Popen(compare_beta_unweighted_command.split(' '))
		compare_beta_unweighted_proc.wait()
		compare_beta_weighted_proc = subprocess.Popen(compare_beta_weighted_command.split(' '))
		compare_beta_weighted_proc.wait()

def compute_core_microbiome(biom, mapfile, categories, log):
	os.mkdir('core_microbiome')
	state_dictionary = map_to_dictionary(mapfile)
	for category in categories:
		for state in state_dictionary[category]:
			core_microbiome_command = 'compute_core_microbiome.py -i %s --mapping_fp %s --num_fraction_for_core_steps 6 -o core_microbiome/core_microbiome_%s_%s --valid_states %s:%s' % (biom, mapfile, category, state, category, state)
			log.write(core_microbiome_command+'\n')
			core_microbiome_proc = subprocess.Popen(core_microbiome_command.split(' '))
			core_microbiome_proc.wait()
			

def get_time():
	now = datetime.datetime.now().ctime()
	return str(now)+'\n'

def hms_string(sec_elapsed):
	h = int(sec_elapsed / (60 * 60))
	m = int((sec_elapsed % (60 * 60)) / 60)
	s = sec_elapsed % 60.
	return "{}:{:>02}:{:>05.2f}".format(h, m, s)

def main():
	# Set current working directory
	cwd = os.getcwd()

	# Create the argument parser
	parser = argparse.ArgumentParser(description="Need description")

	#biom -b --biom
	parser.add_argument("-b", "--biom", dest="biom", required=True, help="The biom file")
	#map -m --map
	parser.add_argument("-m", "--map", dest="map", required=True, help="The mapping file")
	#parameters -p --params
	parser.add_argument("-p", "--params", dest="params", required=True, help="The parameters file")
	#tre -t --tre
	parser.add_argument("-t", "--tre", dest="tre", required=True, help="The tre file")
	#categories -c --categories
	parser.add_argument("-c", "--categories", dest="categories", help="The metadata categories to compute. Must be colon seperated")

	# Parse the arguments
	args = parser.parse_args()

	# Assign variables
	biom_file = args.biom
	mapping_file = args.map
	params_file = args.params
	tre_file =  args.tre

	# <><><><><><><><><><><><><
	categories_dictionary = map_to_dictionary(mapping_file)

	if not options.categories == None:
		categories = []
		input_categories = options.categories.split(':')
		expected_categories = read_mapfile_metadata(mapping_file)
		for category in input_categories:
			if category in expected_categories:
				categories.append(category)
			else:
				print "ERROR: %s not found in mapping file, ommiting\n" % category
		
		#CATEGORIES = True
	else:
		categories = read_mapfile_metadata(mapping_file)

	time_stamp = str(int(time.time()))
	logfile = 'secondary_batch_%s.log' % time_stamp

	with open(logfile, 'w') as log:

		# Print Categories
		log.write("Categories:"+str(categories)+'\n')
		
		# Validate Parameters
		start = time.time()
		log.write("Validate Parameters\t"+get_time())
		valid_params(params_file)
		elapsed = time.time() - start
		hms_elapsed = "Elapsed {}\n\n".format(hms_string(elapsed))
		log.write(hms_elapsed)
		log.flush()
		
		# Check qiime in path
		start = time.time()
		log.write("Check if Qiime is in $PATH\t"+get_time())
		check_qiime(log)
		elapsed = time.time() - start
		hms_elapsed = "Elapsed {}\n\n".format(hms_string(elapsed))
		log.write(hms_elapsed)
		log.flush()
		
		# Summarize BIOM table
		start = time.time()
		log.write("BIOM Table Summary\t"+get_time())
		create_biom_summary(biom_file, log)
		elapsed = time.time() - start
		hms_elapsed = "Elapsed {}\n\n".format(hms_string(elapsed))
		log.write(hms_elapsed)
		log.flush()
		
		# Summarize taxa through plots
		start = time.time()
		log.write("Summarize Taxa\t"+get_time())
		summarize_taxa(biom_file, mapping_file, params_file, categories, log)
		elapsed = time.time() - start
		hms_elapsed = "Elapsed {}\n\n".format(hms_string(elapsed))
		log.write(hms_elapsed)
		log.flush()
		
		# Compute Alpha diversity
		start = time.time()
		log.write("Alpha Diversity\t"+get_time())
		alpha_diversity(biom_file, mapping_file, params_file, tre_file, log)
		elapsed = time.time() - start
		hms_elapsed = "Elapsed {}\n\n".format(hms_string(elapsed))
		log.write(hms_elapsed)
		log.flush()
		
		# Compare Alpha diversity
		start = time.time()
		log.write("Compare Alpha Diversity\t"+get_time())
		compare_alpha_diversity(cwd, mapping_file, categories, log)
		elapsed = time.time() - start
		hms_elapsed = "Elapsed {}\n\n".format(hms_string(elapsed))
		log.write(hms_elapsed)
		log.flush()
		
		# Compute Beta diversity
		start = time.time()
		log.write("Beta Diversity\t"+get_time())
		beta_diversity(biom_file, mapping_file, tre_file, log)
		elapsed = time.time() - start
		hms_elapsed = "Elapsed {}\n\n".format(hms_string(elapsed))
		log.write(hms_elapsed)
		log.flush()
		
		# Compare Beta Diversity
		start = time.time()
		log.write("Compare Beta Diversity\t"+get_time())
		compare_beta(mapping_file, categories, log)
		elapsed = time.time() - start
		hms_elapsed = "Elapsed {}\n\n".format(hms_string(elapsed))
		log.write(hms_elapsed)
		log.flush()
		
		# Compute core microbiome
		start = time.time()
		log.write("Core Microbiome\t"+get_time())
		compute_core_microbiome(biom_file, mapping_file, categories, log)
		elapsed = time.time() - start
		hms_elapsed = "Elapsed {}\n\n".format(hms_string(elapsed))
		log.write(hms_elapsed)
		log.flush()
		
		# OTU category sig
		start = time.time()
		log.write("OTU Category Significance\t"+get_time())
		qiime17_otu_category_sig(biom_file, mapping_file, categories, log)
		elapsed = time.time() - start
		hms_elapsed = "Elapsed {}\n\n".format(hms_string(elapsed))
		log.write(hms_elapsed)
		log.flush()




if __name__ == "__main__":
	main()

