#Qiime Secondary Analysis  

This script creates a batch script for our standard Qiime analysis. It is a *very*
thorough script that runs nearly all possible test. You may need to fine tune
some parameters such as the `compare_beta` iterations 
## Help:
```
usage: secondary_analysis.py [-h] -b BIOM -m MAP -p PARAMS -t TRE
                             [-c CATEGORIES] [--qiime17 QIIME17]
                             [--qiime18 QIIME18]

This script runs through a standard QIIME secondary analysis pipeline. The
required input files are the biom, map, tre, and params.

optional arguments:
  -h, --help            show this help message and exit
  -b BIOM, --biom BIOM  The biom file
  -m MAP, --map MAP     The mapping file
  -p PARAMS, --params PARAMS
                        The parameters file
  -t TRE, --tre TRE     The tre file
  -c CATEGORIES, --categories CATEGORIES
                        The metadata categories to compute. Must be colon
                        seperated
```
## Workflow
The pipeline runs the following steps:  
* Valadate the parameters file  
* Check if Qiime is loaded in the path  
* Convert the mapping file to a dictionary  
* Summarize biom table  
* Summarize taxa through plots  
* Compute/Compare Alpha Diversity  
* Compute/Compare Beta Diversity  
* Compute core microbiome  
* Compute OTU category significance  (Qiime 1.7)

The Qiime scripts run are:
* `source qiime18_path`
* [`print_biom_table_summary.py`](http://biom-format.org/documentation/summarizing_biom_tables.html)
* [`summarize_taxa_through_plots.py`](http://qiime.org/scripts/summarize_taxa_through_plots.html)
* [`alpha_rarefaction.py`](http://qiime.org/scripts/alpha_rarefaction.html)
* [`compare_alpha_diversity.py`](http://qiime.org/scripts/compare_alpha_diversity.html)
* [`beta_diversity_through_plots.py`](http://qiime.org/scripts/beta_diversity_through_plots.html)
* [`make_2d_plots.py`](http://qiime.org/scripts/make_2d_plots.html)
* [`compare_categories.py`](http://qiime.org/scripts/compare_categories.html)
* [`compute_core_microbiome.py`](http://qiime.org/scripts/compute_core_microbiome.html)
* `source qiime17_path`
* [`group_significance.py `](http://qiime.org/scripts/group_significance.html)

## Output:
Output is written to the current working directory.
A batch file for the Qiime 1.9 commands is created
