#Qiime Secondary Analysis  

The script expects Qiime 1.7, 1.8, & 1.9 to be installed. 
I have an alias to each Qiime install's `activate.sh`, these are qiime17, qiime18, and qiime19. 
You can pass the path to the activate file with the 
  
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
  --commands            Just write commands, don't run                        
  --qiime17 QIIME17     The path to the Qiime 1.7 activate.sh or alias
  --qiime18 QIIME18     The path to the Qiime 1.8 activate.sh or alias
  --qiime19 QIIME19     The path to the Qiime 1.9 activate.sh or alias
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
* [`otu_category_significance.py`](http://qiime.org/1.7.0/scripts/otu_category_significance.html)

## Output:
Output is written to the current working directory. 
A batch file for the Qiime 1.7/1.8 commands is created and run. 

