#qiime_secondary_analysis  

The script expects both qiime 1.7 and 1.8 to be install. This script uses our convention of qiime17/qiime18 to load them into the path.
  
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
  --qiime17 QIIME17     The path to the Qiime 1.7 activate.sh or alias
  --qiime18 QIIME18     The path to the Qiime 1.8 activate.sh or alias
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
* `print_biom_table_summary.py`
* `mkdir taxa_summary/`
* `summarize_taxa_through_plots.py`
* `alpha_rarefaction.py`
* `compare_alpha_diversity.py`
* `beta_diversity_through_plots.py`
* `make_2d_plots.py`
* `compare_categories.py`
* `compute_core_microbiome.py`
* `source qiime17_path`
* `otu_category_significance.py`

## Output:
Output is written to the current working directory. 
A batch file for the Qiime 1.7/1.8 commands is created and run. 

