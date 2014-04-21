#qiime_secondary_analysis  
  
##A meta script to take care of the Qiime secondary analysis  
Usage: secondary_analysis.py [options]

This script runs through our standard QIIME secondary analysis pipeline. The required input files are the biom, map, tre, and params.  

The pipeline runs the following steps:
Valadate the parameters file  
Check if Qiime is loaded in the path  
Convert the mapping file to a dictionary  
Summarize biom table  
Summarize taxa through plots  
Compute/Compare Alpha Diversity    
Compute/Compare Beta Diversity  
Compute core microbiome 
Compute OTU category significance  


###Arguments:
**h, --help**
>show this help message and exit  

**-b BIOM, --biom=BIOM**  
>The biom file  

**-m MAP, --map=MAP**
>The mapping file  

**-p PARAMS, --params=PARAMS**
>The parameters file  

**-t TRE, --tre=TRE**
>The tre file  

**-c CATEGORIES, --categories=CATEGORIES**
>The metadata categories to compute. Must be colon seperated  




