# MIQ-calc-from-OTU-tables
MIQ score calculator based on the original miqScore16SPublic tool from Zymo that accepts OTU tables generated from different clustering approaches

This calculator only works with ZymoBIOMICS Microbial Community Standard (Cat. No D6300). You can modify the script's expected data for any other mock type.

Formatting the input File:
Use the Example-input-format.tsv file in the directory. Basically, it is the extracted .qza file from Qiime2 and then the biom was converted to tsv.
1. Make a copy of the file and rename it. 
2. Copy and paste your count data from the Qiime2 extracted tables or any other. You can add as many columns (samples) as you wish.
3. ***IMPORTANT***. The column #OTU ID must be identical to the one in the script or Example-input-format.tsv! if you submit a different OTU table, make sure to copy-paste the same exact taxonomic ranks. The accuracy or the level (genus/phylum) of your local taxonomic assignments is not assessed here, rather the composition is what matters.
In the case where in your OTU table you have more than 8 (for 16S) or 10 (shotgun) microorganisms, you should consider omitting the falsely identified OTUs/bacteria. Changing the expected data within the script to expect fewer OTUs (<8 for 16S or <10 for shotgun) yields varying and unverified results.
Row order does not matter. You don't need to delete the final two rows (yeasts) if you submit only 16S rRNA derived counts. Leave them blank. You can delete the columns from the example.

**Usage:**
1. Activate a python environment with installed pandas package
2. python script.py <data_type> <sample_data_file> <dna_extraction_kit> where:
   
**<data_type>** is either "shotgun" or "amplicon"

**<sample_data_file>** - directory of the provided pre-formatted OTU file

**<dna_extraction_kit>** - any string name of choice

**Output:**

It is printed in the terminal. You can redirect the output to a separate file. Each column is considered a separate sample and the name is printed at the end.

**Example:** 

$ python calc-RMSE-final.py amplicon Example-input-format.tsv Example-name

Should print:

Samples name: Example-name

MIQ score for EZNA-Com: 72, RMSE: 27.90

MIQ score for EurX-Com: 70, RMSE: 29.63

MIQ score for FastSpin-Com: 80, RMSE: 19.45

MIQ score for In-house-LPA-Com: 79, RMSE: 20.97

MIQ score for Zymo-Com: 62, RMSE: 37.94

MIQ score for Zymo-Com-55C: 73, RMSE: 26.56

MIQ score for Zymo-Com-58C: 69, RMSE: 30.85

MIQ score for Zymo-Com-62C: 61, RMSE: 38.40

or 
$ python calc-RMSE-final.py amplicon Example-input-format.tsv Example-name >> example.tsv

should create an example.tsv file
