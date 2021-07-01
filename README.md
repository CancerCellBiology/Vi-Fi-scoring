# Vi-Fi-scoring

Requirements:

Python >= 3.7.4

Dependencies: 

For Vi-Fi scoring:
1. pandas
2. numpy
3. scipy
4. itertools

For finding drugs in ChEMBL and DSigDB:
1. pandas
2. numpy
3. chembl_webresource_client
4. urllib

Vi-Fi plot:
1. pandas
2. matplotlib


Before start:
1) Create Vi-Fi scoring folder and add proper address to the codes
2) Unpack Datasets archive into Vi-Fi folder

Find drugs by genes: takes a list of genes to find drugs in ChEMBL and DSigDB database that target query genes/proteins. For ChEMBL search provide Uniprot ids, for DSigDB HUGO gene names. Returns tables of drugs with their respective targets, which are stored in th Results folder. Use Test_genes as a test file and as example of input data.
ChEMBL_Results- refers to a search result within ChEMBL for drugs with known mechanisms.
DSigDB_kinases- refers to a search result within D1 (approved drugs) and D2 (kinase inhibitors) which yields drugs based on their confirmed direct targets.
DSigDB_mining- refers to a search result within D4 (TTD and CTD) which yields drug based on the targets identified by text mining.
For more information refer to:
https://www.ebi.ac.uk/chembl/
http://dsigdb.tanlab.org/DSigDBv1.0/

Vi-Fi scoring: takes a list of drugs and returns a table with fibrotic and viral scores for drugs present in L1000 database of drug-induces signatures for A549 cells. Also saves fibrotic and viral signatures of each drug to Vi-Fi Signatures folder. Use Test_drugs as a test file and as example of input data or the results from "Find drugs by genes."
For the information on L1000 refer to:
https://maayanlab.cloud/L1000FWD/
Current database version can be downloaded from:
https://maayanlab.cloud/L1000FWD/download_page (CD signatures (up/down gene sets) in the full space in json format)

Vi-Fi combination scores: takes a list of drugs and returns a table with fibrotic and viral scores for each possible drug combination. This algorithm uses saved fibrotic and viral signatures from Vi-Fi Signatures folder to speed the process. Use Test_drugs as a test file and as example of input data. Currently this code works only if Vi-Fi Signatures folder was created and individual scores were calculated for each drug with Vi-Fi scoring.

Vi-Fi plot: code used for creation of 2D Vi-Fi plot. Use "Combined ViFi scoring.xlsx" as an input example.

All code results are published at Scientific Reports: https://www.nature.com/articles/s41598-021-90701-w
