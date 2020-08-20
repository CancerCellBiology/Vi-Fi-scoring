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
2) Unpack Datasets archive into Datasets folder within Vi-Fi scoring folder
3) In Vi-Fi scoring create Vi-Fi Signatures folder to store drug signatures

ChEMBL find drugs: takes a list of genes (as UNIPROT ids) to find drugs in ChEMBL database that target query genes/proteins. Returns a table of drugs with their respective targets. Use Test_genes as a test file and as example of input data.

DSigDB find drugs: takes a list of genes (as HUGO names) to find drugs in DsigDB database that target query genes/proteins. Returns two tables: drugs identified in kinase inhibitors databases and drugs identified in text mining and computational databases. Use Test_genes as a test file and as example of input data.

Vi-Fi scoring: takes a list of drugs and returns a table with fibrotic and viral scores for drugs present in L1000 database of drug-induces signatures. Also saves fibrotic and viral signatures of each drug to Vi-Fi Signatures folder. Use Test_drugs as a test file and as example of input data.

Vi-Fi combination scores: takes a list of drugs and returns a table with fibrotic and viral scores for each possible drug combination. This algorithm uses saved fibrotic and viral signatures from Vi-Fi Signatures folder to speed the process. Use Test_drugs as a test file and as example of input data. Currently this code works only if Vi-Fi Signatures folder was created and individual scores were calculated for each drug with Vi-Fi scoring.

Vi-Fi plot: code used for creation of 2D Vi-Fi plot.
