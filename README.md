# autoRNAseq
A bash script and some dependencies to automate RNA-seq analysis. A dummy file structure is provided as an example. <br>This script has been tested on MacOS Catalina, and should work on Unix-based systems.
<br>*Additional Requirements:* A sequence for adapater trimming (e.g. `TruSeq3-PE.fa`) and a reference cDNA sequence (e.g. `Homo_sapiens.GRCh38.cdna.all.fa`) are required.


## Usage
```
cd container #see the File Structure section of Inputs for more details
chmod 777 PROJECT1.sh #change script permissions to execute
. ./PROJECT1.sh #run script
```

## Inputs
The script depends on a particular file structure and `conda`. 

### File Structure
- container
	- `PROJECT1.sh`
	- PROJECT1
		- PRE
			- `sample1_identifier_L002_R1_001.fastq.gz`
			- `sample1_identifier_L002_R2_001.fastq.gz` 
			- `sample2_identifier_L002_R1_001.fastq.gz`
			- `sample2_identifier_L002_R2_001.fastq.gz`
		- POST
			- `sample3_identifier_L002_R1_001.fastq.gz`
			- `sample3_identifier_L002_R2_001.fastq.gz`
			- `sample4_identifier_L002_R1_001.fastq.gz`
			- `sample4_identifier_L002_R2_001.fastq.gz`
		- `analysis.R`
		- `studydesign.txt`
	- PROJECT2
	- `Homo_sapiens.GRCh38.cdna.all.fa`
	- `TruSeq3-PE.fa`

## R Analysis
The R analysis script (`analysis.R`) contains many properties that must be user-defined. `analysis.R` depends on `studydesign.txt`

### studydesign.txt
`studydesign.txt` contains a tab-seperated table of pertinent information for each sample. Formatting is as follows (based on the sample heirarchy provided in the **File Structure** section of **Inputs**):
sample|group|dir
---|---|---
PREsample1|PRE|sample1_identifier
PREsample2|PRE|sample2_identifier
POSTsample3|POST|sample3_identifier
POSTsample4|POST|sample4_identifier

### analysis.R
A barebones R script that outputs basic PCA analysis and a volcano plot of the RNA-seq data into `Rplots.pdf`. Additionally, outputs `LogFC.csv`. See **Output** for more details.

## Output
### logFC.csv
A comma-seperated-values file of every gene mapped by `Kallisto` and its corresponding P-val and log fold change from PRE to POST. Genes are sorted by `abs(logFC)`.
geneID|logFC|AveExpr|t|P.Value|adj.P.Value|B
---|---|---|---|---|---|---
gene1|1.11|2.22|3.33|.444|.555|6.66

### Rplots.pdf
Contains a volcano plot of the aforementioned genes as well as plots for sanity checking the data. 
