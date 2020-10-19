# autoRNAseq
A bash script and some dependencies to automate RNA-seq analysis. A dummy file structure is provided as an example.

## Inputs
The script depends on a particular file structure and `conda`. 

### file structure
- RNA-seq
	- PROJECT1
		- PRE
			- sample1_identifier_L002_R1_001.fastq.gz 
			- sample1_identifier_L002_R2_001.fastq.gz 
			- sample2_identifier_L002_R1_001.fastq.gz 
			- sample2_identifier_L002_R2_001.fastq.gz 
		- POST
			- sample3_identifier_L002_R1_001.fastq.gz 
			- sample3_identifier_L002_R2_001.fastq.gz 
			- sample4_identifier_L002_R1_001.fastq.gz 
			- sample4_identifier_L002_R2_001.fastq.gz 
		- analysis.R
		- studydesign.txt
	- PROJECT2
	- Homo_sapiens.GRCh38.cdna.all.fa
	- TruSeq3-PE.fa

### studydesign.txt
*to be written*

### analysis.R
*to be written*

#### user defined 

### additional requirements
A sequence for adapater trimming (e.g. `TruSeq3-PE.fa`) and a reference cDNA sequence (e.g. `Homo_sapiens.GRCh38.cdna.all.fa`) are required.

## Output
*to be written*
