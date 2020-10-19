#!/bin/sh
conda activate bio
conda install --file ./requirements.txt

kallisto index -i Homo_sapiens.GRCh38.cdna.all.index Homo_sapiens.GRCh38.cdna.all.fa

PROJECT="PARPi"
cd ./${PROJECT}

mkdir ./processed
mv ./PARPi.R ./processed/PARPi.R
mv ./studydesign.txt ./processed/studydesign.txt

SAMPLE_NAME="LTL610"

SUBSAMPLE="610_Sample_1_S96"
cd ./${SAMPLE_NAME}
trimmomatic PE ./${SUBSAMPLE}_L002_R1_001.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001.fastq.gz \
	./${SUBSAMPLE}_L002_R1_001_trimmed.fastq.gz \
	./${SUBSAMPLE}_L002_R1_001_trimmed_unpaired.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001_trimmed.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001_trimmed_unpaired.fastq.gz \
	ILLUMINACLIP:../../TruSeq3-PE.fa:2:30:10:2:keepBothReads \
	LEADING:3 TRAILING:3 MINLEN:36
cd ..

mkdir ./processed/${SUBSAMPLE}
mkdir ./${SAMPLE_NAME}_fastqc/${SUBSAMPLE}
fastqc ./${SAMPLE_NAME}/*_trimmed.fastq.gz -o ./${SAMPLE_NAME}_fastqc/${SUBSAMPLE}

kallisto quant -i ../Homo_sapiens.GRCh38.cdna.all.index \
	-o ./processed/${SUBSAMPLE} \
	-t 4 \
	./${SAMPLE_NAME}/${SUBSAMPLE}_L002_R1_001_trimmed.fastq.gz \
	./${SAMPLE_NAME}/${SUBSAMPLE}_L002_R2_001_trimmed.fastq.gz \

SUBSAMPLE="610_Sample_2_S97"
cd ./${SAMPLE_NAME}
trimmomatic PE ./${SUBSAMPLE}_L002_R1_001.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001.fastq.gz \
	./${SUBSAMPLE}_L002_R1_001_trimmed.fastq.gz \
	./${SUBSAMPLE}_L002_R1_001_trimmed_unpaired.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001_trimmed.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001_trimmed_unpaired.fastq.gz \
	ILLUMINACLIP:../../TruSeq3-PE.fa:2:30:10:2:keepBothReads \
	LEADING:3 TRAILING:3 MINLEN:36
cd ..

mkdir ./processed/${SUBSAMPLE}
mkdir ./${SAMPLE_NAME}_fastqc/${SUBSAMPLE}
fastqc ./${SAMPLE_NAME}/*_trimmed.fastq.gz -o ./${SAMPLE_NAME}_fastqc/${SUBSAMPLE}

kallisto quant -i ../Homo_sapiens.GRCh38.cdna.all.index \
	-o ./processed/${SUBSAMPLE} \
	-t 4 \
	./${SAMPLE_NAME}/${SUBSAMPLE}_L002_R1_001_trimmed.fastq.gz \
	./${SAMPLE_NAME}/${SUBSAMPLE}_L002_R2_001_trimmed.fastq.gz \

SAMPLE_NAME="LTL610R"

SUBSAMPLE="610R_Sample_3_S98"
cd ./${SAMPLE_NAME}
trimmomatic PE ./${SUBSAMPLE}_L002_R1_001.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001.fastq.gz \
	./${SUBSAMPLE}_L002_R1_001_trimmed.fastq.gz \
	./${SUBSAMPLE}_L002_R1_001_trimmed_unpaired.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001_trimmed.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001_trimmed_unpaired.fastq.gz \
	ILLUMINACLIP:../../TruSeq3-PE.fa:2:30:10:2:keepBothReads \
	LEADING:3 TRAILING:3 MINLEN:36
cd ..

mkdir ./processed/${SUBSAMPLE}
mkdir ./${SAMPLE_NAME}_fastqc/${SUBSAMPLE}
fastqc ./${SAMPLE_NAME}/*_trimmed.fastq.gz -o ./${SAMPLE_NAME}_fastqc/${SUBSAMPLE}

kallisto quant -i ../Homo_sapiens.GRCh38.cdna.all.index \
	-o ./processed/${SUBSAMPLE} \
	-t 4 \
	./${SAMPLE_NAME}/${SUBSAMPLE}_L002_R1_001_trimmed.fastq.gz \
	./${SAMPLE_NAME}/${SUBSAMPLE}_L002_R2_001_trimmed.fastq.gz \

SUBSAMPLE="610R_Sample_4_S99"
cd ./${SAMPLE_NAME}
trimmomatic PE ./${SUBSAMPLE}_L002_R1_001.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001.fastq.gz \
	./${SUBSAMPLE}_L002_R1_001_trimmed.fastq.gz \
	./${SUBSAMPLE}_L002_R1_001_trimmed_unpaired.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001_trimmed.fastq.gz \
	./${SUBSAMPLE}_L002_R2_001_trimmed_unpaired.fastq.gz \
	ILLUMINACLIP:../../TruSeq3-PE.fa:2:30:10:2:keepBothReads \
	LEADING:3 TRAILING:3 MINLEN:36
cd ..

mkdir ./processed/${SUBSAMPLE}
mkdir ./${SAMPLE_NAME}_fastqc/${SUBSAMPLE}
fastqc ./${SAMPLE_NAME}/*_trimmed.fastq.gz -o ./${SAMPLE_NAME}_fastqc/${SUBSAMPLE}

kallisto quant -i ../Homo_sapiens.GRCh38.cdna.all.index \
	-o ./processed/${SUBSAMPLE} \
	-t 4 \
	./${SAMPLE_NAME}/${SUBSAMPLE}_L002_R1_001_trimmed.fastq.gz \
	./${SAMPLE_NAME}/${SUBSAMPLE}_L002_R2_001_trimmed.fastq.gz \

cd ./processed
Rscript ./analyze.R
mv ./processed/PARPi.R ./PARPi.R
mv ./processed/studydesign.txt ./studydesign.txt