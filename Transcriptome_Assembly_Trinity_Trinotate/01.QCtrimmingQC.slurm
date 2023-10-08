#!/bin/bash
#SBATCH --partition=large
#SBATCH --job-name="01.eh.QCtrimmingQC"
#SBATCH --output="01.eh.QCtrimmingQC.o"
#SBATCH --error="01.eh.QCtrimmingQC.e"
#SBATCH --constraint="ethernet"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=20
#SBATCH --time=7-00:00:00
#SBATCH --account="generalGrp"
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=saundecj@wfu.edu
#SBATCH --mem=120GB
#################################### Now do your stuff ###################################

cd /deac/generalGrp/johnsonGrp/saundecj/projects/eiseniaHortensis.20180302

/deac/generalGrp/johnsonGrp/saundecj/bin/FastQC/fastqc \
	ehHead.S1.R1.fastq.gz ehHead.S1.R2.fastq.gz \
	ehHead.S2.R1.fastq.gz ehHead.S2.R2.fastq.gz \
	ehHead.S3.R1.fastq.gz ehHead.S3.R2.fastq.gz \
	ehMid.S4.R1.fastq.gz ehMid.S4.R2.fastq.gz \
	ehMid.S5.R1.fastq.gz ehMid.S5.R2.fastq.gz \
	ehMid.S6.R1.fastq.gz ehMid.S6.R2.fastq.gz \

java -classpath Trimmomatic-0.36/trimmomatic-0.36.jar org.usadellab.trimmomatic.TrimmomaticPE \
	-threads 12 -phred33 \
	ehHead.S1.R1.fastq.gz ehHead.S1.R2.fastq.gz \
	trimmedP.ehHead.S1.R1.fastq.gz trimmedUP.ehHead.S1.R1.fastq.gz \
	trimmedP.ehHead.S1.R2.fastq.gz trimmedUP.ehHead.S1.R2.fastq.gz \
	ILLUMINACLIP:Trimmomatic-0.36/adapters/customAdapters.fa:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50

java -classpath Trimmomatic-0.36/trimmomatic-0.36.jar org.usadellab.trimmomatic.TrimmomaticPE \
	-threads 20 -phred33 \
	ehHead.S2.R1.fastq.gz ehHead.S2.R2.fastq.gz \
	trimmedP.ehHead.S2.R1.fastq.gz trimmedUP.ehHead.S2.R1.fastq.gz \
	trimmedP.ehHead.S2.R2.fastq.gz trimmedUP.ehHead.S2.R2.fastq.gz \
	ILLUMINACLIP:Trimmomatic-0.36/adapters/customAdapters.fa:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50

java -classpath Trimmomatic-0.36/trimmomatic-0.36.jar org.usadellab.trimmomatic.TrimmomaticPE \
	-threads 20 -phred33 \
	ehHead.S3.R1.fastq.gz ehHead.S3.R2.fastq.gz \
	trimmedP.ehHead.S3.R1.fastq.gz trimmedUP.ehHead.S3.R1.fastq.gz \
	trimmedP.ehHead.S3.R2.fastq.gz trimmedUP.ehHead.S3.R2.fastq.gz \
	ILLUMINACLIP:Trimmomatic-0.36/adapters/customAdapters.fa:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50

java -classpath Trimmomatic-0.36/trimmomatic-0.36.jar org.usadellab.trimmomatic.TrimmomaticPE \
	-threads 20 -phred33 \
	ehMid.S4.R1.fastq.gz ehMid.S4.R2.fastq.gz \
	trimmedP.ehMid.S4.R1.fastq.gz trimmedUP.ehMid.S4.R1.fastq.gz \
	trimmedP.ehMid.S4.R2.fastq.gz trimmedUP.ehMid.S4.R2.fastq.gz \
	ILLUMINACLIP:Trimmomatic-0.36/adapters/customAdapters.fa:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50

java -classpath Trimmomatic-0.36/trimmomatic-0.36.jar org.usadellab.trimmomatic.TrimmomaticPE \
	-threads 20 -phred33 \
	ehMid.S5.R1.fastq.gz ehMid.S5.R2.fastq.gz \
	trimmedP.ehMid.S5.R1.fastq.gz trimmedUP.ehMid.S5.R1.fastq.gz \
	trimmedP.ehMid.S5.R2.fastq.gz trimmedUP.ehMid.S5.R2.fastq.gz \
	ILLUMINACLIP:Trimmomatic-0.36/adapters/customAdapters.fa:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50
	
java -classpath Trimmomatic-0.36/trimmomatic-0.36.jar org.usadellab.trimmomatic.TrimmomaticPE \
	-threads 20 -phred33 \
	ehMid.S6.R1.fastq.gz ehMid.S6.R2.fastq.gz \
	trimmedP.ehMid.S6.R1.fastq.gz trimmedUP.ehMid.S6.R1.fastq.gz \
	trimmedP.ehMid.S6.R2.fastq.gz trimmedUP.ehMid.S6.R2.fastq.gz \
	ILLUMINACLIP:Trimmomatic-0.36/adapters/customAdapters.fa:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50

/deac/generalGrp/johnsonGrp/saundecj/bin/FastQC/fastqc \
	trimmedP.ehHead.S1.R1.fastq.gz trimmedUP.ehHead.S1.R1.fastq.gz \
	trimmedP.ehHead.S1.R2.fastq.gz trimmedUP.ehHead.S1.R2.fastq.gz \
	trimmedP.ehHead.S2.R1.fastq.gz trimmedUP.ehHead.S2.R1.fastq.gz \
	trimmedP.ehHead.S2.R2.fastq.gz trimmedUP.ehHead.S2.R2.fastq.gz \
	trimmedP.ehHead.S3.R1.fastq.gz trimmedUP.ehHead.S3.R1.fastq.gz \
	trimmedP.ehHead.S3.R2.fastq.gz trimmedUP.ehHead.S3.R2.fastq.gz \
	trimmedP.ehMid.S4.R1.fastq.gz trimmedUP.ehMid.S4.R1.fastq.gz \
	trimmedP.ehMid.S4.R2.fastq.gz trimmedUP.ehMid.S4.R2.fastq.gz \
	trimmedP.ehMid.S5.R1.fastq.gz trimmedUP.ehMid.S5.R1.fastq.gz \
	trimmedP.ehMid.S5.R2.fastq.gz trimmedUP.ehMid.S5.R2.fastq.gz \
	trimmedP.ehMid.S6.R1.fastq.gz trimmedUP.ehMid.S6.R1.fastq.gz \
	trimmedP.ehMid.S6.R2.fastq.gz trimmedUP.ehMid.S6.R2.fastq.gz
exit 0
