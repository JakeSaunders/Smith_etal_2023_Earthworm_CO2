#!/bin/bash
#SBATCH --partition=large
#SBATCH --job-name="02.eh.trinity"
#SBATCH --output="02.eh.trinity.o"
#SBATCH --error="02.eh.trinity.e"
#SBATCH --constraint="ethernet"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=20
#SBATCH --time=7-00:00:00
#SBATCH --account="generalGrp"
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=saundecj@wfu.edu
#SBATCH --mem=220G
##################################### 	Now do your stuff

/deac/generalGrp/johnsonGrp/saundecj/projects/eiseniaHortensis.20180302

module load java/1.8.0_101
module load trinity/2.4.0-gcc

Trinity --seqType fq \
	--left trimmedP.ehHead.S1.R1.fastq.gz,trimmedP.ehHead.S2.R1.fastq.gz,trimmedP.ehHead.S3.R1.fastq.gz,trimmedP.ehMid.S4.R1.fastq.gz,trimmedP.ehMid.S5.R1.fastq.gz,trimmedP.ehMid.S6.R1.fastq.gz \
	--right trimmedP.ehHead.S1.R2.fastq.gz,trimmedP.ehHead.S2.R2.fastq.gz,trimmedP.ehHead.S3.R2.fastq.gz,trimmedP.ehMid.S4.R2.fastq.gz,trimmedP.ehMid.S5.R2.fastq.gz,trimmedP.ehMid.S6.R2.fastq.gz \
	--CPU 20 --max_memory 220G --output trinity.eiseniaHortensis

module unload trinity/2.4.0-gcc
module unload java/1.8.0_101
exit 0
