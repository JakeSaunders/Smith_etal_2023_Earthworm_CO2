#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name="09.EstAbunBWA"
#SBATCH --output="09.EstAbunBWA.o"
#SBATCH --error="09.EstAbunBWA.e"
#SBATCH --constraint="ethernet"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=20
#SBATCH --time=0-90:00:00
#SBATCH --account="generalGrp"
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=saundecj@wfu.edu
#SBATCH --mem=128GB

###Now do your stuff
echo "change dir"
cd /deac/generalGrp/johnsonGrp/saundecj/projects/eiseniaHortensis.20180302/

# ran in command line
# /deac/generalGrp/johnsonGrp/saundecj/bin/chado_test/chado/bin/gmod_fasta2gff3.pl --fasta_dir  trinotate.hits.fasta --gfffilename  trinotate.hits.gff3 --type CDS --nosequence

echo "load bwa"
module load bwa/0.7.12-icc

echo "start head sample 1"
bwa mem trinotate.hits.fasta trimmedP.ehHead.S1.R1.fastq.gz	trimmedP.ehHead.S1.R2.fastq.gz > eh.head.S1.sam

echo "start head sample 2"
bwa mem trinotate.hits.fasta trimmedP.ehHead.S2.R1.fastq.gz	trimmedP.ehHead.S2.R2.fastq.gz > eh.head.S2.sam

echo "start head sample 3"
bwa mem trinotate.hits.fasta trimmedP.ehHead.S3.R1.fastq.gz	trimmedP.ehHead.S3.R2.fastq.gz > eh.head.S3.sam

echo "start mid sample 4"
bwa mem trinotate.hits.fasta trimmedP.ehMid.S4.R1.fastq.gz trimmedP.ehMid.S4.R2.fastq.gz > eh.mid.S4.sam 

echo "start mid sample 5"
bwa mem trinotate.hits.fasta trimmedP.ehMid.S5.R1.fastq.gz trimmedP.ehMid.S5.R2.fastq.gz > eh.mid.S5.sam 

echo "start mid sample 6"
bwa mem trinotate.hits.fasta trimmedP.ehMid.S6.R1.fastq.gz trimmedP.ehMid.S6.R2.fastq.gz > eh.mid.S6.sam 

echo "Finished"
module unload bwa/0.7.12-icc
exit 0