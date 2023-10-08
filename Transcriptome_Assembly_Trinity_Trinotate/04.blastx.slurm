#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name="04.blastx"
#SBATCH --output="04.blastx.o"
#SBATCH --error="04.blastx.e"
#SBATCH --constraint="ethernet"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=20
#SBATCH --time=4-00:00:00
#SBATCH --account="generalGrp"
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=saundecj@wfu.edu
#SBATCH --mem=128GB

###Now do your stuff

cd /deac/generalGrp/johnsonGrp/saundecj/projects/eiseniaHortensis.20180302

# ran before running this script
# makeblastdb -in ../../omes/proteome/proteome.Cele-Dmel-Hrob-Hsap.faa  -dbtype prot -parse_seqids -out proteome.Cele-Dmel-Hrob-Hsap.blastdb

blastx -num_threads 20 -query outs/cd-hit-est.output -db proteome.Cele-Dmel-Hrob-Hsap.blastdb -max_target_seqs 1 -outfmt 6 -out blastx.fasta

exit 0
#blastx -num_threads 20 -query outs/cd-hit-est.output -db proteome.Cele-Dmel-Hrob-Hsap.blastdb -max_target_seqs 1 -outfmt 6
