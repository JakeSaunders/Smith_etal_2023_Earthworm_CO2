#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name="08.trinotateRNAMMER"
#SBATCH --output="08.trinotateRNAMMER.o"
#SBATCH --error="08.trinotateRNAMMER.e"
#SBATCH --constraint="ethernet"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=20
#SBATCH --time=4-00:00:00
#SBATCH --account="generalGrp"
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=saundecj@wfu.edu
#SBATCH --mem=120GB

###Now do your stuff

cd /deac/generalGrp/johnsonGrp/saundecj/bin/Trinotate-Trinotate-v3.1.1

./util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer ../RNAMMER/rnammer

exit 0
