#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name="05.trinotateHMMER"
#SBATCH --output="05.trinotateHMMER.o"
#SBATCH --error="05.trinotateHMMER.e"
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

../hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan --cpu 20 --domtblout TrinotatePFAM.out Pfam-A.hmm longest_orfs.pep > pfam.log

exit 0
