#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name="07.trinotateTMHMM"
#SBATCH --output="07.trinotateTMHMM.o"
#SBATCH --error="07.trinotateTMHMM.e"
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

../tmhmm-2.0c/bin/tmhmm --short < longest_orfs.pep > tmhmm.out

exit 0
