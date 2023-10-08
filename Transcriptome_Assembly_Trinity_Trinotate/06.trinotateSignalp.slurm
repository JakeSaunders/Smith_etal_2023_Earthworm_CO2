#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name="06.trinotateSignalp"
#SBATCH --output="06.trinotateSignalp.o"
#SBATCH --error="06.trinotateSignalp.e"
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

../signalp-4.1/signalp -f short -n signalp.out longest_orfs.pep

exit 0
