#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name="04.trinotateBlast"
#SBATCH --output="04.trinotateBlast.o"
#SBATCH --error="04.trinotateBlast.e"
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

/deac/generalGrp/johnsonGrp/saundecj/bin/ncbi-blast-2.6.0+/bin/blastx \
	-query Trinity.fasta -db uniprot_sprot.pep \
	-num_threads 20 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6

/deac/generalGrp/johnsonGrp/saundecj/bin/ncbi-blast-2.6.0+/bin/blastp \
	-query longest_orfs.pep -db uniprot_sprot.pep \
	-num_threads 20 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6

exit 0

