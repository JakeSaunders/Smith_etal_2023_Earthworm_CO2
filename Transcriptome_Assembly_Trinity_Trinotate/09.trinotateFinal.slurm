#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name="09.trinotateFinal"
#SBATCH --output="09.trinotateFinal.o"
#SBATCH --error="09.trinotateFinal.e"
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
alias Trinotate='/deac/generalGrp/johnsonGrp/saundecj/bin/Trinotate-Trinotate-v3.1.1/Trinotate'

#1. Load transcripts and coding regions
## generate Gene/Transcript relationships using trinity already done
#./util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta >  Trinity.fasta.gene_trans_map

Trinotate Trinotate.sqlite init \
	--gene_trans_map Trinity.fasta.gene_trans_map \
	--transcript_fasta Trinity.fasta \
	--transdecoder_pep longest_orfs.pep

#2. Loading BLAST homologies
## load protein hits
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
## load transcript hits
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6

#3. Load Pfam domain entries
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out

#4. Load transmembrane domains
Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out

#5. Load signal peptide predictions
Trinotate Trinotate.sqlite LOAD_signalp signalp.out

# Trinotate: Output an Annotation Report
Trinotate Trinotate.sqlite report --incl_pep --incl_trans > trinotate_annotation_report.xls

exit 0
