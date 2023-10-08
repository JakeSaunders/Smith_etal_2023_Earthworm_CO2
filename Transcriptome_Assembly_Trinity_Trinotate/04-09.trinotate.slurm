#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name="04.trinotateBlast"
#SBATCH --output="04.trinotateBlast.o"
#SBATCH --error="04.trinotateBlast.e"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=20
#SBATCH --time=4-00:00:00
#SBATCH --account="generalGrp"
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=saundecj@wfu.edu
#SBATCH --mem=120GB

###Now do your stuff

####	Set working dir and alias
cd /deac/generalGrp/johnsonGrp/saundecj/bin/Trinotate-Trinotate-v3.1.1
alias Trinotate='/deac/generalGrp/johnsonGrp/saundecj/bin/Trinotate-Trinotate-v3.1.1/Trinotate'

#### 	Set up required files
TransDecoder.LongOrfs -t Trinity.fasta
cp Trinity.fasta.transdecoder_dir/longest_orfs.pep .

####	BLAST 
/deac/generalGrp/johnsonGrp/saundecj/bin/ncbi-blast-2.6.0+/bin/blastx \
	-query Trinity.fasta -db uniprot_sprot.pep \
	-num_threads 20 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6

/deac/generalGrp/johnsonGrp/saundecj/bin/ncbi-blast-2.6.0+/bin/blastp \
	-query longest_orfs.pep -db uniprot_sprot.pep \
	-num_threads 20 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
	
####	HMMER
../hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan \
	--cpu 12 --domtblout TrinotatePFAM.out Pfam-A.hmm longest_orfs.pep > pfam.log 

####	SignalIP
../signalp-4.1/signalp -f short -n signalp.out longest_orfs.pep

####	TMHMM
../tmhmm-2.0c/bin/tmhmm --short < longest_orfs.pep > tmhmm.out

####	RNAMMER
./util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer ../RNAMMER/rnammer

####	Trinotate
../trinityrnaseq-Trinity-v2.5.1/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta >  Trinity.fasta.gene_trans_map

# get Trinity.sqlite
wget "https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v3.sqlite.gz" -O Trinotate.sqlite.gz
gunzip Trinotate.sqlite.gz

Trinotate Trinotate.sqlite init \
	--gene_trans_map Trinity.fasta.gene_trans_map \
	--transcript_fasta Trinity.fasta \
	--transdecoder_pep longest_orfs.pep

#### 	Loading into sqlite db
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
Trinotate Trinotate.sqlite LOAD_signalp signalp.out
Trinotate Trinotate.sqlite LOAD_rnammer Trinity.fasta.rnammer.gff


#### 	Make xls report
Trinotate Trinotate.sqlite report --incl_pep --incl_trans > trinotate_annotation_report.xls

exit 0

