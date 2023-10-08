#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name="10.sam2bam"
#SBATCH --output="10.sam2bam.o"
#SBATCH --error="10.sam2bam.e"
#SBATCH --constraint="ethernet"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=20
#SBATCH --time=0-72:00:00
#SBATCH --account="generalGrp"
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=saundecj@wfu.edu
#SBATCH --mem=128GB

###Now do your stuff

cd /deac/generalGrp/johnsonGrp/saundecj/projects/eiseniaHortensis.20180302

module load samtools/1.5-icc

echo "view eh.head.S1"
samtools view -Sb eh.head.S1.sam > eh.head.S1.unsorted.bam
echo "sort"
samtools sort eh.head.S1.unsorted.bam -o eh.head.S1.bam
echo "index"
samtools index eh.head.S1.bam

echo "view eh.head.S2"
samtools view -Sb eh.head.S2.sam > eh.head.S2.unsorted.bam
echo "sort"
samtools sort eh.head.S2.unsorted.bam -o eh.head.S2.bam
echo "index"
samtools index eh.head.S2.bam

echo "view eh.head.S3"
samtools view -Sb eh.head.S3.sam > eh.head.S3.unsorted.bam
echo "sort"
samtools sort eh.head.S3.unsorted.bam -o eh.head.S3.bam
echo "index"
samtools index eh.head.S3.bam

echo "view eh.mid.S4"
samtools view -Sb eh.mid.S4.sam > eh.mid.S4.unsorted.bam
echo "sort"
samtools sort eh.mid.S4.unsorted.bam -o eh.mid.S4.bam
echo "index"
samtools index eh.mid.S4.bam

echo "view eh.mid.S5"
samtools view -Sb eh.mid.S5.sam > eh.mid.S5.unsorted.bam
echo "sort"
samtools sort eh.mid.S5.unsorted.bam -o eh.mid.S5.bam
echo "index"
samtools index eh.mid.S5.bam

echo "view eh.mid.S6"
samtools view -Sb eh.mid.S6.sam > eh.mid.S6.unsorted.bam
echo "sort"
samtools sort eh.mid.S6.unsorted.bam -o eh.mid.S6.bam
echo "index"
samtools index eh.mid.S6.bam

echo "finished"
module unload samtools/1.5-icc
exit 0
