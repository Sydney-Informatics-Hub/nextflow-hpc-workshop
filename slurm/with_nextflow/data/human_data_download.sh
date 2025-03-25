#!/bin/bash -l
#SBATCH --account=$PAWSEY_PROJECT
#SBATCH --job-name=human_data
#SBATCH --output=human_data_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --partition=copy
#SBATCH --time=01:00:00
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END,FAIL

##Data origin
##Manuscript Antagonizing cholecystokinin A receptor in the lung attenuates obesity-induced airway hyperresponsiveness https://www.nature.com/articles/s41467-022-35739-8
##RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine function in airway smooth muscle cells https://pubmed.ncbi.nlm.nih.gov/24926665/

mkdir human
cd human
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039509/SRR1039509_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039509/SRR1039509_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039510/SRR1039510_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039510/SRR1039510_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039511/SRR1039511_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039511/SRR1039511_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039514/SRR1039514_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039514/SRR1039514_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039517/SRR1039517_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039517/SRR1039517_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039518/SRR1039518_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039518/SRR1039518_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039519/SRR1039519_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039519/SRR1039519_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039520/SRR1039520_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039520/SRR1039520_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039521/SRR1039521_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039521/SRR1039521_2.fastq.gz
