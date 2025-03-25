#!/bin/bash -l
#PBS -N nfrnaseq_test
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -l walltime=24:00:00
cd $PBS_O_WORKDIR
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/044/SRR20630344/SRR20630344.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/049/SRR20630349/SRR20630349.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/055/SRR20630355/SRR20630355.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/047/SRR20630347/SRR20630347.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/050/SRR20630350/SRR20630350.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/042/SRR20630342/SRR20630342.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/053/SRR20630353/SRR20630353.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/043/SRR20630343/SRR20630343.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/039/SRR20630339/SRR20630339.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/056/SRR20630356/SRR20630356.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/054/SRR20630354/SRR20630354.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/041/SRR20630341/SRR20630341.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/045/SRR20630345/SRR20630345.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/051/SRR20630351/SRR20630351.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/040/SRR20630340/SRR20630340.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/048/SRR20630348/SRR20630348.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/052/SRR20630352/SRR20630352.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR206/046/SRR20630346/SRR20630346.fastq.gz
