#!/bin/bash

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -f <fastq1> -r <fastq2> -a <forward_adapt> -A <reverse_adapt> -p <Prefix> -o <outputdir> -b <blacklist> -m <min_length> -x <index> -n <normalize_method> -s <step>"
	echo ""
	echo " -f	file         	[required] fastq1"
	echo ""
	echo " -a	string       	[required] forward_adapt"
	echo ""
	echo " -p	string       	[required] output prefix"
	echo ""
	echo " -o	dir          	[required] output dir"
	echo ""
	echo " -b	file         	[optional] blacklist difult:/DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed"
	echo ""
	echo " -m	number         	[optional] min_length, defult:25"
	echo ""
	echo " -x	path and Prefix [optional] index using for mapping defult:/DATA/software/bcbio/genomes/Mmusculus/mm10/bowtie2/mm10"
	echo ""
	echo " -n	string         	[optional] normalize_method, defult:RPKM;RPKM, CPM, BPM, RPGC, None"
	echo ""
	echo " -s	string         	[optional] which step you want to run: QC,cutadapt,mapping,bigwig,callpeak"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}

min_length=25
#index=/DATA/software/bcbio/genomes/Mmusculus/mm10/bowtie2/mm10
#index=/DATA/genomes/mm10_ofr_3/mm10
#normalize_method=RPKM
#blacklist=/DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed

if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "f:a:p:o:b:m:x:n:s:h" optionName
do
	case $optionName in
		f) fastq1="$OPTARG";;
		a) forward_adapt="$OPTARG";;
		p) Prefix="$OPTARG";;
		o) outputdir="$OPTARG";;
		b) blacklist="$OPTARG";;
		m) min_length="$OPTARG";;
		x) index="$OPTARG";;
		n) normalize_method="$OPTARG";;
		s) step="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $fastq1 = "" ]]; then
	echo "the $fastq1 file is needed "
	exit 1
elif [[ ! -f $fastq1 ]]; then
	echo "$fastq1:   is not found"
	exit 2
fi


if [[ $forward_adapt = "" ]]; then
	echo " the $forward_adapt string is needed "
	exit 1
fi


if [[ $Prefix = "" ]]; then
	echo " the $Prefix string is needed "
	exit 1
fi

if [[ $outputdir = "" ]]; then
	echo "the $outputdir file is needed "
	exit 1
elif [[ ! -d $outputdir ]]; then
	#echo "$outputdir:   is not found"
	#exit 2
	mkdir -p $outputdir
fi

# 写一个能够mapping，转成bigwig文件，还能去接头，call peak的脚本，这个脚本还应该有模块化的功能比如，去接头与质控，mapping与转bigwig，call peak就像hicpro一样；

# 1. 质控 2.去接头 3. mapping 4. 转bigwig 5.call peak 

# fastq1
# fastq2
# forward_adapt
# reverse_adapt
# min_length # defult 25
# Prefix
# index # defult is /DATA/software/bcbio/genomes/Mmusculus/mm10/bowtie2/mm10
# normalize_method # defult is RPKM
# blacklist # defult is /DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed
# outputdir

#mkdir -p ${outputdir}/{bigwig,cutadapt,mapping,QC,rawdata}
if [[ ! -d ${outputdir}/QC ]]; then
	mkdir -p ${outputdir}/QC
	#--- 1. 质控
	fastqc ${fastq1} -o ${outputdir}/QC
fi
#--- 2. 去接头 
if [[ ! -d ${outputdir}/cutadapt ]]; then
	mkdir -p ${outputdir}/cutadapt
	cutadapt -a ${forward_adapt} -q 15 --overlap 1 -m ${min_length} -o ${outputdir}/cutadapt/${Prefix}_1.fq ${fastq1}
	fastqc ${outputdir}/cutadapt/${Prefix}_1.fq -o ${outputdir}/cutadapt
	gzip ${outputdir}/cutadapt/${Prefix}_1.fq
fi
#--- 3. mapping
if [[ ! -d ${outputdir}/mapping ]]; then
	mkdir -p ${outputdir}/mapping

	/usr/bin/STAR --genomeDir /DATA3/genomes/STAR_index/mm10 --readFilesIn ${outputdir}/cutadapt/${Prefix}_1.fq.gz --readFilesCommand zcat --outFileNamePrefix ${outputdir}/mapping/${Prefix}_ --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --peOverlapNbasesMin 20 --twopassMode Basic --outFilterMultimapNmax 3 --runThreadN 2
fi

if [[ ! -d ${outputdir}/mapping_qc ]]; then
	samtools view -bh -F 260 -q 30 ${outputdir}/mapping/${Prefix}_Aligned.sortedByCoord.out.bam -o ${outputdir}/mapping/${Prefix}_primary.bam
	/DATA/software/bcbio/tools/bin/sambamba markdup -r -t 2 ${outputdir}/mapping/${Prefix}_primary.bam ${outputdir}/mapping/${Prefix}_rmdup.bam
	
	mkdir ${outputdir}/mapping_qc
	bam_stat.py -i ${outputdir}/mapping/${Prefix}_rmdup.bam > ${outputdir}/mapping_qc/bam_stat.txt
	split_bam.py -i ${outputdir}/mapping/${Prefix}_rmdup.bam -r /DATA3/genomes/annotations/mm10_rRNA.bed -o ${outputdir}/mapping_qc/${Prefix}_rRNA
	read_distribution.py -i ${outputdir}/mapping/${Prefix}_rmdup.bam -r /DATA3/genomes/annotations/mm10.RefSeq.union.bed > ${outputdir}/mapping_qc/reads_distr.txt
	read_GC.py -i ${outputdir}/mapping/${Prefix}_rmdup.bam -o ${outputdir}/mapping_qc/gc_distr

fi


#if [[ -f ${outputdir}/${Prefix}_statistic.csv ]]; then
if [[ ! -f ${outputdir}/${Prefix}_statistic.csv ]]; then
        echo "calculating mapping stats..."
        pwd
        suffix=$(echo ${fastq1} | awk -F "." '{print $NF}')
        if [[ $suffix == 'gz' ]]; then
                total_reads=$(grep "^@" <(zcat ${fastq1}) | wc -l)
        elif [[ $suffix != 'gz' ]]; then
                total_reads=$(grep "^@" ${fastq1} | wc -l)
        fi

        cutadapt=$(grep "^@" <(zcat ${outputdir}/cutadapt/${Prefix}_1.fq) | wc -l)
        mapping=$(samtools view -c -@ 3 -F 260 ${outputdir}/mapping/${Prefix}_Aligned.sortedByCoord.out.bam)
	rmdup=$(samtools view -F 4 ${outputdir}/mapping/${Prefix}_rmdup.bam | wc -l)
	proper_pair=$(grep "proper" ${outputdir}/mapping_qc/bam_stat.txt | grep -o "[0-9]*$")

#--- stastic file
	echo -e "sample_name,total_reads,cutadapt,mapping reads,rm duplicate,proper_pair" >> ${outputdir}/${Prefix}_statistic.csv
    echo -e "${Prefix},${total_reads},${cutadapt},${mapping},${rmdup},${proper_pair}" >> ${outputdir}/${Prefix}_statistic.csv
fi
