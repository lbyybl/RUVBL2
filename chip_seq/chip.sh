#!/bin/bash

# Boyuan_Li

# Tue Feb 19 16:05:23 2019


set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -f <fastq1> -r <fastq2> -a <forward_adapt> -A <reverse_adapt> -p <Prefix> -o <outputdir> -b <blacklist> -m <min_length> -x <index> -n <normalize_method> -s <step>"
	echo ""
	echo " -f	file         	[required] fastq1"
	echo ""
	echo " -r	file         	[required] fastq2"
	echo ""
	echo " -a	string       	[required] forward_adapt"
	echo ""
	echo " -A	string       	[required] reverse_adapt"
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
	echo " -h			help"
	echo ""
	echo ""
}

min_length=25
index=/DATA/software/bcbio/genomes/Mmusculus/mm10/bowtie2/mm10
normalize_method=RPKM
blacklist=/DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed

if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "f:r:a:A:p:o:b:m:x:h" optionName
do
	case $optionName in
		f) fastq1="$OPTARG";;
		r) fastq2="$OPTARG";;
		a) forward_adapt="$OPTARG";;
		A) reverse_adapt="$OPTARG";;
		p) Prefix="$OPTARG";;
		o) outputdir="$OPTARG";;
		b) blacklist="$OPTARG";;
		m) min_length="$OPTARG";;
		x) index="$OPTARG";;
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

if [[ $fastq2 = "" ]]; then
	echo "the $fastq2 file is needed "
	exit 1
elif [[ ! -f $fastq2 ]]; then
	echo "$fastq2:   is not found"
	exit 2
fi

if [[ $forward_adapt = "" ]]; then
	echo " the $forward_adapt string is needed "
	exit 1
fi

if [[ $reverse_adapt = "" ]]; then
	echo " the $reverse_adapt string is needed "
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
	 echo "$outputdir:   is not found"
	exit 2
fi

if [[ $blacklist = "" ]]; then
	echo "the $blacklist file is needed "
	exit 1
elif [[ ! -f $blacklist ]]; then
	echo "$blacklist:   is not found"
	exit 2
fi

# ???????????????mapping?????????bigwig???????????????????????????call peak????????????????????????????????????????????????????????????????????????????????????mapping??????bigwig???call peak??????hicpro?????????

# 1. ?????? 2.????????? 3. mapping 4. ???bigwig 5.call peak 

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
	#--- 1. ??????
	date
	fastqc ${fastq1} ${fastq2} -o ${outputdir}/QC
	date

else
	fastq_name=${fastq1##*/}
	fastq_prefix=${fastq_name%%.*}
	if [[ ! -f ${outputdir}/QC/${fastq_prefix}_fastqc.html ]]; then
		mkdir -p ${outputdir}/QC
		#--- 1. ??????
		date
		fastqc ${fastq1} ${fastq2} -o ${outputdir}/QC
		date
	fi
fi

#--- 2. ????????? 
if [[ ! -d ${outputdir}/cutadapt ]]; then
	
	mkdir -p ${outputdir}/cutadapt
	date
	cutadapt -a ${forward_adapt} -A ${reverse_adapt} -q 15,15 --overlap 1 -m ${min_length} -o ${outputdir}/cutadapt/${Prefix}_1.fq -p ${outputdir}/cutadapt/${Prefix}_2.fq ${fastq1} ${fastq2}
	date
	fastqc ${outputdir}/cutadapt/${Prefix}_1.fq ${outputdir}/cutadapt/${Prefix}_2.fq -o ${outputdir}/cutadapt
	date

elif [[ ! -f ${outputdir}/cutadapt/${Prefix}_1_fastqc.html ]]; then
	mkdir -p ${outputdir}/cutadapt
	date
	cutadapt -a ${forward_adapt} -A ${reverse_adapt} -q 15,15 --overlap 1 -m ${min_length} -o ${outputdir}/cutadapt/${Prefix}_1.fq -p ${outputdir}/cutadapt/${Prefix}_2.fq ${fastq1} ${fastq2}
	date
	fastqc ${outputdir}/cutadapt/${Prefix}_1.fq ${outputdir}/cutadapt/${Prefix}_2.fq -o ${outputdir}/cutadapt
	date
fi
#--- 3. mapping
if [[ ! -d ${outputdir}/mapping ]]; then
	mkdir -p ${outputdir}/mapping
	date
	bowtie2 -q -x ${index} -X 2000 -1 ${outputdir}/cutadapt/${Prefix}_1.fq -2 ${outputdir}/cutadapt/${Prefix}_2.fq --rg-id ${Prefix} --rg SM:${Prefix} --rg PL:illumina -p 10 -S ${outputdir}/mapping/${Prefix}.sam 2> ${outputdir}/mapping/${Prefix}_mapping.stat
	date
	samtools sort -@ 10 ${outputdir}/mapping/${Prefix}.sam >> ${outputdir}/mapping/${Prefix}.bam
	date
	#rm -f ${outputdir}/mapping/${Prefix}.sam
	sambamba markdup -r -t 10 ${outputdir}/mapping/${Prefix}.bam ${outputdir}/mapping/${Prefix}_rmdup.bam
	date
elif [[ ! -f ${outputdir}/mapping/${Prefix}.bam ]]; then
	date
	bowtie2 -q -x ${index} -X 2000 -1 ${outputdir}/cutadapt/${Prefix}_1.fq -2 ${outputdir}/cutadapt/${Prefix}_2.fq --rg-id ${Prefix} --rg SM:${Prefix} --rg PL:illumina -p 10 -S ${outputdir}/mapping/${Prefix}.sam 2> ${outputdir}/mapping/${Prefix}_mapping.stat
	date
	samtools sort -@ 10 ${outputdir}/mapping/${Prefix}.sam >> ${outputdir}/mapping/${Prefix}.bam
	date
	#rm -f ${outputdir}/mapping/${Prefix}.sam
	sambamba markdup -r -t 10 ${outputdir}/mapping/${Prefix}.bam ${outputdir}/mapping/${Prefix}_rmdup.bam
	date
fi


#--- 4. ?????? coordinate???unique???
if [[ ! -d ${outputdir}/unique ]]; then
	mkdir -p ${outputdir}/unique
	date
	/home/boyuanli/bashscript/bin/pro_seq/unique_and_bw.sh -b ${outputdir}/mapping/${Prefix}_rmdup.bam -o ${outputdir}/unique/ -p ${Prefix}_rmdup -B /DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed
	date
elif [[ ! -f ${outputdir}/unique/${Prefix}_rmdup_uniqe.bam ]]; then
	date
	/home/boyuanli/bashscript/bin/pro_seq/unique_and_bw.sh -b ${outputdir}/mapping/${Prefix}_rmdup.bam -o ${outputdir}/unique/ -p ${Prefix}_rmdup -B /DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed
	date
fi
	

#--- stastics
#--- total reads

if [[ ! -f ${outputdir}/${Prefix}_stastic_all.csv ]]; then
	if [[ ! -f ${outputdir}/${Prefix}_stastic.csv ]]; then
		pwd
		suffix=$(echo ${fastq1} | awk -F "." '{print $NF}')
		if [[ $suffix == 'gz' ]]; then
			total_resds=$(grep "^@" <(zcat ${fastq1}) | wc -l)
		elif [[ $suffix != 'gz' ]]; then
			total_resds=$(grep "^@" ${fastq1} | wc -l)
		fi

		cutadapt=$(grep "^@" ${outputdir}/cutadapt/${Prefix}_1.fq | wc -l)
		mapping=$(samtools view -F 4 ${outputdir}/mapping/${Prefix}.bam | wc -l)
		rmdup=$(samtools view -F 4 ${outputdir}/mapping/${Prefix}_rmdup.bam | wc -l)
		final_uniq=$(head -1 ${outputdir}/unique/${Prefix}_rmdup_uniqe.flag | awk '{print $1}')
		#---- stastic file
		echo -e "sample_name,total_resds,cutadapt,mapping reads,rm duplicate,final unique" > ${outputdir}/${Prefix}_stastic.csv
		echo -e "${Prefix},${total_resds},${cutadapt},${mapping},${rmdup},${final_uniq}" >> ${outputdir}/${Prefix}_stastic.csv
	fi
	sed -e "s/^[ \t]*//g" -e 's/ /,/g' ${outputdir}/mapping/${Prefix}_mapping.stat >> ${outputdir}/mapping/${Prefix}_mapping_format.stat
	/usr/bin/Rscript /home/boyuanli/bashscript/bin/chip_seq/stastic_reads_info.r -b ${outputdir}/mapping/${Prefix}_mapping_format.stat -s ${outputdir}/${Prefix}_stastic.csv -o ${outputdir}/${Prefix}_stastic_all.csv
	rm -f ${outputdir}/mapping/${Prefix}.sam ${outputdir}/mapping/${Prefix}_rmdup.bam
	rm -f ${outputdir}/cutadapt/${Prefix}_1.fq ${outputdir}/cutadapt/${Prefix}_2.fq ${outputdir}/${Prefix}_stastic.csv 
	# rm -rf unique
fi
date
/home/boyuanli/bashscript/bin/chip_seq/chip_qc.sh -b ${outputdir}/mapping/${Prefix}.bam -n ${outputdir}/${Prefix}_qc.csv
date
if [[ -f ${outputdir}/${Prefix}_stastic_all.csv ]]; then
	if [[ ! -f ${outputdir}/${Prefix}_qc.csv ]]; then
		#--- ??????bigwig		
		rm -f ${outputdir}/mapping/${Prefix}.sam ${outputdir}/mapping/${Prefix}_rmdup.bam
		rm -f ${outputdir}/cutadapt/${Prefix}_1.fq ${outputdir}/cutadapt/${Prefix}_2.fq ${outputdir}/${Prefix}_stastic.csv 
		# rm -rf unique
	else
		rm -f ${outputdir}/mapping/${Prefix}.sam ${outputdir}/mapping/${Prefix}_rmdup.bam
		rm -f ${outputdir}/cutadapt/${Prefix}_1.fq ${outputdir}/cutadapt/${Prefix}_2.fq ${outputdir}/${Prefix}_stastic.csv 
		# rm -rf unique
	fi
fi
date

