# Author : Kyeongeui Yun (sindy@hunbiome.com)
#/bin/bash
if [ "$#" -lt "5" -o "$#" -gt "6" ]; then
	echo "##### Prepare reseach samples for qiime2 pipeline #####"
	echo "usage : prepare_samples_for_qiime2.sh qiime2_dir reads_dir first_suffix second_suffix type(1 or 2) (customed:sample-id.name)name"
	echo "type 1 : reads_dir/Sample_name/FASTQ"
	echo "type 2 : reads_dir/FASTQ"
	echo "ex) prepare_samples_for_qiime2.sh /work/directory/ /sample/directory _1.fastq.gz _2.fastq.gz 2"
	echo "typing is ::: $#"
	exit 1
fi
if	[ ! -d ${1} ]; then
	mkdir ${1}
fi

if [ $# -eq "6" ]; then
	sample_add=.${6}
fi
if [ $# -eq "5" ]; then
	sample_add=""
fi


MANIFEST=${1}/manifest${sample_add}.txt
output_path=${1}
data_path=${2}
echo -e "\n\n============ Prepareing the QIIME2 ========="
echo "[ input length : "$#" / sample_add  : "${sample_add}" ]"
echo "Metagenome 16s rRNA path : ${output_path}"
echo "FASTQ input path : ${data_path}"
echo "Input Type : ${5}"
echo -e "MANIFEST FILE : ${MANIFEST}\n"

R1_suffix=${3}
R2_suffix=${4}

echo -e "======= Sample_dir/FASTQ format ======\n"
if [ ${5} = "1" ];then
	echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > $MANIFEST
	for filename in ${data_path}/*/*${R1_suffix};do
		sample_name=`echo ${filename##*/} | sed "s/${R1_suffix}//g"`
		sample_name=${sample_name}${sample_add}
		ch_sample_id=`echo ${sample_name} | sed "s/-/-/"`
		R1=${filename}
		R2=`echo ${filename} | sed "s/${R1_suffix}/${R2_suffix}/g"`
		echo -e ${ch_sample_id}'\t'${R1}'\t'${R2}
	
	done >> $MANIFEST


elif [ ${5} = "2" ];then
	echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > $MANIFEST
	for filename in ${data_path}/*${R1_suffix}; do
		sample_name=`echo ${filename##*/} | sed "s/${R1_suffix}//g"`
		sample_name=${sample_name}${sample_add}
		ch_sample_id=${sample_name}
		R1=${filename}
		R2=`echo ${filename} | sed "s/${R1_suffix}/${R2_suffix}/g"`
		echo -e ${ch_sample_id}'\t'${R1}'\t'${R2}
	done >> $MANIFEST
fi

head ${MANIFEST}
