# 2024. 06. 24 sindy@hunbiome.com
#/bin/bash

if [ "$#" -lt "3" -o "$#" ]; then
	echo "##### check the parameters #####"
	echo "run_script.amplicon_v34.basic_analysis.sh <sequencing_path> <output_path> <sampling-type cut-off>"
fi


date=`date "+%Y-%m-%d"`
sequencing_pth=${1}
out_pth=${2}
min_depth=${3}
cpu=10

report_out_pth=${out_pth}/basic_analysis/
script_pth=${CONDA_PREFIX}/opt/metagenome_amplicon/

# prepare the qiime2 run
echo -e "> Make the manifest file"
echo -e " - sequencing path: ${sequencing_pth}"
echo -e " - output path: ${out_pth}"
bash ${script_pth}/metagenome.amplicon-pe.prepare_samples.sh  ${out_pth}/ ${sequencing_pth} _1.fastq.gz _2.fastq.gz 2
echo -e "done"

# DADA2 denoising
echo -e "> Denoising : DADA2 method"
python ${script_pth}/metagenome.amplicon-pe.figaro_q2-denoise.py ${out_pth} 17 21 430 ${min_depth} ${cpu}
echo -e "done"

# Taxonomy classification - SILVA 138v (V3-V4 extracted DB)
echo -e "> Taxonomy classification"
echo -e " - database: SILVA 138v _ pre-trained"
echo -e " - classification method: sklearn classifier"
python ${script_pth}/metagenome.amplicon-pe.q2-clr.py ${out_pth} 34 ${cpu} q2-2023.9
echo -e "done"

# Taxonomy composition table - Not Exclude; Mitochondria, Chloroplast
echo -e "> Taxonomy filtering"
echo -e " - remove the features : mitochondria & chloroplast (contamination features)"
python ${script_pth}/metagenome.q2-clr_abun-relab.py -i ${out_pth}/table_filt.qza -c ${out_pth}/taxonomy.sklearn.qza -o ${out_pth}/res.taxonomy_composition/
echo -e "done"

# Diversity - table removed contamination features
echo -e "> Diversity"
echo -e " - Alpha diversity: Chao1, Observed features, Shannon's index, Simpson's index, Simpson envenness, Pielou's evenness"
echo -e " - Beta diversity: Bray-curtis, unweighted UniFrac"
python ${script_pth}/metagenome.q2-diversity.py ${out_pth}/representative_sequences.qza ${out_pth}/res.taxonomy_composition/taxonomy.table.qza ${out_pth}/
echo -e "done"

# Basic Report
echo -e "> Report"
ln -sf ${out_pth}/stats_preprocessing.txt ${out_pth}/stats.tsv
python ${script_pth}/metagenome.report.py ${out_pth} ${report_out_pth}
echo -e "done"

