#/bin/bash 2024-03-19 sindy@hunbiome.com

if [ "$#" -lt "5" -o "$#" ]; then
        echo "##### check the parameters #####"
        echo "run_script.amplicon_v34.basic_analysis.sh <output_path> <metadata file> <metadata column> <group value> <paired TRUE|FALSE>"
fi


out_pth=${1}
ori_metafile=${2} #${out_pth}/metadata.240320.txt
grp=${3} #metadata group column
sort_grp=${4} #'Control_Female,Control_Male,Patients_Female,Patients_Male'
paired=${5} # TRUE or FALSE
group_n=`head -n1 ${ori_metafile} | tr "\t" "\n" | grep -n ${grp} | awk -F":" '{print $1-1}'`
script_pth=${CONDA_PREFIX}/opt/metagenome_amplicon/
date=`date "+%Y-%m-%d"`
comp_pth=${out_pth}/${date}.comparison_analysis
grp_prefix_pth=`echo $sort_grp | sed 's/,/__vs__/g'`
grp_comp_pth=${comp_pth}/${grp_prefix_pth}

prev_cutoff=30
avg_relab_cutoff=0.001

mkdir ${comp_pth}
lefse_pth=${grp_comp_pth}/LEfSe

#cp ${metafile} ${comp_pth}
metafile=${comp_pth}/analysis_pass.metadata.txt
perl ${script_pth}/table.addColumns.pl -i stats_out.txt 0 ${ori_metafile} 0 ${group_n} > ${metafile}

echo -e "Calculate the diversity score"
echo -e ">> Alpha diversity"
Rscript ${script_pth}/alpha_boxplot_st.R -i ${out_pth}/res.alpha-diversity.normalization.txt -m ${metafile} -n ${grp} -s ${sort_grp} -o ${grp_comp_pth} -p ${paired}

echo -e ">> Beta diversity"
Rscript ${script_pth}/beta_pcoa_st.R -i ${out_pth}/res.beta-diversity/braycurtis_dm.txt -m ${metafile} -n ${grp} -s ${sort_grp} -o ${grp_comp_pth} -t 'Bray-Curtis'
Rscript ${script_pth}/beta_pcoa_st.R -i ${out_pth}/res.beta-diversity/unweighted_unifrac_dm.txt -m ${metafile} -n ${grp} -s ${sort_grp} -o ${grp_comp_pth} -t 'unweighted UniFrac'
echo -e ".done"

echo -e "Calculate the taxonomy composition"
echo -e ">> statistical analysis"
python ${script_pth}/taxonomy_st.py -i ${out_pth}/res.taxonomy_composition/taxonomy-levelAll.RelAbundance.tsv -m ${metafile} -n ${grp} -s ${sort_grp} -o ${grp_comp_pth}/taxonomy_composition/ -p ${paired}
echo -e ">> taxonomy composition; stacked barplot - groups"
echo -e "  - prelavence > ${prev_cutoff}%"
echo -e "  - all-average relative > ${avg_relab_cutoff}%"

Rscript ${script_pth}/taxonomy_groups_stackedBarPlot.R -i ${grp_comp_pth}/taxonomy_composition/Bacteria.statstical.xlsx -s ${sort_grp} -o ${grp_comp_pth}/taxonomy_composition -r ${avg_relab_cutoff} -p ${prev_cutoff}

echo -e ".done"

echo -e "Calculate the taxonomy features"
echo -e ">> LEfSe analysis"
python ${script_pth}/lefse_intput.py -i ${out_pth}/res.taxonomy_composition/taxonomy-levelAll.RelAbundance.tsv -m ${metafile} -n ${grp} -s ${sort_grp} -o ${lefse_pth}/lefse_in.txt

chmod 777 ${lefse_pth}

docker run -v ${lefse_pth}/:/home/linuxbrew/output biobakery/lefse format_input.py /home/linuxbrew/output/lefse_in.txt /home/linuxbrew/output/lefse.in -c 2 -u 1 -o 1000000
docker run -v ${lefse_pth}/:/home/linuxbrew/output biobakery/lefse run_lefse.py /home/linuxbrew/output/lefse.in /home/linuxbrew/output/lefse.res


#source /data/sindy/miniconda3/bin/activate /data/neville/miniconda3/envs/python3.6/lefse_cosmaxbti
#/data/neville/miniconda3/envs/python3.6/lefse_cosmaxbti/bin/lefse-format_input.py ${lefse_pth}/lefse_in.txt ${lefse_pth}/lefse.in -c 2 -u 1 -o 1000000
#/data/neville/miniconda3/envs/python3.6/lefse_cosmaxbti/bin/run_lefse.py ${lefse_pth}/lefse.in ${lefse_pth}/lefse.res
echo -e "full_lineage\thigh_group\tLDA_score\tp-value" > ${lefse_pth}/lefse_out.txt
cut -f1,3- ${lefse_pth}/lefse.res | grep -vP "\t-$" | grep -v  ".__unclassification" >> ${lefse_pth}/lefse_out.txt

echo -e ">> LEfSe; features barplot"
source /data/sindy/miniconda3/bin/activate /data/sindy/miniconda3/envs/qiime2-amplicon-2023.9
Rscript ${script_pth}/lefse_barplot.R -i ${lefse_pth}/lefse_out.txt -s ${sort_grp}
for level in `seq 2 7`;do Rscript ${script_pth}/lefse_barplot.R -i ${lefse_pth}/lefse_out.txt -s ${sort_grp} -r $level;done

rm ${lefse_pth}/lefse.in ${lefse_pth}/lefse.res ${lefse_pth}/lefse_in.txt
echo -e ".done"
