#!/bin/bash
# update : 2024.10.16 sindy@hunbiome.com

import pandas as pd, numpy as np
import os, re, qiime2, biom
from qiime2.plugins import (composition,taxa, composition,emperor,longitudinal,feature_table,fragment_insertion, diversity, demux, quality_control,alignment, demux, emperor, phylogeny, feature_classifier)
from qiime2 import (Artifact, Metadata, Visualization)
from pandas.api.types import CategoricalDtype
from datetime import date
import glob
from skbio import DistanceMatrix
import sys


def create_dir(path):
	if os.path.exists(path)==False:
		os.makedirs(path)
	return 0

def write_excel(excel_path, df, name):
	if not os.path.exists(excel_path):
		with pd.ExcelWriter(excel_path, mode='w', engine='openpyxl') as writer:
			df.to_excel(writer, index=True, sheet_name=name)
	else:
		with pd.ExcelWriter(excel_path, mode='a', engine='openpyxl') as writer:
			df.to_excel(writer, index=True, sheet_name=name)
	return 0

def rename_taxon(li_taxon, str_split):
	li_new_taxon = []
	for feature in li_taxon:
		if ("__un" not in feature.split(str_split)[-1])&("_metagenome" not in feature.split(str_split)[-1])&("__human_" not in feature.split(".")[-1])&("__gut" not in feature.split(".")[-1])&("_gut" not in feature.split(".")[-1]):
			new_feature = feature.split(str_split)[-1]
		else:           
			index_i = len(list(filter(lambda x: x if ("__un" not in x)&("__metagenome" not in x)&("__human_" not in x)&("__gut" not in x)&("_gut" not in x) else '', feature.split(str_split))))-1
			new_feature = str_split.join([feature.split(str_split)[index_i],feature.split(str_split)[-1]])
		li_new_taxon.append(new_feature)
	return (li_new_taxon)

def remove_file(filename):
	if os.path.exists(filename) == True:
		os.remove(filename)

def main():
	try:
		project_pth = sys.argv[1]
		Metagenome_out = sys.argv[2]
		script_pth=os.path.dirname(os.path.abspath(__file__))

		create_dir(Metagenome_out)
		stats_file = os.path.join(Metagenome_out, "Analysis_stats.xlsx")
		diversity_file = os.path.join(Metagenome_out, "Diversity.xlsx")
		taxon_file = os.path.join(Metagenome_out, "Taxonomy_composition.xlsx")
		krona_pth = os.path.join(Metagenome_out, "Taxonomy_composition.krona_plot")
		core_plt_pth = os.path.join(Metagenome_out, "Taxonomy_composition.core-microbiome_in_samples")
		
		remove_file(stats_file)
		remove_file(diversity_file)
		remove_file(taxon_file)
		
		create_dir(krona_pth)
		create_dir(core_plt_pth)
		
		#if os.path.exists(stats_file) == True:
		#	os.remove(stats_file)

		# analysis stats + bacteria filtering stats = stats_out.txt
		df_dada2_stats = pd.read_csv(os.path.join(project_pth,"stats.tsv"), sep="\t", dtype={"sample-id":str})
		df_tax_stats = pd.read_csv(os.path.join(project_pth, "res.taxonomy_composition/taxonomy.table.stats.txt"), sep="\t", dtype={"sample-id":str})
		df_stats = pd.merge(df_dada2_stats, df_tax_stats, on="sample-id", how = 'left')
		df_stats.to_csv(os.path.join(project_pth, "stats_out.txt"), sep="\t", index=False)

		print(df_stats.head())
		write_excel(stats_file, df_stats, "analysis.stats")
	
		li_samples = df_stats[df_stats['analysis']=="success"]["sample-id"].tolist()
		alpha_filename = os.path.join(project_pth,'res.alpha-diversity.normalization.txt')
		dic_metric = {'observed_features':"Observed features", 'chao1':"Chao1 index", 
                      'shannon_entropy':"Shannon's index", 'simpson':"Simpson's index", 'pielou_evenness':"Pielou's evenness"}
		df_alpha = pd.read_csv(alpha_filename,sep="\t", engine="python",dtype={"sample-id":str}).rename(columns = dic_metric)[["sample-id"]+list(dic_metric.values())]
		write_excel(diversity_file, df_alpha[df_alpha["sample-id"].isin(li_samples)], "alpha.normalization")
		dic_metric = {"braycurtis":"Bray-Curtis", "unweighted_unifrac":"unwighted_UniFrac"}
		for metric_name, metric_rename in dic_metric.items():
			q2_dm = Artifact.load(os.path.join(project_pth, "res.beta-diversity", metric_name+"_dm.qza"))
			df_dm = q2_dm.view(DistanceMatrix).to_data_frame()
			write_excel(diversity_file, df_dm.loc[li_samples, li_samples], "beta."+metric_name+"_distance")

		#if os.path.exists(taxon_file):
		#	os.remove(taxon_file)

		dic_name = {"1":"Kingdom","2":"Phylum", "3":"Class", "4":"Order", "5":"Family","6":"Genus", "7":"Species"}
		core_plt_input_file = os.path.join(project_pth, "res.taxonomy_composition/taxonomy-levelAll.RelAbundance.tsv")
		stats_file = os.path.join(project_pth, "stats_out.txt")
		for lev, lev_name in dic_name.items():
			df_relev = pd.read_csv(os.path.join(project_pth, "res.taxonomy_composition/taxonomy-level"+str(lev)+".RelAbundance.tsv"), sep="\t", index_col=0)
			df_s_relev = df_relev[li_samples]
			df_s_relev = df_s_relev[df_s_relev.sum(axis=1)>0]
			if any(df_s_relev.index.str.contains(";__unclassification"))==False :
				df_s_relev.index = df_s_relev.index.str.replace(";__", ";__unclassification")
			df_s_relev = df_s_relev*100/df_s_relev.sum()
			df_prevalence = pd.DataFrame({"sample_found":(df_s_relev>0).sum(axis=1)*100/len(li_samples)})
			df_mean = df_s_relev.T.reset_index().rename(columns = {"index":"sample-id"}).melt(id_vars="sample-id").groupby(["Taxon"]).mean().rename(columns = {"value":"All.MEAN"})
			df_relev_summary = df_prevalence.merge(df_mean, left_index=True, right_index=True).merge(df_s_relev, left_index=True, right_index=True).reset_index().rename(columns = {"Taxon":"Lineage_full"})
			df_relev_summary.insert(1, "Taxon", rename_taxon(df_relev_summary["Lineage_full"].tolist(), ";"))
			mode = "w" if lev=="1" else "a"
			with pd.ExcelWriter(taxon_file,engine='openpyxl', mode=mode) as writer:
				df_relev_summary.to_excel(writer, sheet_name=lev_name, index=False)
			
			command = "Rscript %s/taxonomy_samples_stackedBarPlot.R -i %s -m %s -n analysis -s success -o %s -l %s" %(script_pth, core_plt_input_file, stats_file, core_plt_pth, lev)
			os.system(command)
		
		#Krona plot
		df_abun = pd.read_csv(os.path.join(project_pth, "res.taxonomy_composition/taxonomy-level7.Abundance.tsv"), sep="\t", index_col=0)
		df_abun = df_abun[li_samples]
		if any(df_abun.index.str.contains(";__unclassification"))==False:
			df_abun.index = df_abun.index.str.replace(";__", ";__unclassification")
		df_abun = df_abun.reset_index().rename(columns = {"Taxon":"Lineage_full"})
		df_abun.insert(1, "Taxon", rename_taxon(df_abun["Lineage_full"].tolist(), ";"))
		for colname in li_samples:
			krona_input_file = os.path.join(krona_pth, colname+"_abun.tsv")
			krona_output_file = os.path.join(krona_pth, colname+"_krona.html")
			df_krona_input = df_abun[["Lineage_full", "Taxon", colname]]
			df_krona_input["Taxon"] = df_krona_input["Taxon"].str.replace("\w__", "").str.replace(";__unclassification", "_unclassification").tolist()
			df_krona_input[list(dic_name.values())] = df_krona_input["Lineage_full"].apply(lambda x: pd.Series(str(re.sub("\w__", "", x)).split(";")))
			df_krona_input["Species"] = df_krona_input["Taxon"]
			df_krona_input = df_krona_input.iloc[::, 2::]
			df_krona_input.to_csv(krona_input_file, sep="\t", index=False, header=False)
			command = "ktImportText "+krona_input_file+" -o "+krona_output_file
			os.system(command)

			if os.path.exists(krona_output_file):
				os.remove(krona_input_file)
		
	except Exception as ex:
		print("\n\n################## ERROR [V34-region] ###################")
		print(ex)
		print(">> python default.q2-pipeline_report.py metagenome_analysis_path metagenome_report_path")
if __name__=="__main__":
	main()
