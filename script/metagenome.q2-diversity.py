#!/bin/bash
# update : 2023.03
import os
import pandas as pd
import sys
from qiime2.plugins import (dada2, diversity, feature_table,phylogeny)
from qiime2 import (Artifact, Metadata)
import time
import datetime
import re
from skbio import DistanceMatrix
import skbio
import json
def DirPath(path_loc):
	dir_path = path_loc
	if(dir_path[-1]!="/"):
		dir_path = dir_path+"/"
	return dir_path

def Alpha(q2_table):
	li_metric = ["chao1", "observed_features", "shannon", "simpson", "pielou_e"] #Alpha diversity 측정 metrics 
	dic_metric = {'observed_features' : "Observed features",'chao1': "Chao1 index", 'shannon_entropy':"Shannon's index", 'simpson':"Simpson's index", 'pielou_evenness':"Pielou's evenness"}
	li_alpha_df = []
	li_alpha_rare = []
	print(">> Alpha diversity")
	for metric in li_metric:
		df = pd.DataFrame(diversity.pipelines.alpha(q2_table, metric).alpha_diversity.view(pd.Series)) # Alha diversity 측정
		li_alpha_df.append(df)
	df_alpha = pd.concat(li_alpha_df, axis=1)
	df_alpha.index.name = "sample-id"
	df_alpha = df_alpha[dic_metric.keys()].rename(columns = dic_metric)
	return df_alpha

def Phylogeny(q2_repseq, q2_table_norm, phylogeny_path):
	q2_phylogeny = phylogeny.pipelines.align_to_tree_mafft_fasttree(q2_repseq, n_threads='auto')
	q2_phylogeny.alignment.save(os.path.join(phylogeny_path, "aligned-rep-seqs.qza"))
	q2_phylogeny.masked_alignment.save(os.path.join(phylogeny_path, "masked-aligned-rep-seqs.qza"))
	q2_phylogeny.tree.save(os.path.join(phylogeny_path, "unrooted-tree.qza"))
	q2_phylogeny.rooted_tree.save(os.path.join(phylogeny_path, "rooted-tree.qza"))
	return

def Beta(q2_table, beta_path, phylogeny_path):
	li_metric = ["braycurtis", "unweighted_unifrac"]
	for metric_name in li_metric:
		if metric_name =="braycurtis":
			q2_dm = diversity.pipelines.beta(table = q2_table, metric =metric_name).distance_matrix
		elif metric_name =="unweighted_unifrac":
			q2_rooted_tree = Artifact.load(os.path.join(phylogeny_path, "rooted-tree.qza"))
			q2_dm = diversity.pipelines.beta_phylogenetic(table = q2_table, phylogeny = q2_rooted_tree, metric = metric_name).distance_matrix
		filename = os.path.join(beta_path, metric_name+"_dm")
		q2_dm.save(filename)
		dm_df = q2_dm.view(DistanceMatrix).to_data_frame()
		dm_df.reset_index().rename(columns = {"index":"distance_matrix"}).to_csv(filename+".txt", sep="\t",index=False)

	return

def write_excel(excel_path, df, name):
	if not os.path.exists(excel_path):
		with pd.ExcelWriter(excel_path, mode='w', engine='openpyxl') as writer:
			df.to_excel(writer, index=True, sheet_name=name)
	else:
		with pd.ExcelWriter(excel_path, mode='a', engine='openpyxl') as writer:
			df.to_excel(writer, index=True, sheet_name=name)
	return

def main():
	try:
		rep_seq_filename = sys.argv[1]
		table_filename = sys.argv[2]
		output_path = sys.argv[3]
		
		phylogeny_path = os.path.join(output_path, "res.phylogeny")
		beta_path = os.path.join(output_path, "res.beta-diversity")
		

		if os.path.exists(output_path)==False:
			os.makedirs(output_path)
		
		q2_repseq = Artifact.load(rep_seq_filename)
		q2_table = Artifact.load(table_filename)
		norm_depth  = int(q2_table.view(pd.DataFrame).sum(axis=1).min())
		
		
		dic_run = {"Input.representative_sequences_fileName":rep_seq_filename, "Input.table_fileName": table_filename, "OutputPath":output_path}
		if(os.path.exists(os.path.join(output_path, "res.alpha-diversity.normalization.txt"))==False):
			df_alpha  = Alpha(q2_table)
			df_alpha.to_csv(os.path.join(output_path,"res.alpha-diversity.tsv"), sep="\t", index=True)
			print("No-mito-no-chloro min depth : ", norm_depth)
			if norm_depth >0:
				q2_table_norm  = feature_table.methods.rarefy(q2_table , sampling_depth = norm_depth).rarefied_table
				q2_table_norm.save(os.path.join(output_path,"table_rarefied.qza"))
				df_alpha_norm  = Alpha(q2_table_norm)
				df_alpha_norm.to_csv(os.path.join(output_path, "res.alpha-diversity.normalization.txt"), sep="\t", index=True)	
				dic_run["Res.Normalization_depth"] = str(norm_depth)
				dic_run["Res.table_rafied"] = "TRUE"
				dic_run["Res.Alpha-diversity_normalization"] = "TRUE"

		if(os.path.exists(phylogeny_path)==False):
			os.makedirs(phylogeny_path)
			#q2_table_norm  = Artifact.load(os.path.join(output_path,"table_rarefied.qza"))
			q2_repseq_filt = feature_table.methods.filter_seqs(data = q2_repseq, table = q2_table_norm).filtered_data
			Phylogeny(q2_repseq_filt , q2_table_norm, phylogeny_path)
			dic_run["Res.Phylogeny_OutputPath"] = phylogeny_path
			dic_run["Res.Phylogeny_rarefiedInput"] = "TRUE"
			
		if(os.path.exists(beta_path)==False):
			os.makedirs(beta_path)
			Beta(q2_table_norm, beta_path, phylogeny_path)

			dic_run["Res.Beta-diversity_OutputPath"] = beta_path
			dic_run["Res.Beta-diversity_rarefiedInput"] = "TRUE"
		with open(os.path.join(output_path, "res.diversity_info.txt"), 'w') as f:
			for key, value in dic_run.items():
				f.write('%s:%s\n' % (key, value))
        


	except Exception as ex:
		print(">> python default.q2-pipeline_4.py final_representative_sequence final_table result_output_path")
		print("Error >>>>>>>> ")
		print(ex)

if __name__ == "__main__":
	main()
