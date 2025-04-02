#!/bin/bash
# Author : Kyeongeui Yun (fiow3250@gmail.com)
# update : 2020.11
import pandas as pd
import os
import numpy as np
import sys
from qiime2.plugins import (dada2, diversity, feature_table, feature_classifier, taxa, cutadapt)
from qiime2 import (Artifact, Metadata)
import time
import datetime
import re
sys.path.append(os.path.join(os.environ.get("CONDA_PREFIX"), "opt/figaro/figaro/"))
from figaro import figaro
file_name = "Figaro_default.q2-pipeline_1.py"
code_time = time.time()

now=datetime.datetime.now()
nowDate=now.strftime('%Y%m%d')

def DirPath(path_loc):
	dir_path = path_loc
	if(dir_path[-1]!="/"):
		dir_path = dir_path+"/"
	return dir_path

def Time_sc(now, start):
	time_sc = time.gmtime(now-start)
	print(">> Run Time : %d h %d m %s s" %(time_sc.tm_hour, time_sc.tm_min, time_sc.tm_sec))
	return

def stats_depth(df_stats, filtering_depth):
	df_stats.loc[df_stats[df_stats["non-chimeric"]<filtering_depth].index,"analysis"] = "fail"
	df_stats.loc[df_stats[df_stats["non-chimeric"]>=filtering_depth].index,"analysis"] = "success"
	if (any(df_stats["non-chimeric"]<filtering_depth)==True):
		print(">> Fail samples ... \n\n")
		print(df_stats[df_stats["non-chimeric"]<filtering_depth])
	return df_stats

def main():
	try:
		#data_path = DirPath(sys.argv[1])
		output_path = DirPath(sys.argv[1])
		f_primer_len = int(sys.argv[2])
		r_primer_len = int(sys.argv[3])
		trimming_len = int(sys.argv[4])
		filtering_depth = int(sys.argv[5])
		manifest = output_path+"manifest.txt"
		cpus = int(sys.argv[6]) #cpu maximum 30


		print("\n\n###################################################################################")
		print("Manifest : ",manifest)
		#print("Data path : ",data_path)
		print("Metagenome 16s rRNA output path : ", output_path)
		print("Analaysis R&D Region >> customer")
		if os.path.isdir(output_path)==False:
			os.system("mkdir "+output_path)
		print("###################################################################################")

		try:
			rename_path = output_path+"Rename/"
			os.system("mkdir "+rename_path)
			print("\n\n################# Rename #################")
			print("Rename : Illumina Sequence Name")
			start = time.time()
			df_manifest = pd.read_csv(manifest, sep="\t", dtype=str)
			for i in df_manifest.index.tolist():
				sampleID = df_manifest.loc[i, "sample-id"]
				R1_fq = df_manifest.loc[i, "forward-absolute-filepath"]
				R2_fq = df_manifest.loc[i, "reverse-absolute-filepath"]
				newID = (sampleID+"_S"+str(i+1)+"_L001").replace(".","_")
				os.system("cp "+R1_fq+" "+rename_path+newID+"_R1.fastq.gz")
				os.system("cp "+R2_fq+" "+rename_path+newID+"_R2.fastq.gz")
			Time_sc(time.time(), start)
			
			print("\n\n################# Figaro #################")
			print("Figaro : Researching the dada2 option")
			start = time.time()
			ampliconLength = trimming_len # V34-region : 450
			forwardPrimerLength = f_primer_len #HanLab_V34-region : 17
			reversePrimerLength = r_primer_len #HanLab_V34-region : 20
			sequenceFolder = rename_path
			resultTable, forwardCurve, reverseCurve = figaro.runAnalysis(sequenceFolder, ampliconLength, forwardPrimerLength, reversePrimerLength)
			trunc_f = resultTable[0].forwardTrimPosition
			trunc_r = resultTable[0].reverseTrimPosition
			ee_f = resultTable[0].forwardMaxExpectedError
			ee_r = resultTable[0].reverseMaxExpectedError
			if ee_f>=4:
				ee_f=3
			if ee_r >=4:
				ee_r =3
			if ee_f<2:
				ee_f = 2
			if ee_r<2:
				ee_r = 2
			print("trimming : ", trunc_f, trunc_r)
			print("max_ee : ", ee_f, ee_r)
			df_figaro = pd.DataFrame({"F" : [trunc_f, ee_f], "R":[trunc_r, ee_r]}, index = ["Truncate", "Max Error"])
			df_figaro.to_csv(output_path+"Figaro.txt", sep="\t", index=True)
			Time_sc(time.time(), start)
			
			print("\n\n################# Importing #################")
			print("Importing : FASTQ format file --> .qza QIIME2 format file")
			start = time.time()
			demux_PE = Artifact.import_data('SampleData[PairedEndSequencesWithQuality]',manifest,'PairedEndFastqManifestPhred33V2')
			demux_PE.save(output_path+"demux.qza")
#			demux_PE = Artifact.load(dada2_path+"demux.qza")
			Time_sc(time.time(), start)
			
			print("\n\n################# DADA2 #################")
			print("DADA2 : paired-end-demux.qza --> table, representative-seqeuences, denoising-stats")
			start = time.time()
			dada2_demux = demux_PE
			#v34_ad_f = ["CCTACGGGNGGCWGCAG"]
			#v34_ad_r = ["GGACTACNVGGGTWTCTAAT"]
			#demux_cleanAT = cutadapt.methods.trim_paired(demux_PE, front_f = v34_ad_f, front_r = v34_ad_r, cores=cpus).trimmed_sequences
			#dada2_demux = demux_cleanAT
			#demux_cleanAT.save(dada2_path+"demux_cleanAT.qza")
			left_f = f_primer_len
			left_r = r_primer_len
			demux_PE_dada2 = dada2.methods.denoise_paired(demultiplexed_seqs = dada2_demux,trunc_len_f = trunc_f, trunc_len_r = trunc_r, trim_left_f = left_f, trim_left_r = left_r, max_ee_f = ee_f, max_ee_r = ee_r , n_threads = int(cpus))

			df_stats = demux_PE_dada2.denoising_stats.view(Metadata).to_dataframe()
			#df_stats.index.name = "sample-id"
			#df_stats.reset_index(inplace=True)
			df_fil_stats = stats_depth(df_stats, filtering_depth)
			df_fil_stats.to_csv(output_path+"stats_preprocessing.txt", sep="\t",index=True)
			

			demux_PE_dada2.table.save(output_path+"table.qza")
			demux_PE_dada2.representative_sequences.save(output_path+"representative_sequences.qza")
			print("\n\n################# ASV filtering #################\n")
			print("Step1. Sample filtering : sampling depth < "+str(filtering_depth)+"\n")
			q2_filtable = feature_table.methods.filter_samples(table =demux_PE_dada2.table, metadata = Metadata.load(output_path+"stats_preprocessing.txt"), where = 'analysis=="success"').filtered_table
			print("Step2. The minimum number of samples that a feature must be observed in to be retained : 1\n")
			q2_fil_table = feature_table.methods.filter_features(table= q2_filtable, min_samples =1).filtered_table
			q2_fil_table.save(output_path+"table_filt.qza")


			Time_sc(time.time(), start)
			os.system("rm -rf "+rename_path)

		except Exception as ex:
			print("\n\n################## ERROR ###################")
			print(ex)
			
	except Exception as ex:
		print("\n\n################## ERROR [V34-region] ###################")
		print(ex)
		print(">> python "+file_name+" Metagnome_outputPath f_primer_len r_primer_len trimming_len rare-cut-depth cpu_cores")
		print("V3V4 >> python "+file_name+" Metagnome_outputPath 17 21 450 3000 10")



if __name__=="__main__":
	main()
