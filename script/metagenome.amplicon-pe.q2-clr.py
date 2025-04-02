#!/bin/bash
# Author : Kyeongeui Yun (fiow3250@gmail.com)
# update : 2022.11
import os
import pandas as pd
import sys
from qiime2.plugins import (dada2, diversity, feature_table, feature_classifier, taxa)
from qiime2 import (Artifact, Metadata)
import time
import datetime
import re

db_path=os.path.join(os.path.dirname(os.path.abspath(__file__)), "database/")
def DirPath(path_loc):
	dir_path = path_loc
	if(dir_path[-1]!="/"):
		dir_path = dir_path+"/"
	return dir_path

def Taxonomy(dada2_path,db_name,cpus, con):
	dada2_rep_seq = Artifact.load(dada2_path+"representative_sequences.qza")
	db_classifier = Artifact.load(os.path.join(db_path,db_name))
	# Sklearn 을 이용하여 Taxonomy classification 진행
	taxonomy_out = (feature_classifier.methods.classify_sklearn(reads = dada2_rep_seq, classifier = db_classifier, n_jobs=int(cpus), confidence = con)).classification
	return taxonomy_out

def main():
	try:
		dada2_path=DirPath(sys.argv[1])
		os.system("export TMPDIR="+dada2_path) #주의 : 서버 공간 부족시
		taxonomy_path = dada2_path
		region = sys.argv[2].upper() # FULL mapping or V34 mapping or V4 mapping
		cpus = int(sys.argv[3])
		q2_v = (sys.argv[4])
		con = 0.7 #default
			
		print("################## Taxonomy Classifier ##################")
		print("Metagenome 16s rRNA output path : ", dada2_path) #분석 후 저장될 위치
		print("Taxonomy save path : ", taxonomy_path) #분석후 저장될 위치랑 같음
		print("Region : ", region)
		print("cpus : ", cpus) # Maximum 40
		print("Confidence cut-off : ",con) # default 0.7로 지정되어있음

		if region=="34":
			#qiime2-2020.11 : V34-region에 맞춰 유전자 서열 Extract 한 DB
			#HanLab V3-V4 region : 341F(CCTACGGGNGGCWGCAG), 806RB(GGACTACNVGGGTWTCTAAT) /보통 341F,805R과 다름
			#https://docs.qiime2.org/2020.11/tutorials/feature-classifier/ 참조 (Extract reference reads)
			#gg_db_name = "q2-2020.11-gg-13-8-99-341-806-nb-classifier.qza"
			#s_db_name = "q2-2020.11-silva-138-99-341-806-nb-classifier.qza"
			gg_db_name = "gg-13-8-99-341-805-nb-classifier.qza"
			s_db_name = "silva-138-99-341-805-nb-classifier.qza"
	
		elif region=="4":
			#qiime2-2020.11 : V4-region에 맞춰 유전자 서열 Extract 한 DB
			gg_db_name = "gg-13-8-99-515-806-nb-classifier.qza"
			s_db_name = "silva-138-99-515-806-nb-classifier.qza"
		
		elif region=="FULL":
			gg_db_name = "gg-13-8-99-nb-classifier.qza"
			s_db_name = "silva-138-99-nb-classifier.qza"
		elif region=="34_BF":
			gg_db_name = "gg-13-8-99-341-806-nb-classifier.qza"
			s_db_name = "silva-138-99-341-806-nb-classifier.qza"
		gg_db_name = q2_v+"-"+gg_db_name
		s_db_name = q2_v+"-"+s_db_name
		print("Taxonomy classification DB")
		print("GreenGene 13.8v : ", gg_db_name)
		print("SILVA 138v : ", s_db_name)
		#gg_taxonomy = Taxonomy(dada2_path, gg_db_name, cpus, con)
		#gg_taxonomy.save(taxonomy_path+"q2-2020.11-GG_13_8.qza")
		print(s_db_name)
		s_taxonomy = Taxonomy(dada2_path, s_db_name, cpus, con)
		s_taxonomy.save(os.path.join(taxonomy_path, "taxonomy.sklearn.qza"))

	except Exception as ex:
		print(">> python default.q2-pipeline_2.py Metagenome_outputPath mapping_region(FULL or 34 or 4) cpu qiime2_v(q2-2020.11 or q2-2022.8)")
		print("Error >>>>>>>> ")
		print(ex)

if __name__ == "__main__":
	main()
