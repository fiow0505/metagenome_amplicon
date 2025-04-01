#-*- coding:utf-8 -*-
#!/bin/bash
# update : 2022.11
import sys, os
import pandas as pd
from optparse import OptionParser as OP
from qiime2.plugins import (diversity, feature_table, feature_classifier, taxa)
from qiime2 import (Artifact, Metadata)
import time
import datetime
import re



debug = 0
def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-table", dest="filein", help="ASVs count table")
    parser.add_option("-c", "--classification-table", dest="classification", help="classification table")
    parser.add_option("-t", "--taxon-filtering", dest = "taxon_filt", default = True, 
                      help = "taxonomy filtering <Default = TRUE>") 
    parser.add_option("--include-taxon", dest="include_taxon", default = "d__Bacteria", 
                      help="include taxon <Default = 'Bacteria'>", )
    parser.add_option("--exclude-taxon", dest="exclude_taxon", default = "chloroplast,mitochondria", 
                      help="exclude taxon <Default = 'chloroplast,mitochondria'>")
    parser.add_option("-o", "--output-path", dest="output_path", help = "output_path", default='taxonomy_composition')
    parser.add_option("-d", "--debug", dest="debug", default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)



def Taxonomy_level(q2_fil_table, q2_taxonomy, output_d):
    li_lev_relab = list()
    li_lev_abun = list()
    for lev in range(1,8): #lev 1: Kingdom > .... lev 7: Species
        print("Level is "+str(lev))
        q2_collapse_table = taxa.methods.collapse(table = q2_fil_table, taxonomy = q2_taxonomy, level=lev).collapsed_table
        df_abun_table = q2_collapse_table.view(pd.DataFrame).T
        df_abun_table.index.name = "Taxon"
        df_abun_table.reset_index(inplace=True)
        df_abun_table["Taxon"] = df_abun_table["Taxon"].str.replace(";__", ";__unclassification").tolist()
        df_abun_table.to_csv(os.path.join(output_d, "taxonomy-level"+str(lev)+".Abundance.tsv"), sep="\t", index=False)
        li_lev_abun.append(df_abun_table)
        #샘플마다 미생물 측정 합이 "1"
        q2_re_table = (feature_table.methods.relative_frequency(table = q2_collapse_table)).relative_frequency_table 
        df_re_table = q2_re_table.view(pd.DataFrame).T
        df_re_table.index.name = "Taxon"
        df_re_table.reset_index(inplace = True)
        df_re_table["Taxon"] = df_re_table["Taxon"].str.replace(";__", ";__unclassification").tolist()
        df_re_table.to_csv(os.path.join(output_d, "taxonomy-level"+str(lev)+".RelAbundance.tsv"), sep="\t", index=False)
        li_lev_relab.append(df_re_table)
    pd.concat(li_lev_abun).to_csv(os.path.join(output_d, "taxonomy-levelAll.Abundance.tsv"), sep="\t", index=False)
    pd.concat(li_lev_relab).to_csv(os.path.join(output_d, "taxonomy-levelAll.RelAbundance.tsv"), sep="\t", index=False)
    return 


def main():
    global debug
    options, args = cmdparameter(sys.argv)
    debug = options.debug
    tblname = os.path.abspath(options.filein)
    clname = os.path.abspath(options.classification)
    output_d = options.output_path
    
    filtering_op = options.taxon_filt
    include_taxon = options.include_taxon
    exclude_taxon = options.exclude_taxon

    q2_table = Artifact.load(tblname)
    q2_taxonomy = Artifact.load(clname)
    if os.path.exists(output_d)==False:
        os.makedirs(output_d)
    
    if filtering_op==True:
        q2_fil_table = taxa.methods.filter_table(q2_table, q2_taxonomy, include=include_taxon, exclude = exclude_taxon).filtered_table
        
    else:
        q2_fil_table = q2_table
    q2_fil_table.save(os.path.join(output_d, 'taxonomy.table.qza'))
    df_taxonomy_tbl = pd.DataFrame({"taxonomy.filtered":q2_fil_table.view(pd.DataFrame).sum(axis=1)}).reset_index().rename(columns = {"index":"sample-id"})
    df_taxonomy_tbl.to_csv(os.path.join(output_d, "taxonomy.table.stats.txt"), sep="\t", index=False)
	

    with open(os.path.join(output_d, "taxonomy_abundance_info.txt"), "w") as f:
        f.write("Running Date: "+datetime.date.today().strftime('%Y.%m.%d')+"\n")
        f.write("ASVs table: "+tblname+"\n")
        f.write("classification table: "+clname+"\n")
        f.write("\t> v34 amplicon seq; taxonomy.sklearn - SILVA 138v")
        f.write("\t> full-length seq; taxonomy.vsearch - GTDB 207r")
        f.write("Taxonomy filtering: "+str(filtering_op)+"\n")
        if filtering_op==True:
            f.write("\t> Include taxon: "+include_taxon+"\n")
            f.write("\t> Exclude taxon: "+exclude_taxon+"\n")
        
    Taxonomy_level(q2_fil_table, q2_taxonomy, output_d)
    



if __name__ == "__main__":
    main()
