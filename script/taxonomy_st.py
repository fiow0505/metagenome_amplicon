#2023-08-24; sindy (sindy@hunbiome.com)
#!/usr/bin/env python
import sys, os
from time import localtime, strftime
from optparse import OptionParser as OP
from pandas.api.types import CategoricalDtype
import pandas as pd
from itertools import combinations
import scipy.stats as st, numpy as np

debug = 0
def cmdparameter(argv):
    if len(argv) == 1:
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein", help="Taxonomy relative abundance table")
    parser.add_option("-m", "--map-file", dest="metadata", help="Map file containing group information")
    parser.add_option("-n", "--group", dest="group", help="The column contains group information.")
    parser.add_option("-s","--grpsep", dest="sort_value",help="group information.  ex.Group1,Group2,Group3")
    parser.add_option("-p","--paired", dest="paired",default=False, help="If TRUE, Wilcoxon signed rank test. FALSE - MannWitney U test=Wilcoxon rank sum test")
    parser.add_option("-o", "--output-path", dest="output_path", help = "output_path")
    parser.add_option("-d", "--debug", dest="debug", default=False, action="store_true", help="Debug the program")
    parser.add_option("-a", "--sample_append", dest="sample_append", default="False", help="append the composition table")

    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)

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

def write_excel(excel_path, df, name): # 엑셀 파일로 저장
    if not os.path.exists(excel_path):
        with pd.ExcelWriter(excel_path, mode='w', engine='openpyxl') as writer:
            df.to_excel(writer, index=True, sheet_name=name)
    else:
        with pd.ExcelWriter(excel_path, mode='a', engine='openpyxl') as writer:
            df.to_excel(writer, index=True, sheet_name=name)
    return 0


def main():
    global debug
    options, args = cmdparameter(sys.argv)
    debug = options.debug
    taxa_f = options.filein
    sample_f = options.metadata
    groupID = options.group
    sort_value = options.sort_value
    cobn_paired = options.paired
    output_d = options.output_path
    append_sample = options.sample_append
 
    metadata = pd.read_csv(sample_f, sep="\t", dtype=str)
    colname = groupID
    group_value = metadata[colname].unique().tolist() if sort_value==None else sort_value.split(",")
    df_meta = metadata[["sample-id",colname]].rename(columns = {colname:"group_type"})
    df_meta = df_meta[df_meta["group_type"].isin(group_value)]
    li_samples = df_meta["sample-id"].tolist()

    df_groupby_count = df_meta.groupby("group_type").count().rename(columns = {"sample-id":"N"}).reset_index()[["group_type", "N"]]
    cat_group_order = CategoricalDtype(group_value, ordered=True)
    df_meta["group_type"] = df_meta["group_type"].astype(cat_group_order)
    
    output_d = os.path.join("Taxonomy-RelativeAbundance", "_vs_".join(group_value)) if output_d==None else output_d
    if os.path.exists(output_d)==False:
        os.makedirs(output_d)
    
    if (any(df_groupby_count["N"]==2)==True):
	    s_method = "Mann-Whitney U test"# "T-test"
    else:
        s_method = "Wilcoxon signed rank test" if cobn_paired=="TRUE" else "Mann-Whitney U test"

    print("paired-test %s" %cobn_paired)
    print("statistical method %s" %s_method) 
    
    df_s_relev = pd.read_csv(taxa_f, sep="\t", index_col=0)
    li_samples = list(set(df_s_relev.columns.tolist())&set(li_samples))
    df_s_relev = df_s_relev[li_samples]
    if(df_s_relev[li_samples[0]].sum(axis=0)<=100):
        df_s_relev =  df_s_relev[li_samples]*100
    df_s_relev = df_s_relev[df_s_relev.sum(axis=1)>0]
    #df_s_relev.index = df_s_relev.index.str.replace(";__", ";__unclassification")
    df_s_relev.index.name ="Full_lineage"
    df_prevalence = pd.DataFrame({"prevalence(%)":(df_s_relev>0).sum(axis=1)*100/len(li_samples)})
    
    li_get_taxon = df_s_relev.index.tolist()
    #filter = [';' in taxon for taxon in li_get_taxon]
    #li_get_taxon_filter = [i for indx,i in enumerate(li_get_taxon) if filter[indx] == True]
    df_s_relev_melt = df_s_relev.T.reset_index().rename(columns = {"index":"sample-id"}).melt(id_vars="sample-id")
    df_s_relev_melt_meta = pd.merge(df_s_relev_melt, df_meta[["sample-id", "group_type"]].astype(str), on="sample-id", how = "inner")
    df_mean = df_s_relev_melt_meta.groupby("Full_lineage").mean().rename(columns = {"value":"All.MEAN"})#*100
    df_groupby_mean = df_s_relev_melt_meta.groupby(["Full_lineage","group_type"]).mean().reset_index().pivot(index = "Full_lineage", columns="group_type")
    df_groupby_mean.columns = df_groupby_mean.columns.droplevel()
    df_groupby_mean = df_groupby_mean[group_value]

    
    li_st = []
    li_kruskal_p =[]
    if (len(group_value)>2):
        p_colname = " vs ".join(group_value)
        for taxon in li_get_taxon:
            if (len(taxon.split(";"))==1):
                kruskal_p = ""
            else: 
                kruskal_p = st.kruskal(*[group["value"].tolist() for name, group in df_s_relev_melt_meta[df_s_relev_melt_meta["Full_lineage"]==taxon].groupby("group_type")]).pvalue 
            li_kruskal_p.append(kruskal_p)
            
        df_groupby_mean.loc[li_get_taxon,p_colname] = li_kruskal_p
        df_groupby_mean = df_groupby_mean[group_value+[p_colname]]
        li_st.append(pd.DataFrame({"comparison_set": [p_colname], "statistical_method": "kruskal-wallis"}))

    li_combi = list(combinations(group_value, 2))
    for combi_value in li_combi:
        p_colname = " vs ".join(combi_value)
        if s_method=="T-test":
            li_combi_p = [st.ttest_ind(*[group["value"].tolist()  for name, group in df_s_relev_melt_meta[(df_s_relev_melt_meta["Full_lineage"]==taxon)&((df_s_relev_melt_meta["group_type"]==combi_value[0])|(df_s_relev_melt_meta["group_type"]==combi_value[1]))].groupby("group_type")], alternative="two-sided").pvalue if df_s_relev_melt_meta[(df_s_relev_melt_meta["Full_lineage"]==taxon)&((df_s_relev_melt_meta["group_type"]==combi_value[0])|(df_s_relev_melt_meta["group_type"]==combi_value[1]))].sum().value!=0 else '' for taxon in li_get_taxon]
        elif s_method=="Mann-Whitney U test":
            li_combi_p = [st.mannwhitneyu(*[group["value"].tolist()  for name, group in df_s_relev_melt_meta[(df_s_relev_melt_meta["Full_lineage"]==taxon)&((df_s_relev_melt_meta["group_type"]==combi_value[0])|(df_s_relev_melt_meta["group_type"]==combi_value[1]))].groupby("group_type")], alternative="two-sided").pvalue if df_s_relev_melt_meta[(df_s_relev_melt_meta["Full_lineage"]==taxon)&((df_s_relev_melt_meta["group_type"]==combi_value[0])|(df_s_relev_melt_meta["group_type"]==combi_value[1]))].sum().value!=0 else '' for taxon in li_get_taxon]
        elif s_method =="Wilcoxon signed rank test":
            li_combi_p = [st.wilcoxon(*[group["value"].tolist()  for name, group in df_s_relev_melt_meta[(df_s_relev_melt_meta["Full_lineage"]==taxon)&((df_s_relev_melt_meta["group_type"]==combi_value[0])|(df_s_relev_melt_meta["group_type"]==combi_value[1]))].groupby("group_type")], alternative="two-sided").pvalue if (df_s_relev_melt_meta[(df_s_relev_melt_meta["Full_lineage"]==taxon)&((df_s_relev_melt_meta["group_type"]==combi_value[0])|(df_s_relev_melt_meta["group_type"]==combi_value[1]))].sum().value!=0)&(";" in taxon) else '' for taxon in li_get_taxon]
        elif s_method =="Kruskal-test":
            li_combi_p = [st.kruskal(*[group["value"].tolist()  for name, group in df_s_relev_melt_meta[(df_s_relev_melt_meta["Full_lineage"]==taxon)&((df_s_relev_melt_meta["group_type"]==combi_value[0])|(df_s_relev_melt_meta["group_type"]==combi_value[1]))].groupby("group_type")]).pvalue if (df_s_relev_melt_meta[(df_s_relev_melt_meta["Full_lineage"]==taxon)&((df_s_relev_melt_meta["group_type"]==combi_value[0])|(df_s_relev_melt_meta["group_type"]==combi_value[1]))].sum().value!=0)&(";" in taxon) else '' for taxon in li_get_taxon]
            
        df_groupby_mean.loc[li_get_taxon, p_colname] = li_combi_p
        li_st.append(pd.DataFrame({"comparison_set": [p_colname], "stats_method.two_group": s_method}))
    df_st_info = pd.concat(li_st)
    df_groupby_mean[group_value] = df_groupby_mean[group_value]#*100
    df_groupby_save = pd.merge(df_prevalence, df_groupby_mean, left_index=True, right_index=True)
    df_groupby_save = pd.merge(df_prevalence,df_mean, left_index=True, right_index=True).merge(df_groupby_mean, left_index=True, right_index=True)
    if(append_sample=="TRUE"):
        df_groupby_save = df_groupby_save.merge(df_s_relev, left_index=True, right_index=True)
    df_groupby_save.insert(0, "Taxon", rename_taxon(df_groupby_save.index.tolist(), ";"))
    df_groupby_save = df_groupby_save.reset_index()
    df_groupby_save.insert(2, "Taxon_rank", df_groupby_save["Full_lineage"].str.split(";").apply(lambda x: len(x)).tolist())

    dic_name = {2:"Phylum", 3:"Class", 4:"Order", 5:"Family",6:"Genus", 7:"Species"}
    excel_path = os.path.join(output_d, "Bacteria.statstical.xlsx")
    if os.path.exists(excel_path):
        os.remove(excel_path)
     
    write_excel(excel_path, df_groupby_count, "number_of_samples")
    for lev, df in df_groupby_save.groupby("Taxon_rank"):
        if lev==1:
            pass
        elif lev<=7:
            if (len(group_value)>2):
                df['stats_method.all_group'] = 'Kruskal-wallis'
            df["stats_method.two_group"] = s_method
            write_excel(excel_path, df.sort_values("All.MEAN", ascending=False).reset_index(drop=True), dic_name[lev])
        


    
    
if __name__ == '__main__':
    main()

