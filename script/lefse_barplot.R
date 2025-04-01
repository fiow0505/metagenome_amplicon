#!/usr/bin/env Rscript
# Author : Kyeongeui Yun (sindy@hunbiome.com)
# 2023-07-17 sindy@hunbiome.conm
options(warn = -1) # Turn off warning

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
pkgs = c("optparse", "ggplot2", "tidyverse", "ggprism")
for(p in pkgs) suppressWarnings(suppressMessages(library(p, quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE)))

if (TRUE){
  option_list = list(
    make_option(c("-i", "--input_file"), type="character", default="",
                help="LEfSe result directory"),
    make_option(c("-s", "--sort_value"), type="character", default="",
                help="sort by group value.  ex.Group1,Group2,Group3"),
    make_option(c("-r", "--rank"), type="character", default="all",
                help="Taxonomy rank [default %default]"),
    make_option(c("-c", "--cutoff_lda"), type="numeric", default=2,
                help="LDA cut-off value [default %default]"),
    make_option(c("-E", "--E_format"), type="character", default="png",
               help="The type of output figures. select: png, pdf [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(ggprism)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
rename_taxon <- function(li_taxon, sep){
    ch_taxon = as.character()
    for (taxon in li_taxon){
        li_sp_taxon = unlist(str_split(taxon, sep))
        lev_taxon =length(li_sp_taxon)
        st_num = ifelse(any(grepl(paste(c("_UCG", "__un", "__metagenome", "gut"),collapse="|"), li_sp_taxon))==TRUE, 
                        grep(paste(c("_UCG", "__un", "__metagenome", "gut"),collapse="|"), li_sp_taxon)[1]-1, lev_taxon)
        
        new_taxon = ifelse(st_num==lev_taxon,li_sp_taxon[lev_taxon], paste0(li_sp_taxon[c(st_num,lev_taxon)], collapse = sep))
        
        ch_taxon = c(ch_taxon, new_taxon)
    }
    if(any(duplicated(ch_taxon))==TRUE){
        duplicated_taxon = ch_taxon[duplicated(ch_taxon)]
        duplicated_i = grep(paste0(duplicated_taxon, collapse = "|"), ch_taxon)
        for(i in duplicated_i){
            li_sp_taxon = unlist(str_split(li_taxon[i], sep))
            ch_taxon[i] = paste0(c(li_sp_taxon[length(li_sp_taxon)-1], ch_taxon[i]), collapse = sep)
        }
    }
    return(ch_taxon)}

rename_lefse <- function(li_get_taxon){
    ret = li_get_taxon
    pattern = paste0(c(' ','$','@','#','%','^','&','\\*','\"','\''), collapse = "|")
    ret2 = str_replace_all(ret,pattern,'_')
    pattern = paste0(c("/",'\\(','\\)','\\-', '\\+','=','\\{','\\}','\\[','\\]',',',
                      '\\.',':','\\?','\\<','\\>','\\,'), collapse = "|")
    ret3 = str_replace_all(ret2,pattern,'_')
    ret4 = str_replace_all(ret3,"\\|",';')
    return(ret4)
}

li_levelName = list("2"="Phylum", "3"="Class", "4"="Order", "5"="Family","6"="Genus","7"="Species", "all" = "All_rank")
li_colorvalues = brewer.pal(n=5,"Set1")
#c('#4daf4a','#e41a1c','#377eb8','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf')

inputfile = opts$input_file
sort_value = opts$sort_value
cutLDA = opts$cutoff_lda
cutlevel = opts$rank
levelName = as.character(li_levelName[cutlevel])
output_path = dirname(normalizePath(inputfile))

print(paste("Loading) LEfSe output formatted-file:", inputfile))
df_lefse_out = read.table(inputfile, sep = "\t", header = T, check.names = F) %>% na.omit()  %>% mutate_at("full_lineage", str_replace_all, fixed("."), ";") %>% filter(!grepl("__unclassification", .$full_lineage))
if(any(grepl("Taxon", colnames(df_lefse_out)))==FALSE){df_lefse_out["Taxon"] = rename_taxon(df_lefse_out$full_lineage, ";")}
if(any(grepl("level", colnames(df_lefse_out)))==FALSE){df_lefse_out["level"] = unlist(lapply(df_lefse_out$full_lineage, function(x){length(unlist(str_split(x, fixed(";"))))}))}

if(sort_value!=""){
    group_value = unlist(str_split(sort_value, ","))}else{
    group_value = sort(unique(df_lefse_out$high_group))}

df_lefse_cutLDA = df_lefse_out %>% filter(.$LDA_score>= cutLDA)
if(cutlevel!="all"){
    df_lefse_cutLDA_level=df_lefse_cutLDA %>% filter(.$level==as.numeric(cutlevel)) }else{
    df_lefse_cutLDA_level=df_lefse_cutLDA}

df_lefse_cutLDA_level$high_group = factor(df_lefse_cutLDA_level$high_group, levels = group_value)
df_lefse_cutLDA_level = df_lefse_cutLDA_level %>% dplyr::arrange(high_group, desc(LDA_score), level)

df_lefse_cutLDA_level$Taxon = factor(df_lefse_cutLDA_level$Taxon, levels = rev(df_lefse_cutLDA_level$Taxon))
nfeatures = df_lefse_cutLDA_level %>% nrow()
if(nfeatures!=0){
    plot_height = 1+(0.1*nfeatures) #ifelse(nfeatures<5, 1+0.1*nfeatures,  (nfeatures/5)*1.2)
    plot_colors = li_colorvalues[1:length(group_value)]
    names(plot_colors) = group_value
    
    outputprefix = paste0("LEfSe_formatted.LDAupper", cutLDA, ".",  levelName)
    setwd(output_path)
    print(paste("Making) LEfSe LDA plot:", outputprefix))
    write.table(df_lefse_cutLDA_level, paste0(outputprefix, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F)
    g = ggplot(df_lefse_cutLDA_level, aes(x= Taxon, y = LDA_score, fill = high_group))+
        geom_bar(stat = "identity", show.legend=T)+
        geom_hline(yintercept = cutLDA, linetype = "dashed")+
        theme_prism(base_family = 'serif', base_fontface = 'plain', base_size = 6)+
        theme(plot.title = element_text(hjust = 0), axis.text.y = element_text(family = 'sans',face="italic"),
              legend.position = "top", legend.justification='left',
              legend.title = element_text(size=10), legend.key.size = unit(2, 'mm'))+
        labs(x="", y = "LDA score (log10)", fill = "", title = paste0("Features in ", levelName))+
        coord_flip()+
        scale_fill_manual(values = plot_colors, drop=F)
    ggsave(paste0(outputprefix, ".", opts$E_format), width = 7, height = plot_height, g, limitsize = FALSE)

}   

