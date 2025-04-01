#2023-08-23
# Author : Kyeongeui Yun (sindy@hunbiome.com)
#source /data/shared/enviroment/add.sh
#!/usr/bin/env Rscript

options(warn = -1)
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character",
                help="result of beta distance matrix (ex. braycurtis_dm.txt)"),
    make_option(c("-m", "--metadata"), type="character",
                help="Metadata"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-s", "--sort_value"), type="character", default="",
                help="sort by group value.  ex.Group1,Group2,Group3"),
    make_option(c("-o", "--output_path"), type="character", default="",
                help="output_path;[default %default]"),
    make_option(c("-t", "--title"), type="character", default="",
                help="title name [default %default]"),
    make_option(c("-k", "--korean"), type="logical", default=FALSE,
                help="Korean language [default %default]"),
    make_option(c("-f", "--fontsize"), type="numeric", default = 8,
               help="font size [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=7,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=7,
                help="Figure heidth in mm [default %default]"),
    make_option(c("-E", "--E_format"), type="character", default="png",
               help="The type of output figures. select: png, pdf [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(vegan)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggprism)))
suppressWarnings(suppressMessages(library(RColorBrewer)))

dmfile = tools::file_path_as_absolute(opts$input)
metafile = tools::file_path_as_absolute(opts$metadata)
groupID = opts$group
sort_value = opts$sort_value
output_path = ifelse(opts$output_path=="", str_replace_all(sort_value, ",", "_vs_"), opts$output_path)
language_font = ifelse(opts$korean, "AppleGothic", "serif")
dmfile_name = gsub(".txt", "", basename(dmfile))
title = ifelse(opts$title=="", ifelse(dmfile_name=="braycurtis_dm", "Bray-Curtis", "unweighted UniFrac"), opts$title)

pcoaname = paste0(dmfile_name, ".pcoa")
stname = paste0(dmfile_name, "-st.txt")
output_pcoa = paste0(output_path,"/", pcoaname)
output_st = paste0(output_path,"/", stname)


if(dir.exists(output_path)==FALSE){dir.create(output_path)}
str.pval<- function (p.value) {
    unclass(symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("p<0.001", "p<0.01", "p<0.05", "p<0.1", "p>0.1")))
}
metadata = read.table(metafile, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F, check.names = FALSE, fill=TRUE)[groupID] %>% na.omit()
colnames(metadata) = "group"
if(sort_value!=""){group_value = unlist(str_split(sort_value, ","))}else{group_value = unique(metadata$group)}

filt_row = metadata$group %in% group_value
df_meta = data.frame("row.names" = rownames(metadata)[filt_row], "group" = metadata$group[filt_row])
samples = rownames(df_meta)
df.dist = read.table(dmfile, sep = "\t",header = T, check.names = F, row.names = 1)


set.seed(42)
if(length(setdiff(samples, rownames(df.dist)))!=0){
    print(paste0("Fail samples:", paste0(setdiff(samples, rownames(df.dist)), collapse=", "))) 
    samples = intersect(rownames(df.dist),samples)
}
df.dist = df.dist[samples,samples] %>% as.matrix(labels=T)
dm = df.dist %>% as.dist()
dm_adonis = adonis2(dm ~ metadata[labels(dm),], permutations=999, parallel = 4)[1,]
rownames(dm_adonis) = "distance_matrix"
rownames(dm_adonis) = paste0(group_value, collapse = " vs ")
if(length(group_value)>2){
    li_combn_adonis = list()
    li_combn = combn(group_value,2, simplify = F)
    for(combn_i in 1:length(li_combn)){
        combn_value = li_combn[combn_i] %>% unlist()
        combn_samples = intersect(samples, df_meta %>% filter(group %in% combn_value) %>% rownames())
        combn_dm = df.dist[combn_samples,combn_samples] %>% as.dist()
        combn_dm_adonis = adonis2(combn_dm ~ metadata[labels(combn_dm),], permutations=999, parallel = 4)[1,]
        rownames(combn_dm_adonis) = paste0(combn_value, collapse = " vs ")
        li_combn_adonis[[combn_i]] = combn_dm_adonis
    }
    df_adonis = rbind(dm_adonis, Reduce("rbind", li_combn_adonis))
}else{df_adonis = dm_adonis}

caption = paste0(str.pval(df_adonis$`Pr(>F)`)[1], ", PERMANOVA-test")
res_pcoa = ape::pcoa(df.dist)
meta_pcoa = res_pcoa$vectors[,1:3] %>% as.data.frame()  %>% merge(df_meta, by="row.names")
group_mean = meta_pcoa[c("Axis.1", "Axis.2", "Axis.3", "group")] %>% dplyr::group_by(group) %>% dplyr::summarise(`mean.Axis.1` = mean(`Axis.1`), `mean.Axis.2` = mean(`Axis.2`), `mean.Axis.3` = mean(`Axis.3`))

df_plot = meta_pcoa %>% merge(group_mean, by = "group")
df_plot$group = factor(df_plot$group, levels = group_value)
df_plot %>% write.table(paste0(output_pcoa, ".txt"), sep="\t", quote=F, row.names = F, col.names=T)

xlabel = paste0("PCoA 1 : ",round(100*res_pcoa$values$Relative_eig[1],digits = 1),"%")
ylabel = paste0("PCoA 2 : ",round(100*res_pcoa$values$Relative_eig[2],digits = 1),"%")



li_color = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf', '#66c2a5','#fc8d62', '#8da0cb','#a6d854','#ffd92f','#e5c494','#b3b3b3')[1:length(group_value)]
g <- ggplot(df_plot, aes(x = `Axis.1`, y=`Axis.2`, group=group))+geom_point(aes(color = group))+
stat_ellipse(aes(color=group, fill=group), type = 'norm', level=0.9)+
geom_point(mapping = aes(x=`mean.Axis.1`, y=`mean.Axis.2`), shape=10, cex=0.75)+
geom_segment(aes(x = `mean.Axis.1`, y=`mean.Axis.2`, xend = `Axis.1`, yend= `Axis.2`, color = group), lwd=.1)+
theme_prism(base_family = 'serif', base_fontface = 'plain', base_size = opts$fontsize)+
theme(plot.title = element_text(hjust = 0.5), 
      legend.title = element_blank(), legend.title.align=0.5, 
      legend.key.height= unit(.5, 'cm'), legend.text = element_text(size=opts$fontsize))+
labs(title = title, x= xlabel, y = ylabel, caption = caption)+scale_color_manual(values=li_color)+scale_fill_manual(values = li_color)
ggsave(paste0(output_pcoa, ".", opts$E_format), width = opts$width, height = opts$height, g)
df_adonis %>% rownames_to_column("comparison_set") %>% write.table(output_st, sep="\t", quote = F, row.names = F, col.names = T)


#df_adonis %>% rownames_to_column("comparison_set") %>% write.table(paste0(gsub(".pcoa.png","",output_pcoa), ".twogroup-st.txt"), sep="\t", quote = F, row.names = F, col.names = T)


