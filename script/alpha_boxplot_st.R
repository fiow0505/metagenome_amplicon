# Author : Kyeongeui Yun (fiow3250@gmail.com)
#!/usr/bin/env Rscript
# 2023.09.20
options(warn = -1) # Turn off warning

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/alpha-diversity.txt",
                help="Alpha diversity matrix [default %default]"),
    make_option(c("-m", "--metadata"), type="character", default="result/metadata.txt",
                help="Metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-s", "--sort_value"), type="character", default="",
                help="sort by group value.  ex.Group1,Group2,Group3"),
    make_option(c("-p", "--paired"), type="logical", default=FALSE,
                help = "if TRUE, Wilcoxon signed rank test [default %default]"),
    make_option(c("-t", "--transpose"), type="logical", default=FALSE,
                help="Design file or metadata [default %default]"),
    make_option(c("-a", "--alpha_index"), type="character", default="",
                help="Alpha metrics [default All metrics]"),
    make_option(c("-o", "--output_path"), type="character", default="",
                help="output_path"),
    make_option(c("-x", "--xlabAngle"), type="logical", default=TRUE,
                help="X lab set in angle [default %default]"),
    make_option(c("-k", "--korean"), type="logical", default=FALSE,
                help="Korean language [default %default]"),
    make_option(c("-f", "--fontsize"), type="numeric", default = 8,
               help="font size [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default = 10,
               help="figure width [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=5,
               help="figure height [default %default]"),
    make_option(c("-E", "--E_format"), type="character", default="png",
               help="The type of output figures. select: png, pdf [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

#suppressWarnings(suppressMessages(library(amplicon)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(ggprism)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(rstatix)))
alphafile = opts$input
metafile = opts$metadata
groupID = opts$group
sort_value = opts$sort_value
paired = opts$paired
output_path = ifelse(opts$output_path=="", str_replace_all(sort_value, ",", "_vs_"), opts$output_path)
if(dir.exists(output_path)==FALSE){dir.create(output_path)}
language_font = ifelse(opts$korean, "AppleGothic", "serif")

alpha_div = read.table(alphafile, row.names = 1, header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
if (opts$transpose){
  alpha_div = as.data.frame(t(alpha_div))
}
metadata = read.table(metafile, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F, check.names = FALSE)[groupID]
colnames(metadata) = "group"
df = merge(alpha_div, metadata, by = 'row.names', all=F)
if(opts$alpha_index==""){li_index = colnames(alpha_div)}else{li_index = opts$alpha_index}
df_plot_melt = df %>% melt(id.vars = c("Row.names", "group")) %>% filter(variable %in% li_index) %>% dplyr::rename("metrics" = "variable")
if(sort_value!=""){group_value = unlist(str_split(sort_value, ","))}else{group_value = unique(metadata$group)}
two_group = combn(group_value, 2, simplify = F)
df_plot_in = df_plot_melt %>% filter(group %in% group_value) %>% mutate("group" = factor(.$group, levels = group_value))  %>% mutate("metrics" = factor(.$metrics, levels = li_index))

p = ggplot(data =df_plot_in, aes(x = group, y=value, color = group))+
	geom_boxplot()+geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+
	facet_wrap(~ metrics, scale="free", nrow = 1)+labs(x = "", y="alpha-diversity")+
    #stat_compare_means(method = "kruskal", size=1.5, comparisons = two_group, aes(label = paste0("p=",..p.format.. )),symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1, 1.1),symbols = c("***", "**", "*", "p>0.05", "p>0.05")))
	stat_compare_means(method = "wilcox.test", size=1.5, comparisons = two_group, paired = paired, aes(label = paste0("p=",..p.format.. )))
#,symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1, 1.1),symbols = c("***", "**", "*", "p>0.05", "p>0.05")))


df_st = df_plot_in %>% group_by(metrics) %>% wilcox_test(value~group, paired = paired) %>% mutate("comparison" = paste0(.$group1, " vs ", .$group2)) %>% dcast(metrics~comparison, value.var = "p")
#df_st = df_plot_in %>% group_by(metrics) %>% kruskal_test(value~group) %>% mutate("comparison" = paste0(group_value, collapse = " vs ")) %>% dcast(metrics~comparison, value.var = "p")

if(length(group_value)>=3){
	#p = p+stat_compare_means(method = "kruskal.test", size=1.5, aes(label = paste0("kruskal,p=",..p.format.. )))
	df_kw = df_plot_in %>% group_by(metrics) %>% kruskal_test(value~group) %>% mutate("comparison" = paste0(group_value, collapse = " vs ")) %>% dcast(metrics~comparison, value.var = "p")
	df_st = df_kw %>% merge(df_st, by="metrics")
	}

p = p+theme_prism()+theme(text=element_text(family=language_font, size=opts$fontsize), legend.position = "top")+scale_color_brewer(palette = 'Set1')
if (opts$xlabAngle){
	p = p + theme(axis.text.x=element_text(angle=90,vjust=1, hjust=1))
}
ggsave(paste0(output_path, "/alpha-boxplot_st.", opts$E_format), width = opts$width, height = opts$height, p)

df_summary = plyr::ddply(df_plot_in, c("metrics", "group"), summarise, v = paste0(round(mean(value, na.rm=T),2),"±", round(sd(value, na.rm=T), 2))) %>% dcast(metrics~ group, value.var = "v")
df_summary_st = df_summary %>% merge(df_st, by="metrics")
df_summary_st$metrics = factor(df_summary_st$metrics, levels = li_index)
write.table(df_summary_st %>% dplyr::arrange(metrics) , paste0(output_path, "/alpha-boxplot_st.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
