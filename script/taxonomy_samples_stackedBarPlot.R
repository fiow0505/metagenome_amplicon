#!/usr/bin/env Rscript
# Author : Kyeongeui Yun (fiow3250@gmail.com)
# 2024. 10. 16

options(warn = -1) # Turn off warning

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
if (TRUE){
  option_list = list(
      make_option(c("-i", "--input"), type="character", 
                  help="Taxonomy composition table [silva138_bacteria-no-mitochondria-no-chloroplast.all.tsv]"),
      make_option(c("-m", "--metadata"), type="character", help="Metadata file"),
      make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
      make_option(c("-s", "--sort_value"), type="character", default="",
                help="sort by group value.  ex.Group1,Group2,Group3"),
      make_option(c("-o", "--output_path"), type="character", default="", help="Output path"),
      make_option(c("-l", "--level"), type="character", default = "All", help = "taxonomy level: Kingdom(1) ~ Species(7) [default. %default]"),
      make_option(c("-p", "--prevalence_min"), type="numeric", default = 30, help = "Prevalence(%) [default %default]"),
      make_option(c("-r", "--relev_min"), type="numeric", default = 0.001, help = "min relative abundance(%) [default  %default]"),
      make_option(c("-f", "--fontsize"), type="numeric", default = 7,  help = "png fontsize [default %default]"),
      make_option(c("-w", "--width"), type="numeric", default = 20,  help = "png height [default %default]"),
      make_option(c("-e", "--height"), type="numeric", default=5, help = "png height [default %default]"),
      make_option(c("--xlabel"), help = "xlabel(=sample-id) print [default %default]", type="character", default="TRUE")
      )
    opts = parse_args(OptionParser(option_list=option_list))
}


rename_taxon <- function(li_taxon, sep, op){
    ch_taxon = as.character()
    for (taxon in li_taxon){
        li_sp_taxon = unlist(str_split(taxon, sep))
        lev_taxon =length(li_sp_taxon)
        st_num = ifelse(any(grepl(paste(c("_UCG", "__un", "__metagenome", "gut"),collapse="|"), li_sp_taxon))==TRUE, 
                        grep(paste(c("_UCG", "__un", "__metagenome", "gut"),collapse="|"), li_sp_taxon)[1]-1, lev_taxon)
        
        new_taxon = ifelse(st_num==lev_taxon,li_sp_taxon[lev_taxon], paste0(li_sp_taxon[c(st_num,lev_taxon)], collapse = sep))
        
        ch_taxon = c(ch_taxon, new_taxon)
    }
    if(any(duplicated(ch_taxon))==op){
        duplicated_taxon = ch_taxon[duplicated(ch_taxon)]
        duplicated_i = grep(paste0(duplicated_taxon, collapse = "|"), ch_taxon)
        for(i in duplicated_i){
            li_sp_taxon = unlist(str_split(li_taxon[i], sep))
            ch_taxon[i] = paste0(c(li_sp_taxon[length(li_sp_taxon)-1], ch_taxon[i]), collapse = sep)
        }
    }
    return(ch_taxon)
}
save_data <- function(plt_excel, df, name) {
    library(openxlsx)
    wb <- createWorkbook(plt_excel)
    sheet <- creatSheet(wb, name)
    writeDataTable(wb, name, df, tableStyle = "TableStyleLight9")
    saveWorkbook(wb, plt_excel)
}

  
prev_fun <- function(x){sum(x>0)*100/length(x)}
li_rank = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colFunc <- colorRampPalette(c("purple4", "lightpink","firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))

library(tidyverse)
library(ggplot2)


inputfile = tools::file_path_as_absolute(opts$input)
metafile = tools::file_path_as_absolute(opts$metadata)
groupID = opts$group
sort_value = opts$sort_value
level = opts$level
if(level=="All"){rank = li_rank}else{rank = li_rank[as.numeric(level)]}

# Core-microbiome option
mean_value = opts$relev_min
prev_value = opts$prevalence_min

# Plot option
xlabel_out=opts$xlabel
plt_caption  = paste0("Core taxon: Prevalence>",prev_value,"%, Average>", mean_value,"%")
fontsize = opts$fontsize
plt_width = opts$width
plt_height = opts$height
output_path = opts$output_path
outputfile_name = paste0("samples.core-taxon-", rank, ".png")
pngfile = ifelse(output_path=="", outputfile_name, paste0(c(output_path, outputfile_name), collapse = "/"))


#df_relev =  read.table(inputfile, sep = "\t", header = T, check.names = F)
df_relev =  read.table(inputfile, sep = "\t", header = T, check.names = F)
df_relev = df_relev[lapply(str_split(df_relev[,1], ";"), function(x){length(x)})==level,]
rownames(df_relev) = df_relev[,1]
if(any(grepl(";__unclassification", rownames(df_relev)))==FALSE){
    rownames(df_relev) = str_replace_all(df_relev[,1], ";__", ";__unclassification")
}
df_relev = df_relev[,-1]
rownames(df_relev) = rename_taxon(rownames(df_relev), ";", TRUE)

meta_df = read.table(metafile, sep = "\t", header = T, check.names = F)
rownames(meta_df) = meta_df[,1]
meta_df = meta_df[groupID]
colnames(meta_df) = "group"

if(sort_value!=""){group_value = unlist(str_split(sort_value, ","))}else{group_value = unique(meta_df$group)}
get_meta_df = meta_df %>% filter(group %in% group_value)
get_samples = rownames(get_meta_df)
df_get_relev = df_relev[colnames(df_relev) %in% get_samples]*100
df_get_relev_melt = df_get_relev %>% rownames_to_column("Taxon") %>% reshape2::melt(id.vars = "Taxon") %>%
dplyr::rename("sample-id" = "variable") %>% merge(get_meta_df %>% rownames_to_column("sample-id"), by = "sample-id") 

group_filt_summary = df_get_relev_melt %>% group_by(Taxon, group) %>% dplyr::summarise(group_mean=mean(value), group_prev = prev_fun(value)) %>% filter((group_prev > prev_value)&(group_mean>mean_value))
core_taxon = group_filt_summary %>% dplyr::select(c("Taxon", "group")) %>% group_by(Taxon) %>% count() %>% filter(n==length(group_value)) %>% dplyr::select("Taxon") %>% unlist() %>% as.character()
taxon_order = group_filt_summary %>% filter(Taxon %in% core_taxon) %>% reshape2::dcast(Taxon~group, value.var = "group_mean") %>%
column_to_rownames("Taxon") %>% data.frame("mean" = apply(., 1, mean)) %>% arrange(desc(mean)) %>% rownames()
all_order = taxon_order
print(length(core_taxon))
df_core_melt = df_get_relev_melt %>% filter(Taxon %in% core_taxon)
df_other_melt = df_core_melt %>% group_by(`sample-id`) %>% dplyr::summarise("Taxon" = "other", value = 100-sum(value)) %>% merge(get_meta_df %>% rownames_to_column("sample-id"), by="sample-id") 
plt_input = rbind(df_core_melt, df_other_melt)

color_list <- colFunc(length(core_taxon))
taxon_color <- rev(c("black", color_list))
names(taxon_color) = c(core_taxon, "other")

df_relev_t = plt_input %>% reshape2::dcast(`sample-id`~Taxon, value.var = "value") %>% column_to_rownames("sample-id")
df_relev_t = df_relev_t[all_order]
colnames(df_relev_t) = paste0("col", 1:(length(all_order)))

#df_relev_t = df_relev_t[c(all_order, "other")]
#colnames(df_relev_t) = paste0("col", 1:(length(all_order)+1))
li_order_samples = as.character()
for(group_v in group_value){
    get_samples = (meta_df %>% filter(group==group_v)) %>% rownames()
    get_relev_t = df_relev_t[get_samples,]
    #get_relev_t = get_relev_t[colnames(get_relev_t)[rev(order(apply(get_relev_t, 2, mean)))]]
    order_get_samples = get_relev_t[rev(do.call(order, as.list(get_relev_t))),] %>% rownames()
    li_order_samples = c(li_order_samples, order_get_samples)
}


plt_input2 =  plt_input %>% 
mutate("sample-id" = factor(.$`sample-id`, levels = li_order_samples)) %>%
mutate("Taxon" = factor(.$Taxon, levels = rev(c(taxon_order, "other")))) %>%
mutate("group" = factor(.$group, levels = group_value))
#plt_input2 %>% reshape2::dcast(Taxon ~ `sample-id`, fill.value="value")

n_col = ifelse(length(core_taxon)<45, 1, ceiling(length(core_taxon)/50)+1)

print(n_col)

# samples stacked bar
g1 <- ggplot(plt_input2, aes(x = `sample-id`, y = value, fill = Taxon))+
geom_col(width=1)+
theme_classic()+
theme(text = element_text(family = "serif"), plot.caption = element_text(size=5),
      axis.ticks = element_blank(), panel.spacing.x = unit(1, "mm"), axis.title.x = element_blank(), 
      strip.background.x = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(size=3, angle=90), 
      legend.text = element_text(size=fontsize, lineheight = .5), legend.position = "right", legend.spacing.x = unit(1, "cm"))+
guides(fill= guide_legend(ncol = n_col, title = rank, title.size = 10, size=4,  keywidth = unit(.3, "cm"), keyheight = unit(.3, "cm"), reverse = T, byrow=FALSE))+
labs(y = "Relative abundance [%]", caption  = plt_caption)+
scale_fill_manual(values = taxon_color)+
facet_grid(.~group, scales = "free_x")
#ggsave(paste0(pngfile, ".2.png"), g1, width = plt_width, height = plt_height, dpi = 300)
# group cols
g2<- ggplot(plt_input2)+
geom_bar(mapping = aes(x=`sample-id`, y = 1, fill = group), stat = "identity", width=1)+
theme_void()+
theme(text = element_text(family = "serif"), panel.spacing.x = unit(.1, "mm"))+
guides(fill=guide_legend(ncol = n_col, title = "group", keywidth = unit(.3, "cm"), keyheight = unit(.3, "cm"), reverse = T))+
scale_fill_brewer(palette = 'Set1')+facet_grid(.~group, scales = "free_x")

if(xlabel_out=="FALSE"){
	g1 <- g1+theme(axis.text.x = element_blank())
}
legend <- cowplot::plot_grid(cowplot::get_legend(g1)) #n_col
h1 <- g1+theme(legend.position = "none")
h2 <- g2+theme(legend.position = "none")
plot <- cowplot::plot_grid(h2, h1, align = "v", ncol = 1, axis = "tb", rel_heights = c(.5, 15))
#g<- ggpubr::ggarrange(plot, legend, widths = c(10,8))
#ggsave(pngfile, g, width = plt_width, height = plt_height, dpi = 300, limitsize = FALSE)

plt_legend <- ggplot_gtable(ggplot_build(legend))
plt_legend$widths <- ggplot_gtable(ggplot_build(plot))$widths
plt_legend$heights <- ggplot_gtable(ggplot_build(plot))$heights
ggsave(pngfile, ggpubr::ggarrange(plot, ggplot()+theme_void(), plt_legend, nrow=1, widths = c(1,0.05, 1)), width=plt_width, height=plt_height, dpi=300, limitsize=FALSE)
