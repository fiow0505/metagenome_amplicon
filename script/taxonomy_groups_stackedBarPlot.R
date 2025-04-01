# Author : Kyeongeui Yun (sindy@hunbiome.com)
# 2023. 08. 08
#!/usr/bin/env Rscript
options(warn = -1) # Turn off warning

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
if (TRUE){
  option_list = list(
      make_option(c("-i", "--input"), type="character", 
                  help="taxonomy composition table per group [Bacteria.xlsx]"),
      make_option(c("-s", "--sort_value"), type="character", default="",
                help="sort by group value.  ex.Group1,Group2,Group3"),
      make_option(c("-o", "--output_path"), type="character", default="",
                help="Output path [default %default]"),
      make_option(c("-l", "--level"), type="character", default = "All", help = "taxonomy level: Kingdom(1) ~ Species(7) [default. %default]" ),
      make_option(c("-p", "--prevalence_min"), type="numeric", default = 30, help = "sample found [default. prevalecce > %default]"),
      make_option(c("-r", "--relev_min"), type="numeric", default = 0.001, help = "cutoff average [default. All.mean > %default]"),
      make_option(c("-f", "--fontsize"), type="numeric", default = 8),
      make_option(c("-w", "--width"), type="numeric", default = 6.5),
      make_option(c("-e", "--height"), type="numeric", default=6.5)
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
        
        new_taxon = ifelse(st_num==lev_taxon,li_sp_taxon[lev_taxon], paste0(li_sp_taxon[c(st_num:lev_taxon)], collapse = sep))
        
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


library(tidyverse)
library(ggplot2)
library(readxl)

colFunc <- colorRampPalette(c("lightpink","purple4",  "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))

inputfile = tools::file_path_as_absolute(opts$input)
outputPth = opts$output_path
sort_value = opts$sort_value
width = opts$width
height = opts$height
level = opts$level
fontsize = opts$fontsize
min_prev = opts$prevalence_min
min_relev = opts$relev_min
caption_txt = paste0("*Core microbiome\nPrevalence>", min_prev,"% and Average >", min_relev,"%")
if(sort_value!=""){group_value = unlist(str_split(sort_value, ","))}else{group_value = unique(metadata$group)}
li_rank = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
if(level=="All"){get_rank = li_rank}else{get_rank = li_rank[as.numeric(level)-1]}
print(get_rank)

plt_excel = paste0(outputPth, "/group.core-taxon.xlsx")

li_xlsx = list()
for (rank in get_rank){
    df_relev = read_excel(inputfile, sheet = rank) %>% as.data.frame()
    rownames(df_relev) = df_relev[,1]
    df_relev = df_relev[,-1]

    df_relev_filt = df_relev %>% filter(`prevalence(%)`> min_prev) %>% filter(`All.MEAN`> min_relev)
    get_taxon = ifelse(nrow(df_relev_filt)>30,  df_relev_filt %>% head(30), df_relev_filt) %>% unlist()
    df_relev_filt_get = df_relev_filt %>% filter(Full_lineage %in% get_taxon)
    #df_relev_filt_get$Taxon = rename_taxon(df_relev_filt_get$Full_lineage, ";", TRUE)
	rownames(df_relev_filt_get) = df_relev_filt_get$Taxon

    df_plt = df_relev_filt_get[c("Taxon", "All.MEAN", group_value)]
    rownames(df_plt) = df_plt[,1]
    df_plt = df_plt[,-1]
    
    df_plt1 = df_plt
    df_plt1["Other",] = 100-apply(df_plt1, 2, sum)
    plt1_color_list <- colFunc(nrow(df_plt1)-1)
    plt1_taxon_color <- rev(c("black", plt1_color_list))
    names(plt1_taxon_color) = rownames(df_plt1)
    
    if(apply(df_plt1["Other",], 1, sum)==0){
        df_plt1 = df_plt1[!grepl("Other",rownames(df_plt1)),]
    }
    print(head(df_plt1))
    g1<- df_plt1 %>% rownames_to_column("Taxon") %>% mutate("Taxon" = factor(.$Taxon, levels = rev(.$Taxon))) %>% 
    reshape2::melt(id_vars = "Taxon") %>%
    mutate("variable" = factor(.$variable, levels = c("All.MEAN", group_value))) %>%
    ggplot()+geom_bar(aes(x=`variable`, y = value, fill = Taxon),  stat = "identity")+
    theme_classic()+
    labs(x="", y = "Relative abundance [%]", fill = rank,
         title = "Relative abundance in core microbiome", 
         caption = caption_txt )+
        theme(text = element_text(family = "serif", size = fontsize-2 ), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
              axis.ticks = element_blank(), panel.spacing.x = unit(.1, "mm"), 
              legend.text = element_text(size=fontsize-2, lineheight = .8))+
        scale_fill_manual(values = plt1_taxon_color, guide = guide_legend(reverse = T, ncol=1, keywidth = unit(.3, "cm"), keyheight = unit(.3, "cm")))
    ggsave(paste0(outputPth, "/group.core-taxon-", rank, ".png"),g1, dpi = 300, width = width, height = height)
#	save_data(plt_excel, df_relev_filt_get , rank)
    li_xlsx[[rank]] = df_relev_filt_get

    if(nrow(df_plt1)>11){
        df_plt2 = df_plt %>% head(10)
        df_plt2["Other",] = 100-apply(df_plt2, 2, sum)
        g2<- df_plt2 %>% rownames_to_column("Taxon") %>% mutate("Taxon" = factor(.$Taxon, levels = rev(.$Taxon))) %>% 
        reshape2::melt(id_vars = "Taxon") %>%
        mutate("variable" = factor(.$variable, levels = c("All.MEAN", group_value))) %>%
        ggplot()+geom_bar(aes(x=`variable`, y = value, fill = Taxon),  stat = "identity")+
        theme_classic()+
        labs(x="", y = "Relative abundance [%]", fill = rank,
             title = "Relative abundance : Top 10 & core microbiome", 
             caption = caption_txt)+
            theme(text = element_text(family = "serif", size = fontsize-2),
                  axis.ticks = element_blank(), panel.spacing.x = unit(.1, "mm"), 
                  legend.text = element_text(size=fontsize-2, lineheight = .8))+
            scale_fill_manual(values = plt1_taxon_color, guide = guide_legend(reverse = T, ncol=1, keywidth = unit(.3, "cm"), keyheight = unit(.3, "cm")))
        ggsave(paste0(outputPth, "/group.core-taxon-", rank,".top10.png"),g2, dpi = 300, width = width, height = height)

        
    }
    
}
library(openxlsx)
write.xlsx(li_xlsx, file = plt_excel) 
