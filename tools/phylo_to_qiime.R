setwd("C:/Users/Xiaoping/OneDrive/Desktop/mcm_R/updated_pipeline_May172020/microcosm_analysis/")
rm(list=ls())

library(phyloseq)
library(biomformat)
library(tidyverse)
library(ggplot2)

dir.create("./QIIME2")
outqiime = normalizePath("./QIIME2/")

mydata.16s.pruned = readRDS("./imported_data/pruned_data/mydata_16s_pruned.rds")

filelist.16s = list.files("./Deseq2/16s/")
sigASVs.cultype.16s = data.frame()
for(i in filelist.16s){
        if(grepl("Conventional|Organic|Green", i, ignore.case = T)){
                if(grepl("significantOTU_withfoldChange.csv", i, ignore.case = T)){
                        df = read.csv(file.path("Deseq2/16s/", i))
                        sigASVs.cultype.16s = rbind(sigASVs.cultype.16s, df)
                        
                }
        }
}


sigIDs.16s = as.character(sigASVs.cultype.16s$X)
sig16s = subset_taxa(mydata.16s.pruned, rownames(tax_table(mydata.16s.pruned)) %in% sigIDs.16s)




# otu
otu = as(otu_table(sig16s), "matrix")
otu_biom = make_biom(otu)
write_biom(otu_biom, file.path(outqiime, "sigcultype_otu_16s.biom"))

# tax tabel
tax = as(tax_table(sig16s), "matrix")
tax_cols = colnames(tax)
tax = as.data.frame(tax)
tax[is.na(tax)] = "unassigned"
tax$taxonomy = do.call(paste, c(tax[tax_cols], sep=";"))
tax = tax[, !names(tax) %in% tax_cols, drop = F]
write.table(tax, file.path(outqiime, "tax.tsv"), quote = F, col.names = F, sep = "\t")


# meta data

meta16s = as.data.frame(sample_data(sig16s))

write.table(meta16s, file.path(outqiime, "meta16s_sigcultype.tsv"), sep="\t", col.names = T, row.names = F, quote = F)


########## qiime2 return data

qiim2 = read.table("./QIIME2/qiime2_results_longitudial/important_ASV.tsv", sep = "\t", header = T)
head(qiim2)


# longformat



q2_long = qiim2 %>%
        gather(key = "ASV", value = "important_score", bASV10735:bASV70996)

# append taxonomy
taxdf = as.data.frame(tax_table(sig16s), stringsAsFactors = F)

taxdf[is.na(taxdf)] = "unassigned"

q2_long.mg = merge(q2_long, taxdf, by.x = "ASV", by.y = "row.names")


# top

q2_long.mg %>%
        filter(important_score > 0.01) %>%
        ggplot(aes(x = period, y = important_score, color = cul_type, group = cul_type)) +
        stat_summary(fun = mean, geom = "line") +
        stat_summary(fun = mean, geom = "point") +
        stat_summary(fun.data = mean_se, geom="errorbar", width=0.1)
