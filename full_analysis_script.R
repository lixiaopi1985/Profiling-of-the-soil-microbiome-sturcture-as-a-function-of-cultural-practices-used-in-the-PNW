setwd("your current working path")
rm(list = ls())


library(phyloseq)
library(tidyverse)
library(tools)
library(vegan)
library(ampvis2)
library(microbiome)
library(microbiomeutilities)
library(ggrepel)
library(vegan)
library(ggplot2)
library(tidyverse)
library(ggvegan)
library(ggrepel)
library(ggthemes)
library(RColorBrewer)
library(stringr)
library(latex2exp)
library(pals)
library(lmerTest)
library(car)
library(indicspecies)
library(SpiecEasi)
library(igraph)
library(DESeq2)
library(emmeans)
library(pheatmap)


sessionInfo()

devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

dir.create("./field_soil_analysis_output")


# source data folder contains exported biom files from QIIME2


#######################################################################################################################
# The first section runs the analysis for field soil samples
#
#######################################################################################################################


# import data from QIIME2 output

#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# 16S
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

fdata_16s = import_biom("./source_data/onsite/16s/feature-table-tax.biom")
tree.16s = read_tree("./source_data/onsite/16s/tree.nwk")
phy_tree(fdata_16s) = tree.16s

meta_16s = data.frame(sample_data(fdata_16s))

# change the taxonomy rank
tax_table(fdata_16s) = gsub("D_[0-9]__", "", tax_table(fdata_16s)) # remove prefix to make it look better
colnames(tax_table(fdata_16s)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # fixed the rank names

numeric_cols = c('Meloidogyne_hapla',	'Meloidogyne_chitwoodi',	'Pratylenchus_penetrans',	'Pratylenchus_neglectus',	'Pratylenchus_thornei',	'Paratrichodorus_sp',	'Tylenchorhynchus_sp','Paratylenchus_sp',	'Helicotylenchus_sp','Mesocriconema_sp',	'Ditylenchus_sp',	'Xiphinema_sp',	'Hemicycliophora_sp',	'vert_count',	'bd_count', 'pH')


colnames(meta_16s)[which(colnames(meta_16s) == "MesocriconemaÂ.sp")]="Mesocriconema_sp"

# change pathogen count to numeric
for(i in numeric_cols){
        meta_16s[i] = lapply(meta_16s[i], as.numeric)
}

factor_cols = c("cul_type", "soil_compname","soil_taxclass","soil_taxgrtgroup","soil_taxorder","soil_taxsubgrp","soil_taxsuborder", paste0("y", seq(1998, 2018)))

for(i in factor_cols){
        meta_16s[i] = lapply(meta_16s[i], as.factor)
}

# colum names to lower case
colnames(meta_16s) = tolower(colnames(meta_16s))

cols_16s = colnames(meta_16s)

meta_16s$SampleID = row.names(meta_16s)

meta_16s = meta_16s[c("SampleID", cols_16s)]

sample_data(fdata_16s) = meta_16s

# add bASV to indicate bacterial community
taxa_names(fdata_16s) = paste0("bASV", seq_along(taxa_names(fdata_16s)))



#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# Fungi
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

fdata_its = import_biom("./source_data/onsite/its/feature-table-tax.biom")
tree.its = read_tree("./source_data/onsite/its/tree.nwk")
phy_tree(fdata_its) = tree.its

meta_its = data.frame(sample_data(fdata_its))

tax_table(fdata_its) = gsub("k__|p__|c__|o__|f__|g__|s__", "", tax_table(fdata_its))
colnames(tax_table(fdata_its)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

numeric_cols = c('Meloidogyne_hapla',	'Meloidogyne_chitwoodi',	'Pratylenchus_penetrans',	'Pratylenchus_neglectus',	'Pratylenchus_thornei',	'Paratrichodorus_sp',	'Tylenchorhynchus_sp','Paratylenchus_sp',	'Helicotylenchus_sp','Mesocriconema_sp',	'Ditylenchus_sp',	'Xiphinema_sp',	'Hemicycliophora_sp',	'vert_count',	'bd_count', 'pH')


colnames(meta_its)[which(colnames(meta_its) == "MesocriconemaÂ.sp")]="Mesocriconema_sp"
 
# change pathogen count to numeric
for(i in numeric_cols){
        meta_its[i] = lapply(meta_its[i], as.numeric)
}

factor_cols = c("cul_type", "soil_compname","soil_taxclass","soil_taxgrtgroup","soil_taxorder","soil_taxsubgrp","soil_taxsuborder", paste0("y", seq(1998, 2018)))

for(i in factor_cols){
        meta_its[i] = lapply(meta_its[i], as.factor)
}


colnames(meta_its) = tolower(colnames(meta_its))

cols_its = colnames(meta_its)

meta_its$SampleID = row.names(meta_its)

meta_its = meta_its[c("SampleID", cols_its)]

sample_data(fdata_its) = meta_its

taxa_names(fdata_its) = paste0("fASV", seq_along(taxa_names(fdata_its)))

#((((((((((((((((((((((((()))))))))))))))))))))))))
# save
#((((((((((((((((((((((((()))))))))))))))))))))))))

# dir.create("./field_soil_analysis_output/data/", recursive = T)
# saveOrig.16s = file.path("./field_soil_analysis_output/data/", "fdata_16s.rds")
# saveRDS(fdata_16s, saveOrig.16s)
# 
# 
# saveOrig.its = file.path("./field_soil_analysis_output/data/", "fdata_its.rds")
# saveRDS(fdata_its, saveOrig.its)


##################################################
# Get summary of the reads and reads distribution
##################################################



# total counts

sum(sample_sums(fdata_16s)) # 3724679
sum(sample_sums(fdata_its)) # 7043461

# min depth
min(sample_sums(fdata_16s)) # min 40026
min(sample_sums(fdata_its)) # min 31323


# rare taxa
min(taxa_sums(fdata_16s)) # min 1
min(taxa_sums(fdata_its)) # min 1


# prune data

# remove NA phylum

fdata_16s.prune = subset_taxa(fdata_16s, !is.na(Phylum))
fdata_its.prune = subset_taxa(fdata_its, !is.na(Phylum))

# remove singleton and doubleton and OTU abundance less than 10 across all samples

fdata_16s.prune = prune_taxa(taxa_sums(fdata_16s.prune) > 10, fdata_16s.prune)
fdata_its.prune = prune_taxa(taxa_sums(fdata_its.prune) > 10, fdata_its.prune)


# add crop rotation data


library(tidyverse)


meta_16s.prune = as(sample_data(fdata_16s.prune), "data.frame")
meta_its.prune = as(sample_data(fdata_its.prune), "data.frame")

crops.16s = meta_16s.prune %>%
        select(SampleID, field, "y1998":"y2018") %>%
        mutate_all(list(~na_if(., ""))) %>%
        gather(key="rotYear", value="crops", "y1998":"y2018") %>%
        group_by(field) %>%
        mutate(totalCrops = n_distinct(crops, na.rm = T)) %>%
        mutate(totalRotYrs = n_distinct(rotYear[!is.na(crops)]), totalYrs_Pota = sum(crops=="POTATO", na.rm = T)/n_distinct(SampleID)) %>%
        mutate(percPota = totalYrs_Pota / totalRotYrs ) %>%
        spread(rotYear, crops)


new_meta.16s.prune = merge(meta_16s.prune, crops.16s[, c("SampleID", "totalCrops", "totalRotYrs", "totalYrs_Pota", "percPota")], by = "SampleID")

new_meta.its.prune = merge(meta_its.prune, crops.16s[, c("SampleID", "totalCrops", "totalRotYrs", "totalYrs_Pota", "percPota")], by = "SampleID")

rownames(new_meta.16s.prune) = new_meta.16s.prune$SampleID
sample_data(fdata_16s.prune) = new_meta.16s.prune

rownames(new_meta.its.prune) = new_meta.its.prune$SampleID
sample_data(fdata_its.prune) = new_meta.its.prune

# generate pathogen data

path_meta = new_meta.16s.prune %>%
        select(-c("field_coord", "coll_coord", "cords_note", "quant_reading", "sample_or_control", "soil_collected_date", "mean_day_air_temp"))

pathop = file.path("./field_soil_analysis_output/pathogens/", "meta_pathogens.tsv")
# write.table(path_meta, file=pathop, sep = "\t", row.names = F)


# use microbiomeutility to find the best taxonomy rank for NA
fbst16s = format_to_besthit(fdata_16s.prune)
fbstITS = format_to_besthit(fdata_its.prune)

fbst16s_taxtb = tax_table(fbst16s)[, 1:7]
rownames(fbst16s_taxtb) = gsub("OTU-|:.*$", "", rownames(fbst16s_taxtb))

fbstits_taxtb = tax_table(fbstITS)[, 1:7]
rownames(fbstits_taxtb) = gsub("OTU-|:.*$", "", rownames(fbstits_taxtb))


tax_table(fdata_16s.prune) = fbst16s_taxtb
tax_table(fdata_its.prune) = fbstits_taxtb


# saveRDS(fdata_16s.prune, file.path("./field_soil_analysis_output/data/", "fdata_16s_pruned.rds"))
# saveRDS(fdata_its.prune, file.path("./field_soil_analysis_output/data/", "fdata_its_pruned.rds"))

########################################################################

# Analysis pathogen association

########################################################################

rm(list=ls())

pathogen_data = read.table("./field_soil_analysis_output/pathogens/meta_pathogens.tsv", header = T)


pathogen_cols = tolower(c('Meloidogyne_hapla',	'Meloidogyne_chitwoodi',	'Pratylenchus_penetrans',	'Pratylenchus_neglectus',	'Pratylenchus_thornei',	'Paratrichodorus_sp',	'Tylenchorhynchus_sp','Paratylenchus_sp',	'Helicotylenchus_sp','Mesocriconema_sp',	'Ditylenchus_sp',	'Xiphinema_sp',	'Hemicycliophora_sp',	'vert_count',	'bd_count'))

patho_labels = gsub("_", " ", pathogen_cols)
patho_labels[(length(patho_labels)-1):length(patho_labels)] = c("Verticillium dahliae", "Colletorichum coccodes")
patho_labels = str_to_sentence(patho_labels)

nemas = tolower(c('Meloidogyne_hapla',	'Meloidogyne_chitwoodi',	'Pratylenchus_penetrans',	'Pratylenchus_neglectus',	'Pratylenchus_thornei',	'Paratrichodorus_sp',	'Tylenchorhynchus_sp','Paratylenchus_sp',	'Helicotylenchus_sp','Mesocriconema_sp',	'Ditylenchus_sp',	'Xiphinema_sp',	'Hemicycliophora_sp'))


########## counting pathogen in each field

colPalette = colorRampPalette(pals::alphabet())
index2 = which(!colnames(pathogen_data) %in% pathogen_cols)

data_copy = pathogen_data[c(pathogen_cols, colnames(pathogen_data)[index2])]
colnames(data_copy)[colnames(data_copy) %in% pathogen_cols] = patho_labels

(path_bcul = data_copy %>%
                gather(key="Pathg", value="counts", all_of(patho_labels)) %>%
                group_by(cul_type, Pathg) %>%
                summarise(total_counts = sum(counts)) %>% filter(total_counts>0))

# log 2 transform
(path_bcul2 = path_bcul %>%
                mutate(se = sd(log2(total_counts))/sqrt(n()), log2T = log2(total_counts)) %>%
                mutate(position=rank(-log2T)))


# use total count for original data
(path_m = path_bcul2 %>%
                select(cul_type, Pathg, total_counts) %>%
                spread(Pathg, total_counts))


# make matrix
rowNames = path_m$cul_type
path_mtrx = as.matrix(path_m[, 2:ncol(path_m)])
rownames(path_mtrx) = rowNames
# fill NA with 0
path_mtrx[is.na(path_mtrx)] = 0
path_mtrx[1:3, 1:7]
(chitest = chisq.test(path_mtrx))
chitestp = chitest$p.value

### plot counts

(
        pathogen_plot = ggplot(path_bcul2, aes(x=cul_type, y=log2T, fill=reorder(Pathg, -log2T), group = position)) +
                geom_bar(stat="identity", position=position_dodge(.9), width=.9, alpha = 0.8) +
                geom_errorbar(aes(ymin=log2T - se, ymax = log2T + se), position = position_dodge(.9), size=0.5, width=0.2) +
                theme_bw() +
                annotate('text', x= 3, y = 14, label = TeX('$\\chi^2$: df = 14,$\\textit{p}$ < 0.0001')) +
                scale_fill_manual(name = "Pathogens", values = colPalette(length(patho_labels))) +
                theme(legend.text = element_text(face="italic"), legend.position = "bottom") +
                labs(x="Cropping systems", y = "Pathogen counts (log2 tranformed)") +
                scale_x_discrete(labels = c("Conventional", "Mixed", "Organic"))
)



# CCA


patho_abund = pathogen_data[, pathogen_cols]
rownames(patho_abund) = pathogen_data$SampleID
patho_abund

pathogen_cols
colnames(pathogen_data)

env = pathogen_data %>%
        select(-c("SampleID", "field", "y1998":"y2018", all_of(pathogen_cols))) %>%
        replace_na(replace = list(percPota=0))

rownames(env) = pathogen_data$SampleID

# select variable we are interested

env2 = env %>%
        select(cul_type, ph, soil_compname, totalCrops, totalRotYrs, percPota)


env2
str(patho_abund)
sum(!complete.cases(patho_abund))

str(env2)
sum(!complete.cases(env2))
tail(env2)

#log transform
otu_log = decostand(patho_abund[rowSums(patho_abund) >0, ], method = "log")



env.2 = env2[rownames(otu_log), ]

decorana(otu_log)

range(env.2[env.2$cul_type == "Green", ]$percPota)
range(env.2[env.2$cul_type == "Organic", ]$percPota)
range(env.2[env.2$cul_type == "Conventional", ]$percPota)


patho.cca = vegan::cca(otu_log ~ ., data = env.2)

patho.cca



screeplot(patho.cca)
set.seed(123)
anova(patho.cca) # significant
anova(patho.cca, by="terms") # ph, soil type, crop diversity, not significant
anova(patho.cca, by="margin") #ph, soil type not significant
anova(patho.cca, by="axis") # CCA1 and CCA2

# get better model

final.patho.cca = ordistep(cca(otu_log ~1, data=env.2), scope = formula(patho.cca), direction = "forward", pstep=1000)


vif.cca(final.patho.cca)

# test significance
set.seed(123)
anova.cca(final.patho.cca, step=1000)
anova.cca(final.patho.cca, by="terms", step=1000)
anova.cca(final.patho.cca, by="axis", step = 1000) # now CCA1, CCA2 and CCA3


final.patho.cca # constrained total 0.4624


RsquareAdj(final.patho.cca) # 0.46

#CCA1 eigen value
0.6769/2.3779
#CCA2 eigen value
0.2953/2.3779
#CCA3 eigen value
0.1273/2.3779

# vegan plot

# plot(final.patho.cca, display = "sites",type = "n", scaling = 3, xlim = c(-5, 5), ylim=c(-3, 8))
# points(final.patho.cca, display = "sites", scaling = 3, col="dodgerblue", pch=19)
# points(final.patho.cca, display = "sp", scaling = 3, col=538, pch=8, cex=1.5)
# text(final.patho.cca, display = "bp", scaling = 3, cex=1)
# text(final.patho.cca, display = "sp",scaling = 3, col="maroon", cex=0.8)


cca.df.sp = fortify(final.patho.cca, display = "sp") # get species infor
cca.df.sa = fortify(final.patho.cca, display = "sites") %>% bind_cols(env.2) # get sample info
cca.df.bp = fortify(final.patho.cca, display = "bp") # get biplot info
cca.df.cn = fortify(final.patho.cca, display = "cn") # get centroids info

# autoplot(final.patho.cca, layers = c("species", "sites", "biplot"))


axes = c("CCA1", "CCA2")


## for the arrow, if specifiy bp and cn, then have to consider the plot size scale in order to match plot.cca(display=c("bp", "cn))
# mul = ggvegan:::arrowMul(cca.df.bp[, axes], rbind(cca.df.sa[, axes], cca.df.sp[, axes]))



colPalette = colorRampPalette(brewer.pal(11, "Paired"))


summ = summary(final.patho.cca)


xlab.cca = paste0("CCA1 ( ", round(summ$cont$importance[2,1], 4)*100, "% of total /", round(summ$concont$importance[2,1], 4)*100, "% of constrained )")
ylab.cca = paste0("CCA2 ( ", round(summ$cont$importance[2,2], 4)*100, "% of total /", round(summ$concont$importance[2,2], 4)*100, "% of constrained )")



# biplot

## arrows for biplot
# mul = ggvegan:::arrowMul(cca.df.bp[, axes], rbind(cca.df.sa[, axes]))
# arrows = cca.df.bp[, axes]
# arrows = arrows*mul
# 
# g = ggplot()
# g +
#         xlim(c(-2,3))+
#         ylim(c(-3, 5)) +
#         geom_hline(aes(yintercept=0), linetype="dashed", size=0.8, color="grey") +
#         geom_vline(aes(xintercept=0), linetype="dashed", size=0.8, color="grey") +
#         geom_point(data=cca.df.sa, aes(x=CCA1, y=CCA2, color=Label), size=3, alpha=0.6) +
#         geom_point(data=cca.df.sp, aes(x=CCA1, y=CCA2), shape="triangle", color="darkgreen", size=4, alpha=0.6) +
#         geom_text_repel(data=cca.df.sp, aes(x=CCA1, y=CCA2, label=Label), color="darkgreen", nudge_x = 0.5) +
#         geom_segment(data=arrows, aes(x=0, y=0, xend=CCA1, yend=CCA2), size=0.6, color="maroon", arrow = arrow(length = unit(0.01, "npc"))) +
#         geom_text(data=arrows, aes(x=CCA1-0.1, y=CCA2, label=cca.df.bp$Label), color="maroon") +
#         coord_fixed() +
#         theme_few() +
#         labs(color = "Samples", x=xlab.cca, y=ylab.cca) +
#         geom_label(aes(x=2, y=4, label="Scaling=1")) +
#         scale_color_manual(values = colPalette(nrow(env.2)))


# with centroids

mul = ggvegan:::arrowMul(cca.df.bp[, axes], rbind(cca.df.sa[, axes], cca.df.sp[, axes]))
arrows = cca.df.bp[, axes]
arrows = arrows*mul

levels(cca.df.cn$Label) = c("C", "M", "O")

colMap = c("#fc032c", "#030ffc", "#32a852")




g = ggplot()
g +
        xlim(c(-2,3))+
        ylim(c(-3, 5)) +
        geom_hline(aes(yintercept=0), linetype="dashed", size=0.8, color="grey") +
        geom_vline(aes(xintercept=0), linetype="dashed", size=0.8, color="grey") +
        geom_point(data=cca.df.sa, aes(x=CCA1, y=CCA2, color = cul_type), size=4, alpha=0.6) +
        geom_point(data=cca.df.sp, aes(x=CCA1, y=CCA2, shape=as.factor(Label)), color="black", size=4, alpha=0.8) +
        geom_text_repel(data = cca.df.cn, aes(x=CCA1, y=CCA2, label=Label), color="blue") +
        geom_segment(data=arrows[3, ], aes(x=0, y=0, xend=CCA1, yend=CCA2), arrow = arrow(length = unit(0.01, "npc")), color="maroon") +
        geom_text(data=arrows[3, ], aes(x=CCA1-0.1, y=CCA2, label=cca.df.bp[3,]$Label), color="maroon") +
        coord_fixed() +
        theme_few() +
        labs(color = "Cropping systems", x=xlab.cca, y=ylab.cca) +
        geom_text(aes(x=2, y=4.5, label = "Total inertia = 2.3779\nConstrained inertia = 1.0996")) +
        scale_shape_manual(name = "Pathogens", values = c(19, 17, 15, 9, 8, 10, 12, 7))

# dir.create("./field_soil_analysis_output/pathogens")
ggsave(file.path("./field_soil_analysis_output/pathogens/", "cca_plot.svg"), width = 12, height = 10)
# final.patho.cca



crops = env16s %>%
        select(field, y1998:y2018) %>%
        gather("years", "crops", y1998:y2018) %>%
        distinct()

crops %>%
        group_by(field, crops) %>%
        tally() %>%
        top_n(3) %>%
        arrange(-n, .by_group = T) %>%
        write.csv("./field_soil_analysis_output/pathogens/topcrops.csv", quote = F, row.names = F)

####################################################

# analyzing alpha diversity

# if saved

rm(list = ls())

fdata_16s.prune = readRDS("./field_soil_analysis_output/data/fdata_16s_pruned.rds")
fdata_its.prune = readRDS("./field_soil_analysis_output/data/fdata_its_pruned.rds")


trial = 100
####################################################

#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# Bacteria alpha diversity
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

# rarefy



observed_mat_16s = matrix(nrow = nrow(sample_data(fdata_16s.prune)), ncol = trial)
row.names(observed_mat_16s) = sample_names(fdata_16s.prune)
shannon_mat_16s = matrix(nrow = nrow(sample_data(fdata_16s.prune)), ncol = trial)
row.names(shannon_mat_16s) = sample_names(fdata_16s.prune)
invSimpson_mat_16s = matrix(nrow = nrow(sample_data(fdata_16s.prune)), ncol = trial)
row.names(invSimpson_mat_16s) = sample_names(fdata_16s.prune)

min_lib_16s = min(sample_sums(fdata_16s.prune))

set.seed(100)
for(i in 1:trial){

        cat("------------------------------------------------------------------------------------")
        print(i)
        ra_16s = rarefy_even_depth(fdata_16s.prune, sample.size = min_lib_16s, replace=F)

        est_16s = as.matrix(estimate_richness(ra_16s, measures = c("Observed", "Shannon", "InvSimpson")))
        observed_mat_16s[, i] = as.numeric(est_16s[,1])
        shannon_mat_16s[,i] = as.numeric(est_16s[,2])
        invSimpson_mat_16s[,i] = as.numeric(est_16s[,3])

}

SampleID = row.names(observed_mat_16s)
mean = apply(observed_mat_16s, 1, mean)
sd = apply(observed_mat_16s, 1, sd)
measure = rep("Observed", nsamples(fdata_16s.prune))
observed_stats_16s = data.frame(SampleID, mean, sd, measure)

SampleID = row.names(shannon_mat_16s)
mean = apply(shannon_mat_16s, 1, mean)
sd = apply(shannon_mat_16s, 1, sd)
measure = rep("Shannon", nsamples(fdata_16s.prune))
shannon_stats_16s = data.frame(SampleID, mean, sd, measure)


SampleID = row.names(invSimpson_mat_16s)
mean = apply(invSimpson_mat_16s, 1, mean)
sd = apply(invSimpson_mat_16s, 1, sd)
measure = rep("invSimpson", nsamples(fdata_16s.prune))
invSimpson_stats_16s = data.frame(SampleID, mean, sd, measure)


aldiv_est_16s = rbind(observed_stats_16s, shannon_stats_16s, invSimpson_stats_16s)

a_diversity_16s = merge(data.frame(sample_data(fdata_16s.prune)), aldiv_est_16s, by="SampleID")

a_diversity_16s$percPota[is.na(a_diversity_16s$percPota)] = 0

# dir.create("./field_soil_analysis_output/alpha_diversity")
saveRDS(a_diversity_16s, file.path("./field_soil_analysis_output/alpha_diversity/fdata_16s_alpha.rds"))


#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# fungai alpha diversity
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

observed_mat_its = matrix(nrow = nrow(sample_data(fdata_its.prune)), ncol = trial)
row.names(observed_mat_its) = sample_names(fdata_its.prune)
shannon_mat_its = matrix(nrow = nrow(sample_data(fdata_its.prune)), ncol = trial)
row.names(shannon_mat_its) = sample_names(fdata_its.prune)
invSimpson_mat_its = matrix(nrow = nrow(sample_data(fdata_its.prune)), ncol = trial)
row.names(invSimpson_mat_its) = sample_names(fdata_its.prune)

min_lib_its = min(sample_sums(fdata_its.prune))

set.seed(101)

for(i in 1:trial){

        print(i)
        ra_its = rarefy_even_depth(fdata_its.prune, sample.size = min_lib_its, replace=F)

        est_its = as.matrix(estimate_richness(ra_its, measures = c("Observed", "Shannon", "InvSimpson")))
        observed_mat_its[, i] = as.numeric(est_its[,1])
        shannon_mat_its[,i] = as.numeric(est_its[,2])
        invSimpson_mat_its[,i] = as.numeric(est_its[,3])

}

SampleID = row.names(observed_mat_its)
mean = apply(observed_mat_its, 1, mean)
sd = apply(observed_mat_its, 1, sd)
measure = rep("Observed", nsamples(fdata_its.prune))
observed_stats_its = data.frame(SampleID, mean, sd, measure)

SampleID = row.names(shannon_mat_its)
mean = apply(shannon_mat_its, 1, mean)
sd = apply(shannon_mat_its, 1, sd)
measure = rep("Shannon", nsamples(fdata_its.prune))
shannon_stats_its = data.frame(SampleID, mean, sd, measure)


SampleID = row.names(invSimpson_mat_its)
mean = apply(invSimpson_mat_its, 1, mean)
sd = apply(invSimpson_mat_its, 1, sd)
measure = rep("invSimpson", nsamples(fdata_its.prune))
invSimpson_stats_its = data.frame(SampleID, mean, sd, measure)

alpha_est_its = rbind(observed_stats_its, shannon_stats_its, invSimpson_stats_its)

a_diversity_its = merge(data.frame(sample_data(fdata_its.prune)), alpha_est_its, by="SampleID")
a_diversity_its$percPota[is.na(a_diversity_its$percPota)] = 0

saveRDS(a_diversity_its, file.path("./field_soil_analysis_output/alpha_diversity/fdata_its_alpha.rds"))

# read from saved data


# combined diversity from bacteria and fungi



a_diversity_16s$community = "Bacteria"
a_diversity_its$community = "Fungi"

col_to_keep = c("SampleID", "field", "cul_type", "soil_compname", "soil_taxorder", "soil_taxsubgrp", "soil_taxclass", "soil_taxsuborder", "totalCrops", "totalRotYrs", "totalYrs_Pota", "percPota", "mean", "sd", "measure", "ph", "community")

diversity.all = rbind(a_diversity_16s[, col_to_keep], a_diversity_its[, col_to_keep])

head(diversity.all)

saveRDS(diversity.all, file.path("./field_soil_analysis_output/alpha_diversity/diversity_all.rds"))



# changed to exponential Shannon

diversity.all$mean[diversity.all$measure == "Shannon"] = exp(diversity.all$mean[diversity.all$measure == "Shannon"])

diversity.all$measure = as.character(diversity.all$measure)

diversity.all$measure[diversity.all$measure == "Shannon"] = "e^Shannon"
diversity.all$measure = factor(diversity.all$measure, levels = c("Observed", "e^Shannon", "invSimpson"))

saveRDS(diversity.all, file.path("./field_soil_analysis_output/alpha_diversity/diversity_all_modified.rds"))

################################################
# alpha diversity: factors that may affect it
################################################
rm(list=ls())

diversity.all = readRDS("./field_soil_analysis_output/alpha_diversity/diversity_all_modified.rds")
diversity.all$measure

# look at cultural types
diversity.all %>%
        ggplot(aes(x=cul_type, y = mean, fill = cul_type)) +
        geom_boxplot() +
        geom_jitter(width = 0.2) +
        facet_wrap(community ~ measure,scales = "free_y", nrow = 2) +
        theme_bw() +
        labs(fill = "Cropping systems", y = "Alpha diversity indices", x = "Cropping systems") +
        scale_x_discrete(labels = c("C", "M", "O")) +
        scale_fill_discrete(labels = c("Conventional", "Mixed", "Organic"))

ggsave("./field_soil_analysis_output/alpha_diversity/alpha_cultypes.svg", device = "svg", width = 10, height = 8)
# look at soil types

diversity.all %>%
        ggplot(aes(x= soil_compname, y = mean, fill = soil_compname)) +
        geom_boxplot() +
        geom_jitter() +
        facet_wrap(measure ~ community,scales = "free_y", ncol = 2) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        labs(fill = "Soil types", y = "Alpha diversity indices", x = "Soil types")

ggsave("./field_soil_analysis_output/alpha_diversity/alpha_soiltypes.svg", device = "svg", width = 8, height = 10)
# look at total crops
diversity.all %>%
        ggplot(aes(x = totalCrops, y = mean, color = cul_type)) +
        geom_point() +
        geom_smooth(method = "lm", color = "black", fill = "grey90") +
        facet_wrap(measure ~ community,scales = "free_y", ncol = 2) +
        theme_bw()+
        labs(color = "Cultural practices", x = "Total numbers of crops rotated", y = "Alpha diversity indices")

ggsave("./field_soil_analysis_output/alpha_diversity/alpha_totalCrops.svg", device = "svg", width = 8, height = 10)

# look at potato perc
diversity.all %>%
        ggplot(aes(x = percPota, y = mean, color = cul_type)) +
        geom_point() +
        geom_smooth(method = "lm", color = "black", fill = "grey90") +
        facet_wrap(measure ~ community,scales = "free_y", ncol = 2) +
        theme_bw() +
        labs(color = "Cultural practices", x = "Potato rotation proportion", y = "Alpha diversity indices")

ggsave("./field_soil_analysis_output/alpha_diversity/alpha_potatoPerc.svg", device = "svg", width = 8, height = 10)

# look at years of rotation
diversity.all %>%
        ggplot(aes(x = totalRotYrs, y = mean, color = cul_type)) +
        geom_point() +
        geom_smooth(method = "lm", color = "black", fill = "grey90") +
        facet_wrap(measure ~ community,scales = "free_y", ncol = 2) +
        theme_bw() +
        labs(color = "Cultural practices", x = "Total years of rotations", y = "Alpha diversity indices")

ggsave("./field_soil_analysis_output/alpha_diversity/alpha_yearsofrotation.svg", device = "svg", width = 8, height = 10)
# look at pH
diversity.all %>%
        ggplot(aes(x = ph, y = mean, color = cul_type)) +
        geom_point() +
        geom_smooth(method = "lm", color = "black", fill = "grey90") +
        facet_wrap(measure ~ community,scales = "free_y", ncol = 2) +
        theme_bw() +
        labs(color = "Cultural practices", x = "pH values", y = "Alpha diversity indices")

ggsave("./field_soil_analysis_output/alpha_diversity/alpha_pH.svg", device = "svg", width = 8, height = 10)

# look at each field

# diversity.all %>%
#         ggplot(aes(x = field, y = mean, color = cul_type)) +
#         geom_boxplot() +
#         facet_wrap(measure ~ community,scales = "free_y", ncol = 2) +
#         theme_bw() +
#         theme(axis.text.x = element_text(angle=90, hjust = 1))


############################################
# statistical analyzing alpha diversity
############################################
rm(list = ls())

library(nlme)
diversity.all = readRDS("./field_soil_analysis_output/alpha_diversity/diversity_all_modified.rds")
# 16S
head(diversity.all)

diversity.all$rep = gsub(".*_", "", diversity.all$SampleID)
diversity.all$field = as.factor(diversity.all$field)


head(diversity.all)

contrasts(diversity.all$cul_type) = "contr.sum"
contrasts(diversity.all$soil_compname) = "contr.sum"
contrasts(diversity.all$field) = "contr.sum"

obs.16s = subset(diversity.all, measure == "Observed" & community == "Bacteria")
shan.16s = subset(diversity.all, measure == "e^Shannon" & community == "Bacteria")
invS.16s = subset(diversity.all, measure == "invSimpson" & community == "Bacteria")


obs.its = subset(diversity.all, measure == "Observed" & community == "Fungi")
shan.its = subset(diversity.all, measure == "e^Shannon" & community == "Fungi")
invS.its = subset(diversity.all, measure == "invSimpson" & community == "Fungi")

colnames(obs.16s)

#------------------------------------------
boxplot(mean ~ field, data = obs.16s)
summary(aov(mean ~ rep + Error(field/rep), data = obs.16s))

summary(lmer(mean ~ 1 + (1|field), data = obs.16s))

# update remove totalRotYrs in the mixed model

lme.obs.16s = lmer(mean ~ cul_type + soil_compname  + totalCrops + ph + percPota + (1|field), data = obs.16s)

lmer(mean ~ 1 + (1|field/rep), data=obs.16s)
lmer(mean ~ 1 + (1|field) + (1|field:rep), data=obs.16s)

summary(lme.obs.16s)

plot(lme.obs.16s)
qqnorm(resid(lme.obs.16s))
qqline(resid(lme.obs.16s))
hist(resid(lme.obs.16s))

Anova(lme.obs.16s, type=3, test.statistic = "F")

#-----------------
boxplot(mean ~ field, data = shan.16s)
summary(aov(mean ~ rep + Error(field/rep), data = shan.16s))

summary(lmer(mean ~ 1 + (1|field), data = shan.16s))
55159/(55159 + 29540)

lme.shan.16s = lmer(mean ~ cul_type + soil_compname  + totalCrops + ph + percPota + (1|field), data = shan.16s)

summary(lme.shan.16s)
plot(lme.shan.16s)
qqnorm(resid(lme.shan.16s))
qqline(resid(lme.shan.16s))
hist(resid(lme.shan.16s))

Anova(lme.shan.16s, type=3, test.statistic = "F")


#-----------------------
boxplot(mean ~ field, data = invS.16s)
summary(aov(mean ~ rep + Error(field/rep), data = invS.16s))

summary(lmer(mean ~ 1 + (1|field), data = invS.16s))

36819/(36819+6492)


lme.invS.16s = lmer(mean ~ cul_type + soil_compname  + totalCrops + ph + percPota + (1|field), data = invS.16s)



summary(lme.invS.16s)
plot(lme.invS.16s)
qqnorm(resid(lme.invS.16s))
qqline(resid(lme.invS.16s))
hist(resid(lme.invS.16s))

Anova(lme.invS.16s, type=3, test.statistic = "F")



## its

boxplot(mean ~ field, data = obs.its)
summary(aov(mean ~ rep + Error(field/rep), data = obs.its))

summary(lmer(mean ~ 1 + (1|field), data = obs.its))

1983/(1983+5178)


lme.obs.its = lmer(mean ~ cul_type + soil_compname  + totalCrops + ph + percPota + (1|field), data = obs.its)

summary(lme.obs.its)
plot(lme.obs.its)
qqnorm(resid(lme.obs.its))
qqline(resid(lme.obs.its))
hist(resid(lme.obs.its))

Anova(lme.obs.its, type=3, test.statistic = "F")

#----------

boxplot(mean ~ field, data = shan.its)
summary(aov(mean ~ rep + Error(field/rep), data = shan.its))

summary(lmer(mean ~ 1 + (1|field), data = shan.its))

51.41/(51.41+561.88)


lme.shan.its = lmer(mean ~ cul_type + soil_compname  + totalCrops + ph + percPota + (1|field), data = shan.its)

summary(lme.shan.its)
plot(lme.shan.its)
qqnorm(resid(lme.shan.its))
qqline(resid(lme.shan.its))
hist(resid(lme.shan.its))

Anova(lme.shan.its, type=3, test.statistic = "F")

#---------

boxplot(mean ~ field, data = invS.its)
summary(aov(mean ~ rep + Error(field/rep), data = invS.its))

summary(lmer(mean ~ 1 + (1|field), data = invS.its))

9.307/(9.307+162.435)

lme.invS.its = lmer(mean ~ cul_type + soil_compname  + totalCrops + ph + percPota + (1|field), data = invS.its)

summary(lme.invS.its)
plot(lme.invS.its)
qqnorm(resid(lme.invS.its))
qqline(resid(lme.invS.its))
hist(resid(lme.invS.its))

Anova(lme.invS.its, type=3, test.statistic = "F")


######################################################################
# beta diversity

rm(list = ls())

fdata_16s.prune = readRDS("./field_soil_analysis_output/data/fdata_16s_pruned.rds")
fdata_its.prune = readRDS("./field_soil_analysis_output/data/fdata_its_pruned.rds")
source("./tools/plotCAP2.R")

dir.create("./field_soil_analysis_output/beta_diversity")
#####################################################################



# transform to relative abundance

fdata16s.norm = fdata_16s.prune %>%
        transform_sample_counts(function(x)x/sum(x))
fdataITS.norm = fdata_its.prune %>%
        transform_sample_counts(function(x)x/sum(x))

otu.16s = data.frame(t(otu_table(fdata16s.norm))) # relative abundance
env.16s = data.frame(sample_data(fdata16s.norm))


otu.its = data.frame(t(otu_table(fdataITS.norm))) # relative abundance
env.its = data.frame(sample_data(fdataITS.norm))



######################################
# CAP
###################################

#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# Bacteria
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

cap.16s = capscale(otu.16s ~cul_type+soil_compname+totalCrops+totalRotYrs+percPota+ph, data = env.16s, distance = "bray", add = T)

vif.cca(cap.16s) # totalRotYrs and cultype > 10, so remove totalRotYrs

cap.16s = capscale(otu.16s ~cul_type+soil_compname+totalCrops+percPota+ph, data = env.16s, distance = "bray", add=T)
vif.cca(cap.16s)

cap.16s

# check this model
set.seed(123)
anova(cap.16s, permutations = 1000) # 0.001
anova(cap.16s, by="terms", permutations = 1000) # ph is not significant
anova(cap.16s, by="axis", permutations = 1000) # first 4

# get parsimony model
# by AIC
cap.16s.select = ordistep(capscale(otu.16s ~ 1, data = env.16s, distance = "bray", add=T), scope = formula(cap.16s), direction = "forward", permutations = 1000)

# adjust p value
cap.16s.select.padj = cap.16s.select
cap.16s.select.padj$anova$`Pr(>F)` = p.adjust(cap.16s.select$anova$`Pr(>F)`, method = "BH")
cap.16s.select.padj$anova
cap.16s.select.padj

# by R2 
cap.16s.select.r2 = ordiR2step(capscale(otu.16s~1, data=env.16s, distance = "bray", add=T), scope = formula(cap.16s), direction = "forward", permutations = 1000, R2permutations = 1000)
cap.16s.select.r2 # aggrees with AIC method, otu.16s ~ cul_type + percPota + totalCrops

cap.16s.select.r2.padj = cap.16s.select.r2
cap.16s.select.r2.padj$anova$`Pr(>F)` = p.adjust(cap.16s.select.r2$anova$`Pr(>F)`, method = "BH")
cap.16s.select.r2.padj$anova

cap.16s.final = capscale(otu.16s ~ cul_type + percPota + totalCrops, data = env.16s, distance="bray", add=T)

cap.16s.final
saveRDS(cap.16s.final, "./field_soil_analysis_output/beta_diversity/cap_16_final.rds")

set.seed(104)
anova(cap.16s.final, permutations = 1000, step = 1000)
anova(cap.16s.final, by="terms", permutations = 1000, step = 1000)
anova(cap.16s.final, by="axis", permutations = 1000, step = 1000)

#------------------------------
### plot
#------------------------------

fdata16s.norm.genus = fdata16s.norm %>%
        tax_glom("Genus") %>%
        psmelt()

top_genus = fdata16s.norm.genus %>%
        select(OTU, Sample, Genus, Abundance) %>%
        droplevels.data.frame() %>%
        filter(Abundance > 0.01) %>%
        arrange(desc(Abundance)) %>%
        group_by(Sample) %>%
        top_n(5, wt=Abundance)

labelASV = top_genus[, c("OTU", "Genus")]
labelASV = labelASV[!duplicated(labelASV), ]
rownames(labelASV) = labelASV$OTU
labelASV = select(labelASV, Genus)

# bray curtis distance and NMDS
set.seed(333)
bc_nmds.16s = metaMDS(otu.16s, distance = "bray")

saveRDS(bc_nmds.16s, "./field_soil_analysis_output/beta_diversity/bc_nmds_16s.rds")

bc_nmds.16s = readRDS("./field_soil_analysis_output/beta_diversity/bc_nmds_16s.rds")
(stress.nmds.16s = bc_nmds.16s$stress)

bc_nmds.16s.sa = as.data.frame(scores(bc_nmds.16s, display = "sites"))
bc_nmds.16s.sa = bc_nmds.16s.sa %>% bind_cols(env.16s) # add metadata to sample scores
bc_nmds.16s.sp = as.data.frame(scores(bc_nmds.16s, display = "species"))


sp_label.16s = bc_nmds.16s.sp[rownames(bc_nmds.16s.sp) %in% rownames(labelASV), ]
sp_label_mg.16s = merge(sp_label.16s, labelASV, by=0)

# update 7/22/2020 do not show ASVs
g = ggplot()

(g.16s.nmds.cultype = g + 
                geom_point(data=bc_nmds.16s.sp, aes(x=NMDS1, y=NMDS2), color="grey", alpha=0.3) +
                # geom_point(data=sp_label_mg.16s, aes(x=NMDS1, y=NMDS2), color="darkblue", alpha=0.5) +
                # geom_text_repel(data=sp_label_mg.16s, aes(x=NMDS1, y=NMDS2, label=Genus), color="darkblue", size=3) +
                geom_point(data=bc_nmds.16s.sa, aes(x=NMDS1, y=NMDS2, color=cul_type), size=4, alpha=0.5) +
                stat_ellipse(data=bc_nmds.16s.sa, aes(x=NMDS1, y=NMDS2, color=cul_type), type = "norm") +
                stat_ellipse(data=bc_nmds.16s.sa, aes(x=NMDS1, y=NMDS2, color=cul_type), type = "t", linetype=2) +
                theme_bw() +
                geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                geom_label(aes(x=2.1, y=1.5, label=paste("Stress:", round(stress.nmds.16s, 2)))) +
                labs(color ="Cropping systems") +
                coord_fixed() +
                theme(plot.margin = unit(c(0.5,0.1,1,0), "cm")) +
                scale_color_discrete(labels = c("Conventional", "Mixed", "Organic"))
        
)


# cap plot

taxcols.16s = as.data.frame(tax_table(fdata16s.norm))
taxcols.16s$ASV = rownames(taxcols.16s)


(graph.cap.16s = plotCAP(cap.16s.final, env = env.16s, tax_cols = taxcols.16s, display = 3, splitDisplay = F, sa.xadj = 0.8, sa.yadj = -0.5,sa.colorSamples = T, sa.sample_group = "cul_type", sa.showCentroidGroup = "cul_type", sa.biplot = T,sa.ellipse = T, sa.addSp = T, sa.labelsp = T, sa.tax_level = "Genus") + scale_color_discrete(labels = c("Conventional", "Mixed", "Organic")) + labs(color = "Cropping systems")) 




#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# Fungal
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF



cap.its = capscale(otu.its ~cul_type+soil_compname+totalCrops+totalRotYrs+percPota+ph, data = env.its, distance = "bray", add = T)
cap.its

vif.cca(cap.its)


cap.its = capscale(otu.its ~cul_type+soil_compname+totalCrops+percPota+ph, data = env.its, distance = "bray", add=T)
cap.its

vif(cap.its)

set.seed(111)
anova(cap.its, permutations = 1000) # sig
anova(cap.its, by="terms", permutations = 1000) # ph not significant
anova(cap.its, by="axis", permutations = 1000) # top 3


# select by AIC
cap.its.select = ordistep(capscale(otu.its ~ 1, data=env.its, distance = "bray", add=T), scope = formula(cap.its), direction = "forward", permutations = 1000)

cap.its.select # percPota + cul_type + totalCrops + soil_compname
cap.its.select.padj = cap.its.select
cap.its.select.padj$anova$`Pr(>F)` = p.adjust(cap.its.select$anova$`Pr(>F)`, method = "BH")
cap.its.select.padj$anova # soil_compname borderline --- go use R2 selection instead



# select by R2
cap.its.select.r2 = ordiR2step(capscale(otu.its~1, data=env.its, distance = "bray", add=T), scope = formula(cap.its), direction="forward", permutations = 1000, R2permutations = 1000)
cap.its.select.r2 # otu.its ~ cul_type + percPota + totalCrop


cap.its.select.r2.padj = cap.its.select.r2
cap.its.select.r2.padj$anova$`Pr(>F)` = p.adjust(cap.its.select.r2$anova$`Pr(>F)`, method = "BH")
cap.its.select.r2.padj$anova

cap.its.final = capscale(otu.its ~ cul_type + percPota + totalCrops, data = env.its, distance = "bray", add=T)
cap.its.final

saveRDS(cap.its.final, "./field_soil_analysis_output/beta_diversity/cap_its_final.rds")

set.seed(10000)
anova(cap.its.final, permutations = 1000)
anova(cap.its.final, by="terms", permutations = 1000)
anova(cap.its.final, by="axis", permutations = 1000)

# plots

fdataITS.norm.genus = fdataITS.norm %>%
        tax_glom("Genus") %>%
        psmelt()



top_genus_its = fdataITS.norm.genus %>%
        select(OTU, Sample, Genus, Abundance) %>%
        droplevels.data.frame() %>%
        filter(Abundance > 0.01) %>%
        arrange(desc(Abundance)) %>%
        group_by(Sample) %>%
        top_n(5, wt=Abundance)

labelASV.its = top_genus_its[, c("OTU", "Genus")]
labelASV.its = labelASV.its[!duplicated(labelASV.its), ]
rownames(labelASV.its) = labelASV.its$OTU
labelASV.its = select(labelASV.its, Genus)


set.seed(335)
bc_nmds.its = metaMDS(otu.its, distance = "bray")
saveRDS(bc_nmds.its, "./field_soil_analysis_output/beta_diversity/bc_nmds_its.rds")
# read in
bc_nmds.its = readRDS("./field_soil_analysis_output/beta_diversity/bc_nmds_its.rds")
stressplot(bc_nmds.its)
(stress.nmds.its = bc_nmds.its$stress)


bc_nmds.its.sa = as.data.frame(scores(bc_nmds.its, display = "sites"))
bc_nmds.its.sa = bc_nmds.its.sa %>% bind_cols(env.its) # add metadata to sample scores
bc_nmds.its.sp = as.data.frame(scores(bc_nmds.its, display = "species"))


sp_label.its = bc_nmds.its.sp[rownames(bc_nmds.its.sp) %in% rownames(labelASV.its), ]
sp_label_mg.its = merge(sp_label.its, labelASV.its, by=0)



(g.its.nmds.cultype = g + 
                geom_point(data=bc_nmds.its.sp, aes(x=NMDS1, y=NMDS2), color="grey", alpha=0.3) +
                # geom_point(data=sp_label_mg.its, aes(x=NMDS1, y=NMDS2), color="darkblue", alpha=0.5) +
                # geom_text_repel(data=sp_label_mg.its, aes(x=NMDS1, y=NMDS2, label=Genus), color="darkblue", size=3) +
                geom_point(data=bc_nmds.its.sa, aes(x=NMDS1, y=NMDS2, color=cul_type), size=4, alpha=0.5) +
                stat_ellipse(data=bc_nmds.its.sa, aes(x=NMDS1, y=NMDS2, color=cul_type), type = "norm") +
                stat_ellipse(data=bc_nmds.its.sa, aes(x=NMDS1, y=NMDS2, color=cul_type), type = "t", linetype=2) +
                theme_bw() +
                geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                geom_label(aes(x=1.2, y=1.4, label=paste("Stress:", round(stress.nmds.its, 2)))) +
                labs(color ="Cropping systems") +
                coord_fixed() +
                theme(plot.margin = unit(c(0.5,0.1,1,0), "cm")) +
                scale_color_discrete(labels = c("Conventional", "Mixed", "Organic"))
        
)

# cap plot

taxcols.its = as.data.frame(tax_table(fdataITS.norm))
taxcols.its$ASV = rownames(taxcols.its)


(graph.cap.its = plotCAP(cap.its.final, env = env.its, tax_cols = taxcols.its, display = 3, splitDisplay = F, sa.xadj = 0.8, sa.yadj = -0.5,sa.colorSamples = T, sa.sample_group = "cul_type", sa.showCentroidGroup = "cul_type", sa.biplot = T,sa.ellipse = T, sa.addSp = T, sa.labelsp = T, sa.tax_level = "Genus")+ scale_color_discrete(labels = c("Conventional", "Mixed", "Organic")) + labs(color  = "Cropping systems"))




svg(file.path("./field_soil_analysis_output/beta_diversity/", "nmds_16s.svg"), width = 12, height = 13)
g.16s.nmds.cultype
dev.off()


svg(file.path("./field_soil_analysis_output/beta_diversity/", "nmds_its.svg"), width = 12, height = 13)
g.its.nmds.cultype
dev.off()


svg(file.path("./field_soil_analysis_output/beta_diversity/", "cap_16s.svg"), width = 12, height = 13)
graph.cap.16s
dev.off()


svg(file.path("./field_soil_analysis_output/beta_diversity/", "cap_its.svg"), width = 12, height = 13)
graph.cap.its
dev.off()




svg(file.path("./field_soil_analysis_output/beta_diversity/", "combined_cap.svg"), width = 17, height = 18)
ggarrange(g.16s.nmds.cultype, graph.cap.16s, g.its.nmds.cultype, graph.cap.its, labels = c("A", "B", "C", "D"), common.legend = T, align = "hv", legend = "bottom")
dev.off()


tiff(file.path("./field_soil_analysis_output/beta_diversity/", "combined_cap.tif"), width = 1000, compression="lzw")
ggarrange(g.16s.nmds.cultype, graph.cap.16s, g.its.nmds.cultype, graph.cap.its, labels = c("A", "B", "C", "D"), common.legend = T, align = "hv", legend = "bottom")
dev.off()

##############################################
# permanova
#############################################

bray.16s = vegdist(otu.16s, method = "bray")


# update remove totalRotYrs
set.seed(122)
perm.16s = adonis2(bray.16s ~ cul_type + soil_compname  + totalCrops + percPota + ph, data = env.16s, permutations = 10000, by="terms", add=T)
perm.16s


# check dispersion

disp.16s.cultype = betadisper(bray.16s, env.16s$cul_type)
set.seed(5)
permutest(disp.16s.cultype, permutations = 10000)

disp.16s.soiltype = betadisper(bray.16s, env.16s$soil_compname)
set.seed(5)
permutest(disp.16s.soiltype, permutations = 10000)


disp.16s.totalYrRot = betadisper(bray.16s, env.16s$totalRotYrs)
set.seed(5)
permutest(disp.16s.totalYrRot, permutations = 10000)


disp.16s.totalCrops = betadisper(bray.16s, env.16s$totalCrops)
set.seed(5)
permutest(disp.16s.totalCrops, permutations = 10000)

disp.16s.percPota = betadisper(bray.16s, env.16s$percPota)
set.seed(5)
permutest(disp.16s.percPota, permutations = 10000)



#--------------

bray.its = vegdist(otu.its, method = "bray")
set.seed(123)
perm.its = adonis2(bray.its ~ cul_type + soil_compname + totalCrops + percPota + ph, data = env.its, permutations = 10000, by="terms")
perm.its


disp.its.cultype = betadisper(bray.its, env.its$cul_type)
set.seed(5)
permutest(disp.its.cultype, permutations = 10000)

disp.its.soiltype = betadisper(bray.its, env.its$soil_compname)
set.seed(5)
permutest(disp.its.soiltype, permutations = 10000)


disp.its.totalYrRot = betadisper(bray.its, env.its$totalRotYrs)
set.seed(5)
permutest(disp.its.totalYrRot, permutations = 1000)


disp.its.totalCrops = betadisper(bray.its, env.its$totalCrops)
set.seed(5)
permutest(disp.its.totalCrops, permutations = 1000)

disp.its.percPota = betadisper(bray.its, env.its$percPota)
set.seed(5)
permutest(disp.its.percPota, permutations = 1000)



############################################################################

# taxonomy composition


rm(list = ls())

fdata_16s.prune = readRDS("./field_soil_analysis_output/data/fdata_16s_pruned.rds")
fdata_its.prune = readRDS("./field_soil_analysis_output/data/fdata_its_pruned.rds")

fdata16s.norm = fdata_16s.prune %>% transform_sample_counts(function(x)100*x/sum(x))
fdataITS.norm = fdata_its.prune %>% transform_sample_counts(function(x)100*x/sum(x))

source("./tools/aggAbundBar.R")
source("./tools/top10.R")

dir.create("./field_soil_analysis_output/taxonomy_composition")
#############################################################################



# summarise taxa

howmanytaxa = function(phylo, taxlevel){
        # exclude NA
        tb = tax_table(phylo)[, taxlevel]
        tb2 = tb[!is.na(tb)]
        return(length(unique(tb2)))
}

fdata16s.norm # 21915 taxa (ASVs)

howmanytaxa(fdata16s.norm, "Phylum")  # 40
howmanytaxa(fdata16s.norm, "Class")   # 127
howmanytaxa(fdata16s.norm, "Order")   # 346
howmanytaxa(fdata16s.norm, "Family")  # 579
howmanytaxa(fdata16s.norm, "Genus")   # 1144
howmanytaxa(fdata16s.norm, "Species") # 1457

fdataITS.norm # 3028 ASVs

howmanytaxa(fdataITS.norm, "Phylum")  # 16
howmanytaxa(fdataITS.norm, "Class")   # 46
howmanytaxa(fdataITS.norm, "Order")   # 103
howmanytaxa(fdataITS.norm, "Family")  # 213
howmanytaxa(fdataITS.norm, "Genus")   # 441
howmanytaxa(fdataITS.norm, "Species") # 668

cropping_systems = c(`Conventional`="C", `Green`="M", `Organic`="O")


(bar16s.phyl = aggAbund_barplot(fdata_16s.prune, "Phylum", pals::alphabet2(), legpos = "right", labeller = as_labeller(cropping_systems)))
ggsave("./field_soil_analysis_output/taxonomy_composition/bar16s_phyl.svg", device = "svg", width = 12, height = 10)

bar16s.genus = aggAbund_barplot(fdata_16s.prune, "Genus", pals::alphabet2(), labeller = as_labeller(cropping_systems))
ggsave("./field_soil_analysis_output/taxonomy_composition/bar16s_genus.png", device = "png", width = 12, height = 15)

(barITS.phyl = aggAbund_barplot(fdata_its.prune, "Phylum", glasbey(), legpos = "right", labeller = as_labeller(cropping_systems)))
ggsave("./field_soil_analysis_output/taxonomy_composition/barITS_phyl.svg", device = "svg", width = 12, height = 10)

barITS.genus = aggAbund_barplot(fdata_its.prune, "Genus", glasbey(), labeller = as_labeller(cropping_systems))
ggsave("./field_soil_analysis_output/taxonomy_composition/barITS_genus.png", device = "png", width = 12, height = 15)

######################
# top on average
#####################

phyl.16s = top_taxa(fdata_16s.prune, "Phylum")
phyl.its = top_taxa(fdata_its.prune, "Phylum")

genus.16s = top_taxa(fdata_16s.prune, "Genus")
genus.its = top_taxa(fdata_its.prune, "Genus")



# group by cultural practices
phyl.16s.cul = top_taxa(fdata_16s.prune, "Phylum", groupby = "cul_type")
genus.16s.cul = top_taxa(fdata_16s.prune, "Genus", groupby = "cul_type")


phyl.its.cul = top_taxa(fdata_its.prune, "Phylum", groupby = "cul_type")
genus.its.cul = top_taxa(fdata_its.prune, "Genus", groupby = "cul_type")


#############################################################################
# indicator species network

rm(list = ls())
dir.create("./field_soil_analysis_output/indicspecies_network/")


fdata_16s.prune = readRDS("./field_soil_analysis_output/data/fdata_16s_pruned.rds")
otu16s = as(t(otu_table(fdata_16s.prune)), "matrix")
meta16s = as(sample_data(fdata_16s.prune), "data.frame")


fdata_its.prune = readRDS("./field_soil_analysis_output/data/fdata_its_pruned.rds")
otuITS = as(t(otu_table(fdata_its.prune)), "matrix")
metaITS = as(sample_data(fdata_its.prune), "data.frame")
##############################################################################




# central log ratio normalize
norm_otu16s = clr(otu16s)
norm_otuITS = clr(otuITS)

cultypes16s = meta16s$cul_type
cultypesITS = metaITS$cul_type


length(unique(cultypes16s))
length(unique(cultypesITS))


# use multipatt, 10000 permutations
# set.seed(128)
# indi_sp_16s = indicspecies::multipatt(norm_otu16s, cultypes16s, func = "r.g", control = how(nperm=10000))
# indi_sp_its = indicspecies::multipatt(norm_otuITS, cultypesITS, func = "r.g", control = how(nperm=10000))

# saveRDS(indi_sp_16s, file.path("./field_soil_analysis_output/indicspecies_network/", "indi_sp_16s.rds"))
# saveRDS(indi_sp_its, file.path("./field_soil_analysis_output/indicspecies_network/", "indi_sp_its.rds"))
indi_sp_16s = readRDS("./field_soil_analysis_output/indicspecies_network/indi_sp_16s.rds")


summary(indi_sp_16s, alpha = 1, indvalcomp = T)
indi_sp_16s_df = indi_sp_16s$sign
head(indi_sp_16s_df)
tail(indi_sp_16s_df)

length((which(indi_sp_16s_df$s.Conventional == 1 & indi_sp_16s_df$s.Organic == 1 & indi_sp_16s_df$s.Green == 0 & indi_sp_16s_df$p.value < 0.05))) # 140

length(unique(rownames(indi_sp_16s_df))) # 21915 this is == the total bacterial ASVs

# get significant ASVs
C16s = as.matrix(indi_sp_16s_df[which(indi_sp_16s_df$s.Conventional == 1 & indi_sp_16s_df$p.value < 0.05), ]) 

C16s = as.data.frame(C16s)
C16s$ASVs = rownames(C16s)

length(unique(C16s$ASVs)) #163

#-----------------------
O16s = as.matrix(indi_sp_16s_df[which(indi_sp_16s_df$s.Organic == 1 & indi_sp_16s_df$p.value < 0.05), ])

O16s = as.data.frame(O16s)
O16s$ASVs = rownames(O16s)
length(unique(O16s$ASVs)) #195

#-----------------------
G16s = as.matrix(indi_sp_16s_df[which(indi_sp_16s_df$s.Green == 1 & indi_sp_16s_df$p.value < 0.05), ])

G16s = as.data.frame(G16s)
G16s$ASVs = rownames(G16s)
length(unique(G16s$ASVs)) #387

# combined together, NOTE: the rownames will change by default if two rownames are the same, an 1 will be added to the second one as XXX.1, 
sig_r_16s = rbind(C16s, O16s, G16s) # some of the rownames maybe shared
colnames(sig_r_16s)[1:3] = c("Conventional", "Green", "Organic") 
# range of correlation coefficients
range(sig_r_16s[, "stat"]) #0.3844348 0.9567103

# number of indicators
(indicators16s = sum(indi_sp_16s_df$p.value < 0.05, na.rm = T)) #590
# venn: 9 + 14 + 372 + 1 + 54 + 140 = 590



# precentage of 16s ASVs bound to cultural type
100*indicators16s/ nrow(t(norm_otu16s)) #2.69%
100*indicators16s/ 21915 #2.69%
nrow(t(norm_otu16s)) # 21915



# taxonomy composition in each cultural types

taxa_df_16s = as.data.frame(tax_table(fdata_16s.prune))

# indi_taxa_C.16s = taxa_df_16s[rownames(sig_r_16s[which(sig_r_16s[, "Conventional"] == 1), ]), ]
# indi_taxa_O.16s = taxa_df_16s[rownames(sig_r_16s[which(sig_r_16s[, "Organic"] == 1), ]), ]
# indi_taxa_G.16s = taxa_df_16s[rownames(sig_r_16s[which(sig_r_16s[, "Green"] == 1), ]),]
# 
# # how many ASVs in the venn for C
# 9+140+14 # 163
# length(unique(rownames(indi_taxa_C.16s))) #317
# 
# 
# indi_taxa_C.16s %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)
# 
# indi_taxa_O.16s %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)
# 
# indi_taxa_G.16s %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)
# 
# 
# 
# indi_taxa_C.16s %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)
# 
# indi_taxa_O.16s %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)
# 
# indi_taxa_G.16s %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)
# 
# 
# head(sig_r_16s)
# unique
#-----------------------------------------

# unique_indi_taxa_C.16s = taxa_df_16s[rownames(sig_r_16s[which((sig_r_16s[, "Green"] == 0 & sig_r_16s[, "Conventional"] == 1 & sig_r_16s[, "Organic"] == 0)), ]),]

unique_indi_taxa_C.16s = taxa_df_16s[sig_r_16s$ASVs[which((sig_r_16s[, "Green"] == 0 & sig_r_16s[, "Conventional"] == 1 & sig_r_16s[, "Organic"] == 0))],]

unique_indi_taxa_C.16s
length(unique(rownames(unique_indi_taxa_C.16s))) #9

#-----------------------------------------

# unique_indi_taxa_O.16s = taxa_df_16s[rownames(sig_r_16s[which((sig_r_16s[, "Green"] == 0 & sig_r_16s[, "Conventional"] == 0 & sig_r_16s[, "Organic"] == 1)), ]),]

unique_indi_taxa_O.16s = taxa_df_16s[sig_r_16s$ASVs[which((sig_r_16s[, "Green"] == 0 & sig_r_16s[, "Conventional"] == 0 & sig_r_16s[, "Organic"] == 1))],]

length(unique(rownames(unique_indi_taxa_O.16s))) #54

#--------------------------------------------
# unique_indi_taxa_G.16s = taxa_df_16s[rownames(sig_r_16s[which((sig_r_16s[, "Green"] == 1 & sig_r_16s[, "Conventional"] == 0 & sig_r_16s[, "Organic"] == 0)), ]),]
unique_indi_taxa_G.16s = taxa_df_16s[sig_r_16s$ASVs[which((sig_r_16s[, "Green"] == 1 & sig_r_16s[, "Conventional"] == 0 & sig_r_16s[, "Organic"] == 0))],]

length(unique(rownames(unique_indi_taxa_G.16s))) # 372

unique_indi_taxa_C.16s %>%
        group_by(Phylum) %>%
        tally() %>%
        mutate(freq = 100*n / sum(n)) %>%
        arrange(-freq)

unique_indi_taxa_C.16s %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100*n / sum(n)) %>%
        arrange(-freq)

length(rownames(unique_indi_taxa_O.16s)) #54

unique_indi_taxa_O.16s %>%
        group_by(Phylum) %>%
        tally() %>% 
        mutate(freq = 100*n / sum(n)) %>%
        arrange(-freq)


unique_indi_taxa_O.16s %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100*n / sum(n)) %>%
        arrange(-freq)

length(rownames(unique_indi_taxa_G.16s)) #372
unique_indi_taxa_G.16s %>%
        group_by(Phylum) %>%
        tally() %>%
        mutate(freq = 100*n / sum(n)) %>%
        arrange(-freq)

unique_indi_taxa_G.16s %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100*n / sum(n)) %>%
        arrange(-freq)      

# -------------------------------------
# fungi
rm(list = ls())

fdata_its.prune = readRDS("./field_soil_analysis_output/data/fdata_its_pruned.rds")
otuITS = as(t(otu_table(fdata_its.prune)), "matrix")
metaITS = as(sample_data(fdata_its.prune), "data.frame")
norm_otuITS = clr(otuITS)


indi_sp_its = readRDS("./field_soil_analysis_output/indicspecies_network/indi_sp_its.rds")
indi_sp_its_df = indi_sp_its$sign
head(indi_sp_its_df)

length((which(indi_sp_its_df$s.Conventional == 1 & indi_sp_its_df$s.Organic == 1 & indi_sp_its_df$s.Green == 0 & indi_sp_its_df$p.value < 0.05))) # 95

# get significant
Cits = as.matrix(indi_sp_its_df[which(indi_sp_its_df$s.Conventional == 1 & indi_sp_its_df$p.value < 0.05), ])

Cits = as.data.frame(Cits)
Cits$ASVs = rownames(Cits)
head(Cits)
length(unique(rownames(Cits))) #102

Oits = as.matrix(indi_sp_its_df[which(indi_sp_its_df$s.Organic == 1 & indi_sp_its_df$p.value < 0.05), ])

Oits = as.data.frame(Oits)
Oits$ASVs = rownames(Oits)
length(unique(rownames(Oits))) #129

Gits = as.matrix(indi_sp_its_df[which(indi_sp_its_df$s.Green == 1 & indi_sp_its_df$p.value < 0.05), ])

Gits = as.data.frame(Gits)
Gits$ASVs = rownames(Gits)
length(unique(rownames(Gits))) #89

# combined all cultural types, some rownames were repeated thus suffix by 1: XXX.1
sig_r_its = rbind(Cits, Oits, Gits)

head(sig_r_its)
colnames(sig_r_its)[1:3] = c("Conventional", "Green", "Organic")



# range of correlation coefficients
range(sig_r_its[, "stat"]) #0.4191476 0.9460002


# number of indicators
1+6+79+4+30+95 #215
indicatorsITS = sum(indi_sp_its_df$p.value < 0.05, na.rm = T)
indicatorsITS # 215
nrow(t(norm_otuITS))

length(unique(colnames(norm_otuITS))) #3028

# precentage of ITS ASVs bound to cultural type

100*indicatorsITS / nrow(t(norm_otuITS)) # 7.1%

taxa_df_its = as.data.frame(tax_table(fdata_its.prune))


# indi_taxa_C.its = taxa_df_its[rownames(sig_r_its[which(sig_r_its[, "Conventional"] == 1), ]), ]
# indi_taxa_O.its = taxa_df_its[rownames(sig_r_its[which(sig_r_its[, "Organic"] == 1), ]), ]
# indi_taxa_G.its = taxa_df_its[rownames(sig_r_its[which(sig_r_its[, "Green"] == 1), ]), ]
#---
# indi_taxa_C.its = taxa_df_its[sig_r_its$ASVs[which(sig_r_its[, "Conventional"] == 1)],]
# 
# length(unique(sig_r_its$ASVs[which(sig_r_its[, "Conventional"] == 1)])) #102
# taxa_df_its[sig_r_its$ASVs[which(sig_r_its[, "Conventional"] == 1)],]
# length(rownames(indi_taxa_C.its)) #203
# 1+6+95 # 102
# 
# indi_taxa_C.its %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)
# 
# length(rownames(indi_taxa_O.its)) #228
# indi_taxa_O.its %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)
# 
# length(rownames(indi_taxa_G.its)) #99
# indi_taxa_G.its %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)
# 
# 
# 
# indi_taxa_C.its %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)
# 
# indi_taxa_O.its %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)
# 
# indi_taxa_G.its %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100*n / sum(n)) %>%
#         arrange(-freq)


# unique
# unique_indi_taxa_C.its = taxa_df_its[rownames(sig_r_its[which((sig_r_its[, "Green"] == 0 & sig_r_its[, "Conventional"] == 1 & sig_r_its[, "Organic"] == 0)), ]),] # this was not detected
unique_indi_taxa_C.its = taxa_df_its[sig_r_its$ASVs[which((sig_r_its[, "Green"] == 0 & sig_r_its[, "Conventional"] == 1 & sig_r_its[, "Organic"] == 0))],]


# unique_indi_taxa_O.its = taxa_df_its[rownames(sig_r_its[which((sig_r_its[, "Green"] == 0 & sig_r_its[, "Conventional"] == 0 & sig_r_its[, "Organic"] == 1)), ]),]

unique_indi_taxa_O.its = taxa_df_its[sig_r_its$ASVs[which((sig_r_its[, "Green"] == 0 & sig_r_its[, "Conventional"] == 0 & sig_r_its[, "Organic"] == 1))],]

length(rownames(unique_indi_taxa_O.its)) #30

# unique_indi_taxa_G.its = taxa_df_its[rownames(sig_r_its[which((sig_r_its[, "Green"] == 1 & sig_r_its[, "Conventional"] == 0 & sig_r_its[, "Organic"] == 0)), ]),]
unique_indi_taxa_G.its = taxa_df_its[sig_r_its$ASVs[which((sig_r_its[, "Green"] == 1 & sig_r_its[, "Conventional"] == 0 & sig_r_its[, "Organic"] == 0))],]

length(rownames(unique_indi_taxa_G.its)) #79


unique_indi_taxa_C.its %>%
        group_by(Phylum) %>%
        tally() %>%
        mutate(freq = 100*n / sum(n)) %>%
        arrange(-freq)

unique_indi_taxa_O.its %>%
        group_by(Phylum) %>%
        tally() %>% 
        mutate(freq = 100*n / sum(n)) %>%
        arrange(-freq)

unique_indi_taxa_O.its %>%
        group_by(Genus) %>%
        tally() %>% 
        mutate(freq = 100*n / sum(n)) %>%
        arrange(-freq)

unique_indi_taxa_G.its %>%
        group_by(Phylum) %>%
        tally() %>%
        mutate(freq = 100*n / sum(n)) %>%
        arrange(-freq)

unique_indi_taxa_G.its %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100*n / sum(n)) %>%
        arrange(-freq)


#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# bipartite network for bacterial
source("./tools/assignColors.R")
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

# node table
# bipar16s = data.frame(from = c(rep("Conventional", length(which(sig_r_16s[, "Conventional"]==1))),
#                                rep("Organic", length(which(sig_r_16s[, "Organic"]==1))),
#                                rep("Green", length(which(sig_r_16s[, "Green"]==1)))),
#                       
#                       to = c(rownames(sig_r_16s)[which(sig_r_16s[, "Conventional"]==1)],
#                              rownames(sig_r_16s)[which(sig_r_16s[, "Organic"]==1)],
#                              rownames(sig_r_16s)[which(sig_r_16s[, "Green"]==1)]),
#                       
#                       r = c(sig_r_16s[which(sig_r_16s[, "Conventional"]==1), "stat"],
#                             sig_r_16s[which(sig_r_16s[, "Organic"]==1), "stat"],
#                             sig_r_16s[which(sig_r_16s[, "Green"]==1), "stat"]
#                       )
# )
bipar16s = data.frame(from = c(rep("Conventional", length(which(sig_r_16s[, "Conventional"]==1))),
                               rep("Organic", length(which(sig_r_16s[, "Organic"]==1))),
                               rep("Green", length(which(sig_r_16s[, "Green"]==1)))),
                      
                      to = c(sig_r_16s$ASVs[which(sig_r_16s[, "Conventional"]==1)],
                             sig_r_16s$ASVs[which(sig_r_16s[, "Organic"]==1)],
                             sig_r_16s$ASVs[which(sig_r_16s[, "Green"]==1)]),
                      
                      r = c(sig_r_16s[which(sig_r_16s[, "Conventional"]==1), "stat"],
                            sig_r_16s[which(sig_r_16s[, "Organic"]==1), "stat"],
                            sig_r_16s[which(sig_r_16s[, "Green"]==1), "stat"]
                      )
)

# node attributes
# bipar16s_vattr = data.frame(node=unique(rownames(sig_r_16s)), indicgroup=0)
bipar16s_vattr = data.frame(node=unique(sig_r_16s$ASVs), indicgroup=0)

for(i in as.character(bipar16s_vattr$node)){
        bipar16s_vattr[bipar16s_vattr$node==i, "indicgroup"] <- paste(
                colnames(sig_r_16s)[which(sig_r_16s[i, 1:3]==1)], 
                collapse = "_")
}

# taxa_df_16s = as.data.frame(tax_table(fdata_16s.prune))
bipar16s_vattr = cbind(bipar16s_vattr, taxa_df_16s[as.character(bipar16s_vattr$node), ])


head(bipar16s)
head(bipar16s_vattr)

unique(bipar16s_vattr$indicgroup)


# create bipartite network with igraph
binet16s = graph.data.frame(bipar16s, directed = F)
# V(binet16s)$type = V(binet16s)$name %in% bipar16s[,1]
binet16s = simplify(binet16s, remove.multiple = T, remove.loops = T)

# set label for nodes, remove label
V(binet16s)$label = V(binet16s)$name
V(binet16s)$label = gsub("bASV.*", "", V(binet16s)$label)
V(binet16s)$annotate = gsub("bASV.*", "", V(binet16s)$label)


#set node size
V(binet16s)$size = c(rep(5, 3), rep(2, length(V(binet16s)) - 3))

# set node shapes
V(binet16s)$shape = c(rep("square", 3), rep("circle", length(V(binet16s))-3))

# node colors

taxa_df_16s_2 = assignColors(taxa_df_16s, "Phylum")

V(binet16s)$color = V(binet16s)$name


V(binet16s)$color = taxa_df_16s_2[V(binet16s)$name, ]$colrs
V(binet16s)$color[1:3] = "white"
V(binet16s)$frame.color = V(binet16s)$color

V(binet16s)$Phylum = as.character(taxa_df_16s_2[V(binet16s)$name, ]$Phylum)
V(binet16s)$Phylum[1:3] = V(binet16s)$name[1:3]
V(binet16s)$Genus = as.character(taxa_df_16s_2[V(binet16s)$name, ]$Genus)
V(binet16s)$Genus[1:3] = V(binet16s)$name[1:3]

head(V(binet16s))

write_graph(binet16s, file = file.path("./field_soil_analysis_output/indicspecies_network/", "binet16s.gml"), format = "gml")


#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# bipartite network for ITS
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF



# node table
# from cultype to ASVs
# biparits = data.frame(from = c(rep("Conventional", length(which(sig_r_its[, "Conventional"]==1))),
#                                rep("Organic", length(which(sig_r_its[, "Organic"]==1))),
#                                rep("Green", length(which(sig_r_its[, "Green"]==1)))),
#                       
#                       to = c(rownames(sig_r_its)[which(sig_r_its[, "Conventional"]==1)],
#                              rownames(sig_r_its)[which(sig_r_its[, "Organic"]==1)],
#                              rownames(sig_r_its)[which(sig_r_its[, "Green"]==1)]),
#                       
#                       r = c(sig_r_its[which(sig_r_its[, "Conventional"]==1), "stat"],
#                             sig_r_its[which(sig_r_its[, "Organic"]==1), "stat"],
#                             sig_r_its[which(sig_r_its[, "Green"]==1), "stat"]
#                       )
# )

biparits = data.frame(from = c(rep("Conventional", length(which(sig_r_its[, "Conventional"]==1))),
                               rep("Organic", length(which(sig_r_its[, "Organic"]==1))),
                               rep("Green", length(which(sig_r_its[, "Green"]==1)))),
                      
                      to = c(sig_r_its$ASVs[which(sig_r_its[, "Conventional"]==1)],
                             sig_r_its$ASVs[which(sig_r_its[, "Organic"]==1)],
                             sig_r_its$ASVs[which(sig_r_its[, "Green"]==1)]),
                      
                      r = c(sig_r_its[which(sig_r_its[, "Conventional"]==1), "stat"],
                            sig_r_its[which(sig_r_its[, "Organic"]==1), "stat"],
                            sig_r_its[which(sig_r_its[, "Green"]==1), "stat"]
                      )
)

# node attributes
# biparits_vattr = data.frame(node=unique(rownames(sig_r_its)), indicgroup=0)
biparits_vattr = data.frame(node=unique(sig_r_its$ASVs), indicgroup=0)

for(i in as.character(biparits_vattr$node)){
        biparits_vattr[biparits_vattr$node==i, "indicgroup"] <- paste(
                colnames(sig_r_its)[which(sig_r_its[i, 1:3]==1)], 
                collapse = "_")
}

taxa_df_its = as.data.frame(as(tax_table(fdata_its.prune), "matrix"))
biparits_vattr = cbind(biparits_vattr, taxa_df_its[as.character(biparits_vattr$node), ])

head(biparits_vattr)

# create bipartite network with igraph
binetits = graph.data.frame(biparits, directed = F)
# V(binetits)$type = V(binetits)$name %in% biparits[,1]
binetits = simplify(binetits, remove.multiple = T, remove.loops = T)

# set label for nodes, remove label
V(binetits)$label = V(binetits)$name
V(binetits)$label = gsub("fASV.*", "", V(binetits)$label)

# set node size, the 3 hub - size 5, 
V(binetits)$size = c(rep(5, 3), rep(2, length(V(binetits)) - 3))
# set node shapes
V(binetits)$shape = c(rep("square", 3), rep("circle", length(V(binetits))-3))


# node colors by phylum

taxa_df_its_2 = assignColors(taxa_df_its, "Phylum")


V(binetits)$color = V(binetits)$name

V(binetits)$color = taxa_df_its_2[V(binetits)$name, ]$colrs
V(binetits)$color[1:3] = "white"
V(binetits)$frame.color = V(binetits)$color

V(binetits)$Phylum = as.character(taxa_df_its_2[V(binetits)$name, ]$Phylum)
V(binetits)$Phylum[1:3] = V(binetits)$name[1:3]
V(binetits)$Genus = as.character(taxa_df_its_2[V(binetits)$name, ]$Genus)
V(binetits)$Genus[1:3] = V(binetits)$name[1:3]

write_graph(binetits, file = file.path("./field_soil_analysis_output/indicspecies_network/", "binetits.gml"), format = "gml")



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# description of the network
source("./tools/summIndic.R")
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



conven16s   = IndicTaxaSumm(bipar16s_vattr, "Conventional")
conven16s_g = IndicTaxaSumm(bipar16s_vattr, grp = "Conventional", taxlevel = "Genus")


organic16s   = IndicTaxaSumm(bipar16s_vattr, "Organic")
organic16s_g = IndicTaxaSumm(bipar16s_vattr, grp = "Organic", taxlevel = "Genus")


green16s = IndicTaxaSumm(bipar16s_vattr, "Green")
green16s_g = IndicTaxaSumm(bipar16s_vattr, grp = "Green", taxlevel = "Genus")


conven_org16s = IndicTaxaSumm(bipar16s_vattr, "Conventional_Organic")
conven_org16s_g = IndicTaxaSumm(bipar16s_vattr, grp = "Conventional_Organic", taxlevel = "Genus")


conven_green16s = IndicTaxaSumm(bipar16s_vattr, "Conventional_Green")
conven_green16s_g = IndicTaxaSumm(bipar16s_vattr, grp = "Conventional_Green", taxlevel = "Genus")


green_org16s = IndicTaxaSumm(bipar16s_vattr, "Green_Organic")
green_org16s_g = IndicTaxaSumm(bipar16s_vattr, grp = "Green_Organic", taxlevel = "Genus")



# fungal

convenITS   = IndicTaxaSumm(biparits_vattr, "Conventional")
convenITS_g = IndicTaxaSumm(biparits_vattr, grp = "Conventional", taxlevel = "Genus")

organicITS   = IndicTaxaSumm(biparits_vattr, "Organic")
organicITS_g   = IndicTaxaSumm(biparits_vattr, grp = "Organic", taxlevel = "Genus")

greenITS   = IndicTaxaSumm(biparits_vattr, "Green")
greenITS_g   = IndicTaxaSumm(biparits_vattr, grp = "Green", taxlevel = "Genus")

conven_orgITS   = IndicTaxaSumm(biparits_vattr, "Conventional_Organic")
conven_orgITS_g   = IndicTaxaSumm(biparits_vattr, grp = "Conventional_Organic", taxlevel = "Genus")

conven_greenITS   = IndicTaxaSumm(biparits_vattr, "Conventional_Green")
conven_greenITS_g   = IndicTaxaSumm(biparits_vattr, grp = "Conventional_Green", taxlevel = "Genus")

green_orgITS   = IndicTaxaSumm(biparits_vattr, "Green_Organic")
green_orgITS_g   = IndicTaxaSumm(biparits_vattr, grp = "Green_Organic", taxlevel = "Genus")



#############################
# venn diagram
############################

library(venn)


# 16S



# including all possible

conv16s_asv = rownames(bipar16s_vattr[bipar16s_vattr$indicgroup == "Conventional" | bipar16s_vattr$indicgroup == "Conventional_Organic" | bipar16s_vattr$indicgroup == "Conventional_Green",])

org16s_asv = rownames(bipar16s_vattr[bipar16s_vattr$indicgroup == "Organic" | bipar16s_vattr$indicgroup == "Conventional_Organic" | bipar16s_vattr$indicgroup == "Green_Organic",])

green16s_asv = rownames(bipar16s_vattr[bipar16s_vattr$indicgroup == "Green" | bipar16s_vattr$indicgroup == "Green_Organic" | bipar16s_vattr$indicgroup == "Conventional_Green",])


conv16s_asv_p = as.character(bipar16s_vattr$Phylum[bipar16s_vattr$indicgroup == "Conventional" | bipar16s_vattr$indicgroup == "Conventional_Organic" | bipar16s_vattr$indicgroup == "Conventional_Green"])

org16s_asv_p = as.character(bipar16s_vattr$Phylum[bipar16s_vattr$indicgroup == "Organic" | bipar16s_vattr$indicgroup == "Conventional_Organic" | bipar16s_vattr$indicgroup == "Green_Organic"])

green16s_asv_p = as.character(bipar16s_vattr$Phylum[bipar16s_vattr$indicgroup == "Green" | bipar16s_vattr$indicgroup == "Green_Organic" | bipar16s_vattr$indicgroup == "Conventional_Green"])




conv16s_asv_g = as.character(bipar16s_vattr$Genus[bipar16s_vattr$indicgroup == "Conventional" | bipar16s_vattr$indicgroup == "Conventional_Organic" | bipar16s_vattr$indicgroup == "Conventional_Green"])

org16s_asv_g = as.character(bipar16s_vattr$Genus[bipar16s_vattr$indicgroup == "Organic" | bipar16s_vattr$indicgroup == "Conventional_Organic" | bipar16s_vattr$indicgroup == "Green_Organic"])

green16s_asv_g = as.character(bipar16s_vattr$Genus[bipar16s_vattr$indicgroup == "Green" | bipar16s_vattr$indicgroup == "Green_Organic" | bipar16s_vattr$indicgroup == "Conventional_Green"])

#FFFFFFFFFFFFFFFFFFFFFFF
# Fungi


convITS_asv = rownames(biparits_vattr[biparits_vattr$indicgroup == "Conventional" | biparits_vattr$indicgroup == "Conventional_Organic" | biparits_vattr$indicgroup == "Conventional_Green",])
orgITS_asv = rownames(biparits_vattr[biparits_vattr$indicgroup == "Organic" | biparits_vattr$indicgroup == "Conventional_Organic" | biparits_vattr$indicgroup == "Green_Organic",])
greenITS_asv = rownames(biparits_vattr[biparits_vattr$indicgroup == "Green" | biparits_vattr$indicgroup == "Green_Organic" | biparits_vattr$indicgroup == "Conventional_Green",])




convITS_asv_p = as.character(biparits_vattr$Phylum[biparits_vattr$indicgroup == "Conventional" | biparits_vattr$indicgroup == "Conventional_Organic" | biparits_vattr$indicgroup == "Conventional_Green"])

orgITS_asv_p = as.character(biparits_vattr$Phylum[biparits_vattr$indicgroup == "Organic" | biparits_vattr$indicgroup == "Conventional_Organic" | biparits_vattr$indicgroup == "Green_Organic"])

greenITS_asv_p = as.character(biparits_vattr$Phylum[biparits_vattr$indicgroup == "Green" | biparits_vattr$indicgroup == "Green_Organic" | biparits_vattr$indicgroup == "Conventional_Green"])

convITS_asv_g = as.character(biparits_vattr$Genus[biparits_vattr$indicgroup == "Conventional" | biparits_vattr$indicgroup == "Conventional_Organic" | biparits_vattr$indicgroup == "Conventional_Green"])

orgITS_asv_g = as.character(biparits_vattr$Genus[biparits_vattr$indicgroup == "Organic" | biparits_vattr$indicgroup == "Conventional_Organic" | biparits_vattr$indicgroup == "Green_Organic"])

greenITS_asv_g = as.character(biparits_vattr$Genus[biparits_vattr$indicgroup == "Green" | biparits_vattr$indicgroup == "Green_Organic" | biparits_vattr$indicgroup == "Conventional_Green"])


unique(conv16s_asv_p) #13
unique(org16s_asv_p) #14
unique(green16s_asv_p) #15




venn_colors = c("cyan3", "khaki1", "hotpink3")
venn_colors_2 = c("forestgreen", "indianred1", "royalblue1")

svg("./field_soil_analysis_output/indicspecies_network/16S_asv.svg")
venn::venn(list(conventional = conv16s_asv, organic = org16s_asv, mixed = green16s_asv), zcolor = venn_colors_2, box=F, borders = F, plotsize = 20)
dev.off()

svg("./field_soil_analysis_output/indicspecies_network/16S_phylum.svg")
venn::venn(list(conventional = conv16s_asv_p, organic = org16s_asv_p, mixed = green16s_asv_p), zcolor = venn_colors_2, box=F, borders = F, plotsize = 20)
dev.off()

svg("./field_soil_analysis_output/indicspecies_network/16S_genus.svg")
venn::venn(list(conventional = conv16s_asv_g, organic = org16s_asv_g, mixed = green16s_asv_g), zcolor = venn_colors_2, box=F, borders = F, plotsize = 20)
dev.off()





svg("./field_soil_analysis_output/indicspecies_network/ITS2_asv.svg")
venn::venn(list(conventional = convITS_asv, organic = orgITS_asv, mixed = greenITS_asv), zcolor = venn_colors_2, box=F, borders = F, plotsize = 20)
dev.off()

svg("./field_soil_analysis_output/indicspecies_network/ITS2_phylum.svg")
venn::venn(list(conventional = convITS_asv_p, organic = orgITS_asv_p, mixed = greenITS_asv_p), zcolor = venn_colors_2, box=F, borders = F, plotsize = 20)
dev.off()

svg("./field_soil_analysis_output/indicspecies_network/ITS2_genus.svg")
venn::venn(list(conventional = convITS_asv_g, organic = orgITS_asv_g, mixed = greenITS_asv_g), zcolor = venn_colors_2, box=F, borders = F, plotsize = 20)
dev.off()


##pie chart
source("./tools/summIndic.R")


(b1 = indic_pieChart(conv16s_asv_p, n=0.02, inc = 0.1, x_baselabel = 1.65))
(b2 = indic_pieChart(org16s_asv_p, n=0.02, inc = 0.1, x_baselabel = 1.65))
(b3 = indic_pieChart(green16s_asv_p, n=0.02, inc = 0.1, x_baselabel = 1.65))


(f1 = indic_pieChart(convITS_asv_p, n=0.02, inc = 0.1, x_baselabel = 1.65))
(f2 = indic_pieChart(orgITS_asv_p, n=0.02, inc = 0.1, x_baselabel = 1.65))
(f3 = indic_pieChart(greenITS_asv_p, n=0.02, inc = 0.1, x_baselabel = 1.65))



#(((((((((((((((((((((((())))))))))))))))))))))))

# Phylum shared

# unique to conventional
setdiff(setdiff(conv16s_asv_p, org16s_asv_p), green16s_asv_p)
# unique to organic
p_uo16s = setdiff(setdiff(org16s_asv_p, conv16s_asv_p), green16s_asv_p)
p_uo16s
# unique to green
p_ug16s = setdiff(setdiff(green16s_asv_p, conv16s_asv_p), org16s_asv_p)
write.table(p_ug16s, "./field_soil_analysis_output/indicspecies_network/venn_list/p_ug16s.txt", row.names = F, quote = F, col.names = F)


# unique to conventional and organic
p_uco16s = setdiff(intersect(conv16s_asv_p, org16s_asv_p), green16s_asv_p)
p_uco16s
write.table(p_uco16s, "./field_soil_analysis_output/indicspecies_network/venn_list/p_uco16s.txt", row.names = F, quote = F, col.names = F)


# unique to organic and green
p_ugo16s = setdiff(intersect(org16s_asv_p,green16s_asv_p), conv16s_asv_p)
write.table(p_ugo16s, "./field_soil_analysis_output/indicspecies_network/venn_list/p_ugo16s.txt", row.names = F, quote = F, col.names = F)


# shared by all three
p_uall16s = Reduce(intersect, list(conv16s_asv_p, org16s_asv_p, green16s_asv_p))
write.table(p_uall16s, "./field_soil_analysis_output/indicspecies_network/venn_list/p_uall16s.txt", row.names = F, quote = F, col.names = F)


#(((((((((((((((((((((())))))))))))))))))))))
# Genus shared
# unique to conventional
setdiff(setdiff(conv16s_asv_g, org16s_asv_g), green16s_asv_g)

# unique to organic
uo16s = setdiff(setdiff(org16s_asv_g, conv16s_asv_g), green16s_asv_g)
write.table(uo16s, "./field_soil_analysis_output/indicspecies_network/venn_list/uo16s.txt", row.names = F, quote = F, col.names = F)

# unique to green
ug16s = setdiff(setdiff(green16s_asv_g, conv16s_asv_g), org16s_asv_g)
write.table(ug16s, "./field_soil_analysis_output/indicspecies_network/venn_list/ug16s.txt", row.names = F, quote = F, col.names = F)


# unique to conventional and organic
uco16s = setdiff(intersect(conv16s_asv_g, org16s_asv_g), green16s_asv_g)
write.table(uco16s, "./field_soil_analysis_output/indicspecies_network/venn_list/uco16s.txt", row.names = F, quote = F, col.names = F)

# unique to conventional and green
ucg16s = setdiff(intersect(conv16s_asv_g, green16s_asv_g), org16s_asv_g)
write.table(ucg16s, "./field_soil_analysis_output/indicspecies_network/venn_list/ucg16s.txt", row.names = F, quote = F, col.names = F)

# unique to green and organic
ugo16s = setdiff(intersect(green16s_asv_g, org16s_asv_g), conv16s_asv_g)
write.table(ugo16s, "./field_soil_analysis_output/indicspecies_network/venn_list/ugo16s.txt", row.names = F, quote = F, col.names = F)

# shared by 3
uall16s = Reduce(intersect, list(conv16s_asv_g, org16s_asv_g, green16s_asv_g))
write.table(uall16s, "./field_soil_analysis_output/indicspecies_network/venn_list/uall16s.txt", row.names = F, quote = F, col.names = F)


#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

#(((((((((((((((((((((((())))))))))))))))))))))))

# Phylum shared

# unique to conventional
setdiff(setdiff(convITS_asv_p, orgITS_asv_p), greenITS_asv_p)
# unique to organic
p_uoits = setdiff(setdiff(orgITS_asv_p, convITS_asv_p), greenITS_asv_p)
p_uoits
# unique to green
p_ugits = setdiff(setdiff(greenITS_asv_p, convITS_asv_p), orgITS_asv_p)
p_ugits
# None


# unique to conventional and organic
# None


# unique to organic and green
# None



# shared by all three
p_uallits = Reduce(intersect, list(convITS_asv_p, orgITS_asv_p, greenITS_asv_p))
write.table(p_uallits, "./field_soil_analysis_output/indicspecies_network/venn_list/p_uallits.txt", row.names = F, quote = F, col.names = F)

#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
#(((((((((((((((((((((())))))))))))))))))))))
# Genus shared
# unique to conventional
setdiff(setdiff(convITS_asv_g, orgITS_asv_g), greenITS_asv_g)
# Geosmithia

# unique to organic
uoits = setdiff(setdiff(orgITS_asv_g, convITS_asv_g), greenITS_asv_g)
write.table(uoits, "./field_soil_analysis_output/indicspecies_network/venn_list/uoits.txt", row.names = F, quote = F, col.names = F)

# unique to green
ugits = setdiff(setdiff(greenITS_asv_g, convITS_asv_g), orgITS_asv_g)
write.table(ugits, "./field_soil_analysis_output/indicspecies_network/venn_list/ugits.txt", row.names = F, quote = F, col.names = F)


# unique to conventional and organic
ucoits = setdiff(intersect(convITS_asv_g, orgITS_asv_g), greenITS_asv_g)
write.table(ucoits, "./field_soil_analysis_output/indicspecies_network/venn_list/ucoits.txt", row.names = F, quote = F, col.names = F)

# unique to conventional and green
ucgits = setdiff(intersect(convITS_asv_g, greenITS_asv_g), orgITS_asv_g)
write.table(ucgits, "./field_soil_analysis_output/indicspecies_network/venn_list/ucgits.txt", row.names = F, quote = F, col.names = F)

# unique to green and organic
ugoits = setdiff(intersect(greenITS_asv_g, orgITS_asv_g), convITS_asv_g)
write.table(ugoits, "./field_soil_analysis_output/indicspecies_network/venn_list/ugoits.txt", row.names = F, quote = F, col.names = F)

# shared by 3
uallits = Reduce(intersect, list(convITS_asv_g, orgITS_asv_g, greenITS_asv_g))
write.table(uallits, "./field_soil_analysis_output/indicspecies_network/venn_list/uallits.txt", row.names = F, quote = F, col.names = F)


# phylum shared by 3 managements
Reduce(intersect, list(convITS_asv_p, orgITS_asv_p, greenITS_asv_p))

#####################################################################################

# DESeq2

rm(list = ls())
source("./tools/plotDeseq2.R")


# dir.create("./field_soil_analysis_output/Deseq2")
# dir.create("./field_soil_analysis_output/Deseq2/sigout/16s", recursive = T)
# dir.create("./field_soil_analysis_output/Deseq2/sigout/its", recursive = T)

fdata_16s.prune = readRDS("./field_soil_analysis_output/data/fdata_16s_pruned.rds")
fdata_its.prune = readRDS("./field_soil_analysis_output/data/fdata_its_pruned.rds")


env.16s = data.frame(sample_data(fdata_16s.prune), stringsAsFactors = F)
env.its = data.frame(sample_data(fdata_its.prune), stringsAsFactors = F)



f16s.dds = phyloseq_to_deseq2(fdata_16s.prune, ~ cul_type)
fits.dds = phyloseq_to_deseq2(fdata_its.prune, ~ cul_type)


# geomean
gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
#####################################################################################


geoMean.16s = apply(counts(f16s.dds), 1, gm_mean)
f16s.dds.mean = estimateSizeFactors(f16s.dds, geoMean = geoMean.16s)
f16s.dds.mean = DESeq(f16s.dds.mean, fitType = "local")


geoMean.its = apply(counts(fits.dds), 1, gm_mean)
fits.dds.mean = estimateSizeFactors(fits.dds, geoMean = geoMean.its)
fits.dds.mean = DESeq(fits.dds.mean, fitType = "local")


png("./field_soil_analysis_output/Deseq2/disper_f16s_dds.png")
plotDispEsts(f16s.dds.mean)
dev.off()


png("./field_soil_analysis_output/Deseq2/disper_fits_dds.png")
plotDispEsts(fits.dds.mean)
dev.off()


#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# bacteria
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# conventional vs organic
ddsout_c_o = plotDeseq2(f16s.dds.mean, fdata_16s.prune, c("cul_type", "Conventional", "Organic"), tax_level = "Genus", alpha = 0.01, plottype = "hbar", show_n = 10, summar = F, outputpath = "./field_soil_analysis_output/Deseq2/sigout/16s/")

# conventional vs green
ddsout_c_g = plotDeseq2(f16s.dds.mean, fdata_16s.prune, c("cul_type", "Conventional", "Green"), tax_level = "Genus", alpha = 0.01, plottype = "hbar", show_n = 10, summar = F, outputpath = "./field_soil_analysis_output/Deseq2/sigout/16s/", legendPos = NULL)

# organic vs green
ddsout_o_g = plotDeseq2(f16s.dds.mean, fdata_16s.prune, c("cul_type", "Organic", "Green"), tax_level = "Genus", alpha = 0.01, plottype = "hbar", show_n = 10, summar = F, outputpath = "./field_soil_analysis_output/Deseq2/sigout/16s/")

#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# Fungi
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

# conventional vs organic
ddsout_c_o_its = plotDeseq2(fits.dds.mean, fdata_its.prune, c("cul_type", "Conventional", "Organic"), tax_level = "Genus", alpha = 0.01, plottype = "hbar", show_n = 10, summar = F, outputpath = "./field_soil_analysis_output/Deseq2/sigout/its/", legendPos = NULL)

# conventional vs green
ddsout_c_g_its = plotDeseq2(fits.dds.mean, fdata_its.prune, c("cul_type", "Conventional", "Green"), tax_level = "Genus", alpha = 0.01, plottype = "hbar", show_n = 10, summar = F, outputpath = "./field_soil_analysis_output/Deseq2/sigout/its/", legendPos = NULL)

# organic vs green
ddsout_o_g_its = plotDeseq2(fits.dds.mean, fdata_its.prune, c("cul_type", "Organic", "Green"), tax_level = "Genus", alpha = 0.01, plottype = "hbar", show_n = 10, summar = F, outputpath = "./field_soil_analysis_output/Deseq2/sigout/its/")

#############
#output plots
#############
svg("./field_soil_analysis_output/Deseq2/c_vs_o.svg", width = 15, height = 16)
ddsout_c_o$graph
dev.off()


svg("./field_soil_analysis_output/Deseq2/c_vs_g.svg", width = 15, height = 16)
ddsout_c_g$graph
dev.off()


svg("./field_soil_analysis_output/Deseq2/o_vs_g.svg", width = 15, height = 16)
ddsout_o_g$graph
dev.off()

# Fungi
svg("./field_soil_analysis_output/Deseq2/c_vs_o_its.svg", width = 15, height = 16)
ddsout_c_o_its$graph
dev.off()


svg("./field_soil_analysis_output/Deseq2/c_vs_g_its.svg", width = 15, height = 16)
ddsout_c_g_its$graph
dev.off()


svg("./field_soil_analysis_output/Deseq2/o_vs_g_its.svg", width = 15, height = 16)
ddsout_o_g_its$graph
dev.off()

library(ggpubr)

svg("./field_soil_analysis_output/Deseq2/combined_16s.svg", width=15, height = 17)
ggarrange(ddsout_c_o$graph, ddsout_c_g$graph, ddsout_o_g$graph, common.legend = T, align = "hv", nrow = 3, labels = c("A", "", ""))
dev.off()


svg("./field_soil_analysis_output/Deseq2/combined_its.svg", width=15, height = 17)
ggarrange(ddsout_c_o_its$graph, ddsout_c_g_its$graph, ddsout_o_g_its$graph, common.legend = T, align = "hv", nrow = 3, labels = c("B", "", ""))
dev.off()

################################################
# field data plot tree
rm(list=ls())
source("./tools/PlotTree.R")
library(ggtree)
fdata_16s.prune = readRDS("./field_soil_analysis_output/data/fdata_16s_pruned.rds")

################################################

# loading signficant ASVs

all_sig.16s = c()
for(i in list.files("./field_soil_analysis_output/Deseq2/sigout/16s/")){
        if(endsWith(i, "level0.01.csv")){
                all_sig.16s = c(all_sig.16s, i)
        } 
}

all_sig.16s

cultype_df.16s = lapply(all_sig.16s, addTime, filepath = "./field_soil_analysis_output/Deseq2/sigout/16s/", patterns = "^.+?_vs_.+?")

cultype_df2.16s = do.call(rbind, cultype_df.16s)

nrow(cultype_df2.16s) # 1565

# select fold change >= 1
cultype_df3.16s = cultype_df2.16s[which(abs(cultype_df2.16s$log2FoldChange) >= 1), ]

ASVs_sig.16s = as.character(cultype_df3.16s$X)
head(ASVs_sig.16s)

allsig.16s = prune_taxa(ASVs_sig.16s, fdata_16s.prune)



# plot tree
tax_table(allsig.16s)[, "Genus"] = gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "A-N-P-Rhizobium", tax_table(allsig.16s)[, "Genus"])

UseColors = colorRampPalette(brewer.pal(9, "Set1"))

Phylum16s = c(tax_table(allsig.16s)[, "Phylum"])
assignColors = data.frame(Phylum = unique(Phylum16s), Colors = UseColors(length(unique(Phylum16s))))
assignCol = as.character(assignColors[, "Colors"])
assignCol
names(assignCol) = unique(Phylum16s)

head(cultype_df2.16s)
# change timestamp
cultype_df3.16s$timestamp = gsub("Conventional", "C", cultype_df3.16s$timestamp)
cultype_df3.16s$timestamp = gsub("Organic", "O", cultype_df3.16s$timestamp)
cultype_df3.16s$timestamp = gsub("Green", "G", cultype_df3.16s$timestamp)

# create matrix
clean_meta_16s = cultype_df3.16s %>%
        select(X, timestamp, log2FoldChange) %>%
        spread(timestamp, log2FoldChange)

# NA sets to 0
clean_meta_16s[is.na(clean_meta_16s)] = 0
# make matrix  
clean_meta_mtrx_16s = as.matrix(clean_meta_16s[, c("C_vs_G", "C_vs_O", "O_vs_G")])
# chane row names
rownames(clean_meta_mtrx_16s) = clean_meta_16s$X

n = length(unique(tax_table(allsig.16s)[, "Phylum"]))
head(clean_meta_mtrx_16s)
clean_meta_mtrx_16s[1:2, 1:3]

# plot the tree
gtree16s = ggtree(allsig.16s, ladderize = F, branch.length = "none")+ geom_tiplab(aes(label="")) + geom_tippoint(aes(color=Phylum), size = 3)

colnames(clean_meta_mtrx_16s) = c("C_vs_M", "C_vs_O", "O_vs_M")
clean_meta_mtrx_16s[1:2, 1:3]

# plot heatmap
# library(ggnewscale)

# svg("./field_soil_analysis_output/Deseq2/gheamap_16s.svg", width = 12, height = 15)

svg("./field_soil_analysis_output/Deseq2/gheamap_16s_colorscaled.svg", width = 12, height = 15)
(tree1 = gheatmap(gtree16s, clean_meta_mtrx_16s, legend_title = "Log2 Fold Change", colnames = T, offset= -4.5, colnames_offset_y = -8) + 
                scale_color_manual(values = assignCol) +
                scale_fill_gradient2(low = "steelblue", mid = "grey90", high = "red") +
                guides(color = guide_legend(order=0), fill = guide_legend(order = 1)))

dev.off()

# summary

# C_vs_G
# cultype_df3.16s %>%
#         filter(timestamp == "C_vs_G" & log2FoldChange > 0) %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq) %>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/enriched_comp_16s_C_vs_G.csv", row.names = F, quote = F)
# 
# cultype_df3.16s %>%
#         filter(timestamp == "C_vs_G" & log2FoldChange > 0) %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq)%>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/enriched_comp_16s_C_vs_G_genus.csv", row.names = F, quote = F)

cultype_df3.16s %>%
        filter(timestamp == "C_vs_G" & log2FoldChange > 0) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)

cultype_df3.16s %>%
        filter(timestamp == "C_vs_G" & log2FoldChange < 0) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)

#----

cultype_df3.16s %>%
        filter(timestamp == "C_vs_G" & log2FoldChange > 0) %>%
        group_by(Family) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)

cultype_df3.16s %>%
        filter(timestamp == "C_vs_G" & log2FoldChange < 0) %>%
        group_by(Family) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)

# cultype_df3.16s %>%
#         filter(timestamp == "C_vs_G" & log2FoldChange < 0) %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq) %>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/depleted_comp_16s_C_vs_G.csv", row.names = F, quote = F)
# 
# cultype_df3.16s %>%
#         filter(timestamp == "C_vs_G" & log2FoldChange < 0) %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq)%>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/depleted_comp_16s_C_vs_G_genus.csv", row.names = F, quote = F)


# C_vs_O
# cultype_df3.16s %>%
#         filter(timestamp == "C_vs_O" & log2FoldChange > 0) %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq) %>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/enriched_comp_16s_C_vs_O.csv", row.names = F, quote = F)
# 
# cultype_df3.16s %>%
#         filter(timestamp == "C_vs_O" & log2FoldChange > 0) %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq) %>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/enriched_comp_16s_C_vs_O_genus.csv", row.names = F, quote = F)

cultype_df3.16s %>%
        filter(timestamp == "C_vs_O" & log2FoldChange >= 2) %>%
        group_by(Phylum) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)


cultype_df3.16s %>%
        filter(timestamp == "C_vs_O" & log2FoldChange >= 2) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)

cultype_df3.16s %>%
        filter(timestamp == "C_vs_O" & log2FoldChange < 0) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)


cultype_df3.16s %>%
        filter(timestamp == "C_vs_O" & log2FoldChange < 0) %>%
        group_by(Family) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)

# cultype_df3.16s %>%
#         filter(timestamp == "C_vs_O" & log2FoldChange < 0) %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq) %>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/depleted_comp_16s_C_vs_O.csv", row.names = F, quote = F)
# 
# cultype_df3.16s %>%
#         filter(timestamp == "C_vs_O" & log2FoldChange < 0) %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq) %>% write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/depleted_comp_16s_C_vs_O_genus.csv", row.names = F, quote = F)


# O_vs_G
# cultype_df3.16s %>%
#         filter(timestamp == "O_vs_G" & log2FoldChange > 0) %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq)%>% write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/enriched_comp_16s_O_vs_G.csv", row.names = F, quote = F)
# 
# cultype_df3.16s %>%
#         filter(timestamp == "O_vs_G" & log2FoldChange > 0) %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq) %>% write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/enriched_comp_16s_O_vs_G_genus.csv", row.names = F, quote = F)

cultype_df3.16s %>%
        filter(timestamp == "O_vs_G" & log2FoldChange > 0) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)

cultype_df3.16s %>%
        filter(timestamp == "O_vs_G" & log2FoldChange < 0) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)



# cultype_df3.16s %>%
#         filter(timestamp == "O_vs_G" & log2FoldChange < 0) %>%
#         group_by(Phylum) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq) %>% write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/depleted_comp_16s_O_vs_G.csv", row.names = F, quote = F)
# 
# cultype_df3.16s %>%
#         filter(timestamp == "O_vs_G" & log2FoldChange < 0) %>%
#         group_by(Genus) %>%
#         tally() %>%
#         mutate(freq = 100* n/sum(n)) %>%
#         arrange(-freq)%>%
#         arrange(-freq) %>% write.csv("./field_soil_analysis_output/Deseq2/sigout/16s/summary_sigASV_in_16s/depleted_comp_16s_O_vs_G_genus.csv", row.names = F, quote = F)

#--------------------------------------------
# family level

# fam16s = c(tax_table(allsig.16s)[, "Family"])
# assignColors = data.frame(Family = unique(fam16s), Colors = UseColors(length(unique(fam16s))))
# assignCol = as.character(assignColors[, "Colors"])
# assignCol
# names(assignCol) = unique(fam16s)
# 
# gtree16s_fam = ggtree(allsig.16s, ladderize = F, branch.length = "none")+ geom_tiplab(aes(label="")) + geom_tippoint(aes(color=Family), size = 3)
# svg("./field_soil_analysis_output/Deseq2/gheamap_16s_fam.svg", width = 12, height = 15)
# (tree2_fam = gheatmap(gtree16s_fam, clean_meta_mtrx_16s, legend_title = "Log2 Fold Change", colnames = T, offset= -5, colnames_offset_y = -8) + guides(color = guide_legend(order=0), fill = guide_legend(order = 1)) + theme(legend.position = "bottom"))
# dev.off()


# scale_color_manual(values = glasbey())+ theme(legend.position = "right")

# save graphs

# svg("./field_soil_analysis_output/Deseq2/DESeq2_tree2_16S.svg", width = 20, height = 25)
# (C16s_graph = plotggtree(allsig.16s, cultype_df3.16s, otucol = "X", colFUN = assignCol, colselected = c("C_vs_G", "C_vs_O", "O_vs_G"), addColnames = T, hmap_offset = 1, hmap_width = 0.8, colnames_offset_y = -0.8, xlimTree = NULL)+ guides(color = guide_legend(order=0), fill = guide_legend(order = 1)))
# dev.off()
# 



#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
#Fungi tree
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

# analyzing tree ITS

# loading signficant ASVs
rm(list=ls())
fdata_its.prune = readRDS("./field_soil_analysis_output/data/fdata_its_pruned.rds")
source("./tools/PlotTree.R")

all_sig.its = c()
for(i in list.files("./field_soil_analysis_output/Deseq2/sigout/its/")){
        if(endsWith(i, "level0.01.csv")){
                all_sig.its = c(all_sig.its, i)
        } 
}

all_sig.its

cultype_df.its = lapply(all_sig.its, addTime, filepath = "./field_soil_analysis_output/Deseq2/sigout/its/", patterns = "^.+?_vs_.+?")

cultype_df2.its = do.call(rbind, cultype_df.its)
nrow(cultype_df2.its) #314
# select fold change >= 1
cultype_df3.its = cultype_df2.its[which(abs(cultype_df2.its$log2FoldChange) >= 1), ]

ASVs_sig.its = as.character(cultype_df3.its$X)
head(ASVs_sig.its)

allsig.its = prune_taxa(ASVs_sig.its, fdata_its.prune)


UseColors = colorRampPalette(brewer.pal(9, "Set1"))

PhylumITS = c(tax_table(allsig.its)[, "Phylum"])
assignColors = data.frame(Phylum = unique(PhylumITS), Colors = UseColors(length(unique(PhylumITS))))
assignCol = as.character(assignColors[, "Colors"])
assignCol
names(assignCol) = unique(PhylumITS)

head(cultype_df3.its)
# change timestamp
cultype_df3.its$timestamp = gsub("Conventional", "C", cultype_df3.its$timestamp)
cultype_df3.its$timestamp = gsub("Organic", "O", cultype_df3.its$timestamp)
cultype_df3.its$timestamp = gsub("Green", "G", cultype_df3.its$timestamp)

# create matrix
clean_meta_its = cultype_df3.its %>%
        select(X, timestamp, log2FoldChange) %>%
        spread(timestamp, log2FoldChange)

# NA sets to 0
clean_meta_its[is.na(clean_meta_its)] = 0
# make matrix  
clean_meta_mtrx_its = as.matrix(clean_meta_its[, c("C_vs_G", "C_vs_O", "O_vs_G")])
# chane row names
rownames(clean_meta_mtrx_its) = clean_meta_its$X

n = length(unique(tax_table(allsig.its)[, "Phylum"]))
head(clean_meta_mtrx_its)
colnames(clean_meta_mtrx_its) = c("C_vs_M", "C_vs_O", "O_vs_M")

# plot tree
gtreeits = ggtree(allsig.its, ladderize = F, branch.length = "none")+ geom_tiplab(aes(label="")) + geom_tippoint(aes(color=Phylum), size = 3)

# plot heatmap
svg("./field_soil_analysis_output/Deseq2/gheamap_its_colorscaled.svg", width = 12, height = 15)
(tree2 = gheatmap(gtreeits, clean_meta_mtrx_its, legend_title = "Log2 Fold Change", colnames = T, offset= -4.2, colnames_offset_y = -2) + 
                scale_color_manual(values = assignCol) +
                scale_fill_gradient2(low="steelblue", mid = "grey90", high = "red") +
                guides(color = guide_legend(order=0), fill = guide_legend(order = 1)))
dev.off()


# save graphs

# svg("./field_soil_analysis_output/Deseq2/DESeq2_tree2_ITS.svg", width = 20, height = 23)
# (plotggtree(allsig.its, cultype_df2.its, otucol = "X", colFUN = assignCol, colselected = c("C_vs_G", "C_vs_O", "O_vs_G"), addColnames = T, hmap_offset = 1, hmap_width = 0.8, colnames_offset_y = -0.8, xlimTree = NULL)+ guides(color = guide_legend(order=0), fill = guide_legend(order = 1)))
# dev.off()
# summary

# C_vs_G
cultype_df3.its %>%
        filter(timestamp == "C_vs_G" & log2FoldChange >= 2) %>%
        group_by(Phylum) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq) 
# %>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/enriched_comp_its_C_vs_G.csv", row.names = F, quote = F)

cultype_df3.its %>%
        filter(timestamp == "C_vs_G" & log2FoldChange >= 2) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)
# %>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/enriched_comp_its_C_vs_G_genus.csv", row.names = F, quote = F)


cultype_df3.its %>%
        filter(timestamp == "C_vs_G" & log2FoldChange <= -2) %>%
        group_by(Phylum) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq) 
# %>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/depleted_comp_its_C_vs_G.csv", row.names = F, quote = F)

cultype_df3.its %>%
        filter(timestamp == "C_vs_G" & log2FoldChange <= -2) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)
# %>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/depleted_comp_its_C_vs_G_genus.csv", row.names = F, quote = F)


# C_vs_O
cultype_df3.its %>%
        filter(timestamp == "C_vs_O" & log2FoldChange >=2) %>%
        group_by(Phylum) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq) 
# %>%
        # write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/enriched_comp_its_C_vs_O.csv", row.names = F, quote = F)

cultype_df3.its %>%
        filter(timestamp == "C_vs_O" & log2FoldChange >= 2) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq) 
# %>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/enriched_comp_its_C_vs_O_genus.csv", row.names = F, quote = F)


cultype_df3.its %>%
        filter(timestamp == "C_vs_O" & log2FoldChange <= -2) %>%
        group_by(Phylum) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq) 
# %>%
#         write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/depleted_comp_its_C_vs_O.csv", row.names = F, quote = F)

cultype_df3.its %>%
        filter(timestamp == "C_vs_O" & log2FoldChange <= -2) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) 
# %>%
#         arrange(-freq) %>% write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/depleted_comp_its_C_vs_O_genus.csv", row.names = F, quote = F)


# O_vs_G
cultype_df3.its %>%
        filter(timestamp == "O_vs_G" & log2FoldChange >= 2) %>%
        group_by(Phylum) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)
# %>% write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/enriched_comp_its_O_vs_G.csv", row.names = F, quote = F)

cultype_df3.its %>%
        filter(timestamp == "O_vs_G" & log2FoldChange >= 2) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq) 
# %>% write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/enriched_comp_its_O_vs_G_genus.csv", row.names = F, quote = F)


cultype_df3.its %>%
        filter(timestamp == "O_vs_G" & log2FoldChange <= -2) %>%
        group_by(Phylum) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq) 
# %>% write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/depleted_comp_its_O_vs_G.csv", row.names = F, quote = F)

cultype_df3.its %>%
        filter(timestamp == "O_vs_G" & log2FoldChange <= -2) %>%
        group_by(Genus) %>%
        tally() %>%
        mutate(freq = 100* n/sum(n)) %>%
        arrange(-freq)

# %>% write.csv("./field_soil_analysis_output/Deseq2/sigout/its/summary_sigASV_in_ITS/depleted_comp_its_O_vs_G_genus.csv", row.names = F, quote = F)





# summarize
source("./tools/plotDeseq2.R")

#BBBBBBBBBBBBBB
# Bacterial
#BBBBBBBBBBBBB

# conv vs organic

summ_deseq2(ddsout_c_o)
summ_deseq2(ddsout_c_o, "Genus")


# conv vs green
summ_deseq2(ddsout_c_g)
summ_deseq2(ddsout_c_g, "Genus")

# organic vs green
summ_deseq2(ddsout_o_g)
summ_deseq2(ddsout_o_g, "Genus")


#FFFFFFFFFFFFFFFFFFFFFFF
# Fungi
#FFFFFFFFFFFFFFFFFFFFFF
# conv vs organic

summ_deseq2(ddsout_c_o_its)
summ_deseq2(ddsout_c_o_its, "Genus")


# conv vs green
summ_deseq2(ddsout_c_g_its)
summ_deseq2(ddsout_c_g_its, "Genus")

# organic vs green
summ_deseq2(ddsout_o_g_its)
summ_deseq2(ddsout_o_g_its, "Genus")


##########################################################################################

##########################################################################################

##########################################################################################



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# microcosm analysis

rm(list = ls())

dir.create("./microcosm_soil_analysis_output")
# import data

# # # import biom files from qiime2
mydata_16s = import_biom("../source_data/mcm/16s/feature-table-tax.biom ")
tree.16s = read_tree("../source_data/mcm/16s/tree.nwk")
phy_tree(mydata_16s) = tree.16s
meta_16s = data.frame(sample_data(mydata_16s))

mydata_its = import_biom("../source_data/mcm/its/feature-table-tax.biom")
tree.its = read_tree("../source_data/mcm/its/tree.nwk")
phy_tree(mydata_its) =  tree.its
meta_its = data.frame(sample_data(mydata_its))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# data processing

#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# Bacteria
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB


# # taxonomy header converted to bASV+number
tax_table(mydata_16s) = gsub("D_[0-9]__", "", tax_table(mydata_16s)) # remove prefix
colnames(tax_table(mydata_16s)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa_names(mydata_16s) = paste0("bASV", seq(1, nrow(otu_table(mydata_16s))))


pathogen_cols = c('Meloidogyne_hapla',
                  'Meloidogyne_chitwoodi',
                  'Pratylenchus_penetrans',
                  'Pratylenchus_neglectus',
                  'Pratylenchus_thornei',
                  'Paratrichodorus_sp',
                  'Tylenchorhynchus_sp',
                  'Paratylenchus_sp',
                  'Helicotylenchus_sp',
                  'Mesocriconema_sp',
                  'Ditylenchus_sp',
                  'Xiphinema_sp',
                  'Hemicycliophora_sp',
                  'vert_count',
                  'bd_count')

numeric_cols = c(pathogen_cols, "Period_numeric", "Soil_pH_Value")
 
colnames(meta_16s)[which(colnames(meta_16s) == "MesocriconemaÂ.sp")] ="Mesocriconema_sp"

# # reorder
#
cols.16s = colnames(meta_16s)
newcols.16s = c("SampleID", cols.16s)
meta_16s$SampleID = rownames(meta_16s)
meta_16s = meta_16s[newcols.16s]


# change pathogen count to numeric
for(i in numeric_cols){
        meta_16s[i] = lapply(meta_16s[i], as.numeric)
}

factor_cols = c("Period", "Soil_depth", "cul_type", "soil_compname","soil_taxclass","soil_taxgrtgroup","soil_taxorder","soil_taxsubgrp","soil_taxsuborder")

for(i in factor_cols){
        meta_16s[i] = lapply(meta_16s[i], as.factor)
}


# change field and microcosm naming
meta_16s$Field = gsub("Boardman", "B", meta_16s$Field)
meta_16s$Field = gsub("Organic", "O", meta_16s$Field)

meta_16s$microcosms = gsub("Boardman", "B", meta_16s$microcosms)
meta_16s$microcosms = gsub("Organic", "O", meta_16s$microcosms)


# all lower case column names
colnames(meta_16s) = tolower(colnames(meta_16s))

# calculate rotation years

yrs = paste0("y", seq(1998, 2018))

totalCrops = apply(meta_16s[, yrs], 1, function(x){

        x = as.vector(unlist(x))
        x[x==""] = NA
        # some rotation has - in them
        allcrops = unlist(str_split(x, "-"))
        uniqueCrops = unique(na.omit(allcrops))
        if(uniqueCrops[1] == "unknown"){
                return(-999)
        } else {
                return(length(uniqueCrops))
        }
})



totalRotYrs = apply(meta_16s[, yrs], 1, function(x){

        # return not ''
        x[x==""] = NA
        allyears = na.omit(unlist(x))
        if(unique(allyears)[1] == "unknown"){
                return(-999)
        }else{
                return(length(allyears))
        }

})


totalRotYrs_Pota = apply(meta_16s[, yrs], 1, function(x){

        x[x==""] = NA
        # some rotation has - in them
        allcrops = na.omit(unlist(str_split(x, "-")))
        if(unique(allcrops)[1] == "unknown"){
                return(-999)
        } else {
                return(sum(tolower(allcrops)=="potato"))
        }

})


percPota = totalRotYrs_Pota / totalRotYrs

meta_16s = cbind(meta_16s, totalCrops, totalRotYrs, totalRotYrs_Pota, percPota)

# new meta data
sample_data(mydata_16s) = meta_16s

dir.create("./microcosm_soil_analysis_output/data")
saveRDS(mydata_16s, "./microcosm_soil_analysis_output/data/mydata_16s.rds")

#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# fungi
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

tax_table(mydata_its) = gsub("k__|p__|c__|o__|f__|g__|s__", "", tax_table(mydata_its))
colnames(tax_table(mydata_its)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

taxa_names(mydata_its) = paste0("fASV", seq(1, nrow(otu_table(mydata_its))))


# change to numbers or factors
pathogen_cols = c('Meloidogyne_hapla',
                  'Meloidogyne_chitwoodi',
                  'Pratylenchus_penetrans',
                  'Pratylenchus_neglectus',
                  'Pratylenchus_thornei',
                  'Paratrichodorus_sp',
                  'Tylenchorhynchus_sp',
                  'Paratylenchus_sp',
                  'Helicotylenchus_sp',
                  'Mesocriconema_sp',
                  'Ditylenchus_sp',
                  'Xiphinema_sp',
                  'Hemicycliophora_sp',
                  'vert_count',
                  'bd_count')
#
numeric_cols = c(pathogen_cols, "Period_numeric", "Soil_pH_Value")


colnames(meta_its)[which(colnames(meta_its) == "MesocriconemaÂ.sp")] = "Mesocriconema_sp"

# reorder put SampleID first
cols.its = colnames(meta_its)
newcols.its = c("SampleID", cols.its)
meta_its$SampleID = rownames(meta_its)
meta_its = meta_its[newcols.its]
#
#
# change those pathogen counts to numeric
for(i in numeric_cols){
        meta_its[i] = lapply(meta_its[i], as.numeric)
}

factor_cols = c("Period", "Soil_depth", "cul_type", "soil_compname","soil_taxclass","soil_taxgrtgroup","soil_taxorder","soil_taxsubgrp","soil_taxsuborder")

# change those to factors
for(i in factor_cols){
        meta_its[i] = lapply(meta_its[i], as.factor)
}


meta_its$Field = gsub("Boardman", "B", meta_its$Field)
meta_its$Field = gsub("Organic", "O", meta_its$Field)

meta_its$microcosms = gsub("Boardman", "B", meta_its$microcosms)
meta_its$microcosms = gsub("Organic", "O", meta_its$microcosms)

# all lower case column names, except for the first
colnames(meta_its) = tolower(colnames(meta_its))

totalCrops = apply(meta_its[, yrs], 1, function(x){

        x = as.vector(unlist(x))
        x[x==""] = NA
        # some rotation has - in them
        allcrops = unlist(str_split(x, "-"))
        uniqueCrops = unique(na.omit(allcrops))
        if(uniqueCrops[1] == "unknown"){
                return(-999)
        } else {
                return(length(uniqueCrops))
        }
})


totalRotYrs = apply(meta_its[, yrs], 1, function(x){

        # return not ''
        x[x==""] = NA
        allyears = na.omit(unlist(x))
        if(unique(allyears)[1] == "unknown"){
                return(-999)
        }else{
                return(length(allyears))
        }

})


totalRotYrs_Pota = apply(meta_its[, yrs], 1, function(x){

        x[x==""] = NA
        # some rotation has - in them
        allcrops = na.omit(unlist(str_split(x, "-")))
        if(unique(allcrops)[1] == "unknown"){
                return(-999)
        } else {
                return(sum(tolower(allcrops)=="potato"))
        }

})

percPota = totalRotYrs_Pota / totalRotYrs
meta_its = cbind(meta_its, totalCrops, totalRotYrs, totalRotYrs_Pota, percPota)

sample_data(mydata_its) = meta_its

# save to orignal files folder
saveRDS(mydata_its, "./microcosm_soil_analysis_output/data/mydata_its.rds")



#################################################################################
# summarize reads
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")
#################################################################################

mydata_16s #124696 taxa 574 samples
mydata_its # 16417 taxa 576 samples


sample_data(mydata_16s)


# total reads
sum(sample_sums(mydata_16s)) # 15,166,470
sum(sample_sums(mydata_its)) # 25,541,570


table(tax_table(mydata_16s)[, "Phylum"], exclude= F) # 1670 NA, 
table(tax_table(mydata_its)[, "Phylum"], exclude = F) # 5012 NA, 1380 unidentified


# make rarefaction plot

# check how many sample would be removed check rarefaction curve
myamp.16s = phyloseq_to_ampvis2(mydata_16s)
myamp.its = phyloseq_to_ampvis2(mydata_its)


curv16s = amp_rarecurve(myamp.16s, color_by = "microcosms", stepsize = 100) 
g1 = curv16s + theme(legend.position = "none", plot.title = element_text(hjust=0.5)) + labs(title="16S rarefaction curve", y="Number of observed ASVs")

curvits = amp_rarecurve(myamp.its, color_by = "microcosms", stepsize = 100) 
g2 = curvits + theme(legend.position = "none", plot.title = element_text(hjust=0.5)) + labs(title="ITS2 rarefaction curve", y="Number of observed ASVs")


g1g2 = ggpubr::ggarrange(plotlist = list(g1, g2), labels = c("A", "B"), nrow = 1)

g1g2

dir.create("./microcosm_soil_analysis_output/rarecurves")

svg("./microcosm_soil_analysis_output/rarecurves/rarecurves_mcm.svg", width = 12, height = 10)
g1g2
dev.off()




##################################################################################################

# prune data

# remove NA, and rare species for downstream (some algorithm maybe sensitive to rare species)
# set at 10

mydata_16s.prune = subset_taxa(mydata_16s, !is.na(Phylum))
mydata_16s.prune = prune_taxa(taxa_sums(mydata_16s.prune)>10, mydata_16s.prune)

mydata_its.prune = subset_taxa(mydata_its, !is.na(Phylum))
mydata_its.prune = prune_taxa(taxa_sums(mydata_its.prune)>10, mydata_its.prune)


# for 16s we will remove samples less than 7000 reads
# for ITS we will remove smaples less than 10000 reads

mydata_16s.prune = prune_samples(sample_sums(mydata_16s.prune) >= 7000, mydata_16s.prune)
min(sample_sums(mydata_16s.prune)) # 7019

mydata_its.prune = prune_samples(sample_sums(mydata_its.prune) >= 10000, mydata_its.prune)
min(sample_sums(mydata_its.prune)) # 10476



# use microbiomeutility to find the best taxonomy rank for NA
mbst16s = format_to_besthit(mydata_16s.prune)
mbstITS = format_to_besthit(mydata_its.prune)

tax_table(mbst16s)[, 1:7]

mbst16s_taxtb = tax_table(mbst16s)[, 1:7]
rownames(mbst16s_taxtb) = gsub("OTU-|:.*$", "", rownames(mbst16s_taxtb))

mbstits_taxtb = tax_table(mbstITS)[, 1:7]
rownames(mbstits_taxtb) = gsub("OTU-|:.*$", "", rownames(mbstits_taxtb))


tax_table(mydata_16s.prune) = mbst16s_taxtb
tax_table(mydata_its.prune) = mbstits_taxtb


saveRDS(mydata_16s.prune, file.path("./microcosm_soil_analysis_output/data/", "mydata_16s_pruned.rds"))
saveRDS(mydata_its.prune, file.path("./microcosm_soil_analysis_output/data/", "mydata_its_pruned.rds"))


#####################################################################################################

## alpha diversity analysis

rm(list=ls())
dir.create("./microcosm_soil_analysis_output/alpha_diversity")

mydata_16s.prune = readRDS("./microcosm_soil_analysis_output/data/mydata_16s_pruned.rds")
mydata_its.prune = readRDS("./microcosm_soil_analysis_output/data/mydata_its_pruned.rds")


mydata_16s.prune # 574 - 555
mydata_its.prune # 576 - 567
574-555
576-567

meta16s = as(sample_data(mydata_16s.prune), "data.frame")
metaITS = as(sample_data(mydata_its.prune), "data.frame")

trial = 100

#####################################################################################################


# estimate indices




# 16s

observed_mat_16s = matrix(nrow = nrow(sample_data(mydata_16s.prune)), ncol = trial)
row.names(observed_mat_16s) = sample_names(mydata_16s.prune)

shannon_mat_16s = matrix(nrow = nrow(sample_data(mydata_16s.prune)), ncol = trial)
row.names(shannon_mat_16s) = sample_names(mydata_16s.prune)

invSimpson_mat_16s = matrix(nrow = nrow(sample_data(mydata_16s.prune)), ncol = trial)
row.names(invSimpson_mat_16s) = sample_names(mydata_16s.prune)


min_lib_16s = min(sample_sums(mydata_16s.prune))

set.seed(101)
for(i in 1:trial){

        print(i)
        ra_16s = rarefy_even_depth(mydata_16s.prune, sample.size = min_lib_16s, replace=F)

        est_16s = as.matrix(estimate_richness(ra_16s, measures = c("Observed", "Shannon", "InvSimpson")))
        observed_mat_16s[, i] = as.numeric(est_16s[,1])
        shannon_mat_16s[,i] = as.numeric(est_16s[,2])
        invSimpson_mat_16s[,i] = as.numeric(est_16s[,3])

}

SampleID = row.names(observed_mat_16s)
mean = apply(observed_mat_16s, 1, mean)
sd = apply(observed_mat_16s, 1, sd)
measure = rep("Observed", nsamples(mydata_16s.prune))
observed_stats_16s = data.frame(SampleID, mean, sd, measure)

SampleID = row.names(shannon_mat_16s)
mean = apply(shannon_mat_16s, 1, mean)
sd = apply(shannon_mat_16s, 1, sd)
measure = rep("Shannon", nsamples(mydata_16s.prune))
shannon_stats_16s = data.frame(SampleID, mean, sd, measure)


SampleID = row.names(invSimpson_mat_16s)
mean = apply(invSimpson_mat_16s, 1, mean)
sd = apply(invSimpson_mat_16s, 1, sd)
measure = rep("invSimpson", nsamples(mydata_16s.prune))
invSimpson_stats_16s = data.frame(SampleID, mean, sd, measure)

alpha_est_16s = rbind(observed_stats_16s, shannon_stats_16s, invSimpson_stats_16s)

a_diversity_16s = merge(meta16s, alpha_est_16s, by.x = "sampleid", by.y = "SampleID")


saveRDS(a_diversity_16s, "./microcosm_soil_analysis_output/alpha_diversity/alpha_diversity_16s.rds")


#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# Fungi
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF


observed_mat_its = matrix(nrow = nrow(sample_data(mydata_its.prune)), ncol = trial)
row.names(observed_mat_its) = sample_names(mydata_its.prune)
shannon_mat_its = matrix(nrow = nrow(sample_data(mydata_its.prune)), ncol = trial)
row.names(shannon_mat_its) = sample_names(mydata_its.prune)
invSimpson_mat_its = matrix(nrow = nrow(sample_data(mydata_its.prune)), ncol = trial)
row.names(invSimpson_mat_its) = sample_names(mydata_its.prune)

min_lib_its = min(sample_sums(mydata_its.prune))

set.seed(100)
for(i in 1:trial){

        print(i)
        ra_its = rarefy_even_depth(mydata_its.prune, sample.size = min_lib_its, replace=F)

        est_its = as.matrix(estimate_richness(ra_its, measures = c("Observed", "Shannon", "InvSimpson")))
        observed_mat_its[, i] = as.numeric(est_its[,1])
        shannon_mat_its[,i] = as.numeric(est_its[,2])
        invSimpson_mat_its[,i] = as.numeric(est_its[,3])

}

SampleID = row.names(observed_mat_its)
mean = apply(observed_mat_its, 1, mean)
sd = apply(observed_mat_its, 1, sd)
measure = rep("Observed", nsamples(mydata_its.prune))
observed_stats_its = data.frame(SampleID, mean, sd, measure)

SampleID = row.names(shannon_mat_its)
mean = apply(shannon_mat_its, 1, mean)
sd = apply(shannon_mat_its, 1, sd)
measure = rep("Shannon", nsamples(mydata_its.prune))
shannon_stats_its = data.frame(SampleID, mean, sd, measure)


SampleID = row.names(invSimpson_mat_its)
mean = apply(invSimpson_mat_its, 1, mean)
sd = apply(invSimpson_mat_its, 1, sd)
measure = rep("invSimpson", nsamples(mydata_its.prune))
invSimpson_stats_its = data.frame(SampleID, mean, sd, measure)

alpha_est_its = rbind(observed_stats_its, shannon_stats_its, invSimpson_stats_its)
a_diversity_its = merge(data.frame(sample_data(mydata_its.prune)), alpha_est_its, by.x = "sampleid", by.y = "SampleID")

saveRDS(a_diversity_its, "./microcosm_soil_analysis_output/alpha_diversity/alpha_diversity_its.rds")


###################################################

# statistical analysis of alpha diversity

rm(list=ls())

alpha16s = readRDS("./microcosm_soil_analysis_output/alpha_diversity/alpha_diversity_16s.rds")
alphaITS = readRDS("./microcosm_soil_analysis_output/alpha_diversity/alpha_diversity_its.rds")
###################################################

# change shannon to exponential

alpha16s$mean[alpha16s$measure == "Shannon"] = exp(alpha16s$mean[alpha16s$measure == "Shannon"])
levels(alpha16s$measure) = c("Observed richness", expression(e^Shannon), "inverse Simpson")


alphaITS$mean[alphaITS$measure == "Shannon"] = exp(alphaITS$mean[alphaITS$measure == "Shannon"])
levels(alphaITS$measure) = c("Observed richness", expression(e^Shannon), "inverse Simpson")

str(alpha16s)

# change microcosm and plots to factors
alpha16s$fmicrocosms = factor(alpha16s$microcosms)
alpha16s$fplots = factor(alpha16s$plots)

alphaITS$fmicrocosms = factor(alphaITS$microcosms)
alphaITS$fplots = factor(alphaITS$plots)


# making plots




# bacterial

# culture practice
(ad.16s.l = ggplot(alpha16s, aes(x=period, y=mean, color=cul_type, group=cul_type))+
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                # facet_wrap( measure ~ soil_depth, scales="free_y") +
                facet_grid( measure ~ soil_depth, scales = "free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                scale_color_discrete(labels = c("Conventional", "Mixed", "Organic")) +
                labs(x = "Time (Week)", y = expression(paste(alpha, " diversity measures")), color="Cropping systems"))

# 2 soil types

(ad.16s.l.soiltype = ggplot(alpha16s, aes(x=period, y=mean, color=soil_compname, group=soil_compname))+
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                facet_wrap( measure ~ soil_depth, scales="free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                labs(x = "Time (Week)", y = expression(paste(alpha, " diversity measures")), color="Soil types"))



# soil pH
(ad.16s.l.pH = alpha16s %>%
        mutate(pH_range = ifelse(soil_ph_value > 7, "Alkaline(>7)", ifelse(soil_ph_value < 7, "Acidic(<7)", "Neutral(=7)"))) %>%
ggplot(aes(x=period, y=mean, color=pH_range, group=pH_range)) +
        stat_summary(fun=mean, geom="line", size=0.6)+
        stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
        stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
        scale_x_discrete(labels=c("0", "1", "3", "6")) +
        facet_wrap(measure ~ soil_depth, scales = "free_y") +
        theme_bw() +
        theme(text = element_text(size=17)) +
        labs(x = "Time (Week)", y = expression(paste(alpha, " diversity measures")), color="pH ranges"))


(ad.16s.l.percPota = alpha16s %>%
                mutate(potato_intensity = ifelse(percPota < .2 & percPota >= 0, "light (0~0.19)", ifelse(percPota >= .2 & percPota < .4, "medium (0.2~0.39)", ifelse(percPota >= .4 & percPota < .6, "strong (0.4~0.59)", "very strong (0.6~1)")))) %>%
                ggplot(aes(x=period, y=mean, color=potato_intensity, group=potato_intensity)) +
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                facet_wrap(measure ~ soil_depth, scales = "free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                labs(x = "Time (Week)", y = expression(paste(alpha, " diversity measures")), color="potat intensity"))


(ad.16s.l.percPota = alpha16s %>%
                mutate(potato_intensity = ifelse(percPota < .2 & percPota >= 0, "light (0~0.19)", ifelse(percPota >= .2 & percPota < .4, "medium (0.2~0.39)", ifelse(percPota >= .4 & percPota < .6, "strong (0.4~0.59)", "very strong (0.6~1)")))) %>%
                ggplot(aes(x=period, y=mean, color=potato_intensity, group=potato_intensity)) +
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                facet_wrap(measure ~ soil_depth, scales = "free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                labs(x = "Time (Week)", y = expression(paste(alpha, " diversity measures")), color="potat intensity"))


(ad.16s.l.totalcrops = alpha16s %>%
                mutate(crop_diversity = ifelse(totalCrops <= 3 & totalCrops >= 2, "light (2~3)", ifelse(totalCrops >= 4 & totalCrops <= 5, "medium (4~5)", "strong (6~7)"))) %>%
                ggplot(aes(x=period, y=mean, color=crop_diversity, group=crop_diversity)) +
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                facet_wrap(measure ~ soil_depth, scales = "free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                labs(x = "Time (Week)", y = expression(paste(alpha, " diversity measures")), color="Crop diversity"))


(ad.16s.l.rotyears = alpha16s %>%
                mutate(rot_years = ifelse(totalRotYrs <= 8 & totalRotYrs >= 6, "short (6~8)", ifelse(totalRotYrs >= 9 & totalRotYrs <= 12, "medium (9~12)", "long (13~21)"))) %>%
                ggplot(aes(x=period, y=mean, color=rot_years, group=rot_years)) +
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                facet_wrap(measure ~ soil_depth, scales = "free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                labs(x = "Time (Week)", y = expression(paste(alpha, " diversity measures")), color= "Rotation years"))




#########################################
# fungal
(ad.its.l = ggplot(alphaITS, aes(x=period, y=mean, color=cul_type, group=cul_type))+
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                # facet_wrap( measure ~ soil_depth, scales="free_y") +
                facet_grid( measure ~ soil_depth, scales = "free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                labs(x = "Time (Week)", y = NULL, color="Cropping systems") +
         scale_color_discrete(labels = c("Conventional", "Mixed", "Organic")))

(ad.its.l.soiltype = ggplot(alphaITS, aes(x=period, y=mean, color=soil_compname, group=soil_compname))+
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                facet_wrap( measure ~ soil_depth, scales="free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                labs(x = "Time (Week)", y = NULL, color="Soil types"))


# soil pH
(ad.its.l.pH = alphaITS %>%
                mutate(pH_range = ifelse(soil_ph_value > 7, "Alkaline(>7)", ifelse(soil_ph_value < 7, "Acidic(<7)", "Neutral(=7)"))) %>%
                ggplot(aes(x=period, y=mean, color=pH_range, group=pH_range)) +
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                facet_wrap(measure ~ soil_depth, scales = "free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                labs(x = "Time (Week)", y = expression(paste(alpha, " diversity measures")), color="pH ranges"))


(ad.its.l.percPota = alphaITS %>%
                mutate(potato_intensity = ifelse(percPota < .2 & percPota >= 0, "light (0~0.19)", ifelse(percPota >= .2 & percPota < .4, "medium (0.2~0.39)", ifelse(percPota >= .4 & percPota < .6, "strong (0.4~0.59)", "very strong (0.6~1)")))) %>%
                ggplot(aes(x=period, y=mean, color=potato_intensity, group=potato_intensity)) +
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                facet_wrap(measure ~ soil_depth, scales = "free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                labs(x = "Time (Week)", y = expression(paste(alpha, " diversity measures")), color="potat intensity"))


(ad.its.l.totalcrops = alphaITS %>%
                mutate(crop_diversity = ifelse(totalCrops <= 3 & totalCrops >= 2, "light (2~3)", ifelse(totalCrops >= 4 & totalCrops <= 5, "medium (4~5)", "strong (6~7)"))) %>%
                ggplot(aes(x=period, y=mean, color=crop_diversity, group=crop_diversity)) +
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                facet_wrap(measure ~ soil_depth, scales = "free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                labs(x = "Time (Week)", y = expression(paste(alpha, " diversity measures")), color="Crop diversity"))


(ad.its.l.rotyears = alphaITS %>%
                mutate(rot_years = ifelse(totalRotYrs <= 8 & totalRotYrs >= 6, "short (6~8)", ifelse(totalRotYrs >= 9 & totalRotYrs <= 12, "medium (9~12)", "long (13~21)"))) %>%
                ggplot(aes(x=period, y=mean, color=rot_years, group=rot_years)) +
                stat_summary(fun=mean, geom="line", size=0.6)+
                stat_summary(fun=mean, geom="point", size=2, alpha=0.5) +
                stat_summary(fun.data = mean_se, geom="errorbar",width=0.1) +
                scale_x_discrete(labels=c("0", "1", "3", "6")) +
                facet_wrap(measure ~ soil_depth, scales = "free_y") +
                theme_bw() +
                theme(text = element_text(size=17)) +
                labs(x = "Time (Week)", y = expression(paste(alpha, " diversity measures")), color= "Rotation years"))



lineplot_combined = ggpubr::ggarrange(ad.16s.l, ad.its.l, labels = c("A", "B"), common.legend = T, align = "h", legend = "bottom")

lineplot_combined_soiltypes = ggpubr::ggarrange(ad.16s.l.soiltype, 
                                                ad.its.l.soiltype, labels = c("A", "B"), common.legend = T, align = "h", legend = "bottom")

lineplot_combined_pH = ggpubr::ggarrange(ad.16s.l.pH, 
                                                ad.its.l.pH, labels = c("A", "B"), common.legend = T, align = "h", legend = "bottom")


lineplot_combined_percPota = ggpubr::ggarrange(ad.16s.l.percPota, 
                                         ad.its.l.percPota, labels = c("A", "B"), common.legend = T, align = "h", legend = "bottom")

lineplot_combined_crops = ggpubr::ggarrange(ad.16s.l.totalcrops, 
                                               ad.its.l.totalcrops, labels = c("A", "B"), common.legend = T, align = "h", legend = "bottom")


lineplot_combined_rotyears = ggpubr::ggarrange(ad.16s.l.rotyears, 
                                            ad.its.l.rotyears, labels = c("A", "B"), common.legend = T, align = "h", legend = "bottom")



#
# svg("./microcosm_soil_analysis_output/alpha_diversity/lineplot_combined.svg", width = 22, height = 17)
# lineplot_combined
# dev.off()

svg("./microcosm_soil_analysis_output/alpha_diversity/lineplot_combined_fixscale2.svg", width = 22, height = 17)
lineplot_combined
dev.off()


svg("./microcosm_soil_analysis_output/alpha_diversity/lineplot_combined_soiltypes.svg", width = 22, height = 17)
lineplot_combined_soiltypes
dev.off()


svg("./microcosm_soil_analysis_output/alpha_diversity/lineplot_combined_pH.svg", width = 22, height = 17)
lineplot_combined_pH
dev.off()



svg("./microcosm_soil_analysis_output/alpha_diversity/lineplot_combined_percPota.svg", width = 22, height = 17)
lineplot_combined_percPota
dev.off()

svg("./microcosm_soil_analysis_output/alpha_diversity/lineplot_combined_crops.svg", width = 22, height = 17)
lineplot_combined_crops
dev.off()

svg("./microcosm_soil_analysis_output/alpha_diversity/lineplot_combined_rotyears.svg", width = 22, height = 17)
lineplot_combined_rotyears
dev.off()

##################
# mixed models
##################

alpha16s$field = as.factor(alpha16s$field)
contrasts(alpha16s$period) = "contr.sum"
contrasts(alpha16s$soil_depth) = "contr.sum"
contrasts(alpha16s$cul_type) = "contr.sum"
contrasts(alpha16s$soil_compname) = "contr.sum"
contrasts(alpha16s$field) = "contr.sum"


alphaITS$field = as.factor(alphaITS$field)
contrasts(alphaITS$period) = "contr.sum"
contrasts(alphaITS$soil_depth) = "contr.sum"
contrasts(alphaITS$cul_type) = "contr.sum"
contrasts(alphaITS$soil_compname) = "contr.sum"
contrasts(alphaITS$field) = "contr.sum"

#BBBBBBBBBBBBBBBBBBBBBBBB
# Bacterial communites
#BBBBBBBBBBBBBBBBBBBBBBB

obs.16s =  subset(alpha16s, measure == "Observed richness")
shan.16s = subset(alpha16s, measure == "e^Shannon")
invS.16s = subset(alpha16s, measure == "inverse Simpson")

head(obs.16s)

##############
# cultype

# random effect changed to fields 2020/12/23
#############
lmm.obs.16s = lmer(mean ~ 1+soil_depth*period_numeric*cul_type + (1+period_numeric|field), data = obs.16s)
plot(lmm.obs.16s) # heterogeneity not bad
qqnorm(resid(lmm.obs.16s)) # normality not bad
qqline(resid(lmm.obs.16s))

hist(resid(lmm.obs.16s)) # looks normal

summ.lmm.obs.16s = summary(lmm.obs.16s)
summ.lmm.obs.16s
write.csv(summ.lmm.obs.16s$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/coef_lmm_obs_16s_fieldasrandom.csv")

Anova(lmm.obs.16s, type = 3, test.statistic = "F")

emmeans(lmm.obs.16s, pairwise ~ soil_depth)


############
# soil types
############

lmm.obs.16s.stype = lmer(mean ~ 1+soil_depth*period_numeric*soil_compname + (1+period_numeric|field), data = obs.16s)

plot(lmm.obs.16s.stype)
qqnorm(resid(lmm.obs.16s.stype)) # normality not bad
qqline(resid(lmm.obs.16s.stype))

hist(resid(lmm.obs.16s.stype)) # looks normal

(summ.lmm.obs.16s.stype = summary(lmm.obs.16s.stype))
summ.lmm.obs.16s.stype
write.csv(summ.lmm.obs.16s.stype$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_obs_16s_soiltype.csv")
a = Anova(lmm.obs.16s.stype, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_obs_16s_soiltypes.csv")
emmeans(lmm.obs.16s.stype, pairwise ~ soil_compname)


##########
# potato intensity
##########


lmm.obs.16s.pota = lmer(mean ~ 1+soil_depth*period_numeric*percPota + (1+period_numeric|field), data = obs.16s)

plot(lmm.obs.16s.pota)
qqnorm(resid(lmm.obs.16s.pota)) # normality not bad
qqline(resid(lmm.obs.16s.pota))

hist(resid(lmm.obs.16s.pota)) # looks normal

(summ.lmm.obs.16s.pota = summary(lmm.obs.16s.pota))
summ.lmm.obs.16s.pota
write.csv(summ.lmm.obs.16s.pota$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_obs_16s_potato.csv")
a = Anova(lmm.obs.16s.pota, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_obs_16s_potato.csv")

# emmeans(lmm.obs.16s.pota, pairwise ~ percPota)

######################
# crop diversity
#####################
lmm.obs.16s.crops = lmer(mean ~ 1+soil_depth*period_numeric*totalCrops + (1+period_numeric|field), data = obs.16s)

plot(lmm.obs.16s.crops)
qqnorm(resid(lmm.obs.16s.crops)) # normality not bad
qqline(resid(lmm.obs.16s.crops))

hist(resid(lmm.obs.16s.crops)) # looks normal

(summ.lmm.obs.16s.crops = summary(lmm.obs.16s.crops))
summ.lmm.obs.16s.crops
write.csv(summ.lmm.obs.16s.crops$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_obs_16s_crops.csv")
a = Anova(lmm.obs.16s.crops, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_obs_16s_crops.csv")

######################
# total years
#####################
lmm.obs.16s.rotyr = lmer(mean ~ 1+soil_depth*period_numeric*totalRotYrs + (1+period_numeric|fplots), data = obs.16s)

plot(lmm.obs.16s.rotyr)
qqnorm(resid(lmm.obs.16s.rotyr)) # normality not bad
qqline(resid(lmm.obs.16s.rotyr))

hist(resid(lmm.obs.16s.rotyr)) # looks normal

summ.lmm.obs.16s.rotyr = summary(lmm.obs.16s.rotyr)
summ.lmm.obs.16s.rotyr
write.csv(summ.lmm.obs.16s.rotyr$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_obs_16s_rotyr.csv")
a = Anova(lmm.obs.16s.rotyr, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_obs_16s_rotyears.csv")

#######################
# soil pH
######################
lmm.obs.16s.pH = lmer(mean ~ 1+soil_depth*period_numeric*soil_ph_value + (1+period_numeric|field), data = obs.16s)

plot(lmm.obs.16s.pH) # not bad
qqnorm(resid(lmm.obs.16s.pH)) # normality not bad
qqline(resid(lmm.obs.16s.pH))
hist(resid(lmm.obs.16s.pH)) # looks normal

(summ.lmm.obs.16s.pH = summary(lmm.obs.16s.pH))

write.csv(summ.lmm.obs.16s.pH$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_obs_16s_pH.csv")
a = Anova(lmm.obs.16s.pH, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_obs_16s_pH.csv")



# shannon
lmm.shan.16s = lmer(mean ~ 1+soil_depth*period_numeric*cul_type + (1+period_numeric|field), data = shan.16s)

plot(lmm.shan.16s) # heterogeneity seems fine
qqnorm(resid(lmm.shan.16s)) # normality seems fine
qqline(resid(lmm.shan.16s)) 
hist(resid(lmm.shan.16s)) # looks normal


Anova(lmm.shan.16s, type = 3, test.statistic = "F")
emmeans(lmm.shan.16s, pairwise ~ soil_depth)

summ.lmm.shan.16s = summary(lmm.shan.16s)
summ.lmm.shan.16s

write.csv(summ.lmm.shan.16s$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/coef_lmm_shan_16s_fieldasrandom.csv")


############
# soil types
############
lmm.shan.16s.stype = lmer(mean ~ 1+soil_depth*period_numeric*soil_compname + (1+period_numeric|field), data = shan.16s)

plot(lmm.shan.16s.stype)
qqnorm(resid(lmm.shan.16s.stype)) # normality not bad
qqline(resid(lmm.shan.16s.stype))

hist(resid(lmm.shan.16s.stype)) # looks normal

(summ.lmm.shan.16s.stype = summary(lmm.shan.16s.stype))
summ.lmm.shan.16s.stype
write.csv(summ.lmm.shan.16s.stype$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_shan_16s_soiltype.csv")
a = Anova(lmm.shan.16s.stype, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_shan_16s_soiltypes.csv")

emmeans(lmm.shan.16s.stype, pairwise ~ soil_compname)


##########
# potato intensity
##########


lmm.shan.16s.pota = lmer(mean ~ 1+soil_depth*period_numeric*percPota + (1+period_numeric|field), data = shan.16s)

plot(lmm.shan.16s.pota)
qqnorm(resid(lmm.shan.16s.pota)) # normality not bad
qqline(resid(lmm.shan.16s.pota))

hist(resid(lmm.shan.16s.pota)) # looks normal

(summ.lmm.shan.16s.pota = summary(lmm.shan.16s.pota))
summ.lmm.shan.16s.pota
write.csv(summ.lmm.shan.16s.pota$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_shan_16s_potato.csv")
a = Anova(lmm.shan.16s.pota, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_shan_16s_potato.csv")

# emmeans(lmm.obs.16s.pota, pairwise ~ percPota)

######################
# crop diversity
#####################
lmm.shan.16s.crops = lmer(mean ~ 1+soil_depth*period_numeric*totalCrops + (1+period_numeric|field), data = shan.16s)

plot(lmm.shan.16s.crops)
qqnorm(resid(lmm.shan.16s.crops)) # normality not bad
qqline(resid(lmm.shan.16s.crops))

hist(resid(lmm.shan.16s.crops)) # looks normal

(summ.lmm.shan.16s.crops = summary(lmm.shan.16s.crops))
summ.lmm.shan.16s.crops
write.csv(summ.lmm.shan.16s.crops$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_shan_16s_crops.csv")
a = Anova(lmm.shan.16s.crops, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_shan_16s_crop.csv")

######################
# total years
#####################
lmm.shan.16s.rotyr = lmer(mean ~ 1+soil_depth*period_numeric*totalRotYrs + (1+period_numeric|fplots), data = shan.16s)

plot(lmm.shan.16s.rotyr)
qqnorm(resid(lmm.shan.16s.rotyr)) # normality not bad
qqline(resid(lmm.shan.16s.rotyr))

hist(resid(lmm.shan.16s.rotyr)) # looks normal

summ.lmm.shan.16s.rotyr = summary(lmm.shan.16s.rotyr)
summ.lmm.shan.16s.rotyr
write.csv(summ.lmm.shan.16s.rotyr$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_shan_16s_rotyr.csv")
a = Anova(lmm.shan.16s.rotyr, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_shan_16s_rotyrs.csv")



#######################
# soil pH
######################
lmm.shan.16s.pH = lmer(mean ~ 1+soil_depth*period_numeric*soil_ph_value + (1+period_numeric|field), data = shan.16s)

plot(lmm.shan.16s.pH) # not bad
qqnorm(resid(lmm.shan.16s.pH)) # normality not bad
qqline(resid(lmm.shan.16s.pH))
hist(resid(lmm.shan.16s.pH)) # looks normal

(summ.lmm.shan.16s.pH = summary(lmm.shan.16s.pH))

write.csv(summ.lmm.shan.16s.pH$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_shan_16s_pH.csv")
a = Anova(lmm.shan.16s.pH, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_shan_16s_pH.csv")



# inverse simpson

# failed to converge if use period as factor: answer: https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer

# use period_numeric as numerical
lmm.invS.16s = lmer(mean ~ 1+soil_depth*period_numeric*cul_type + (1+period_numeric|field), data = invS.16s)
plot(lmm.invS.16s) # not bad
qqnorm(resid(lmm.invS.16s)) # not bad
qqline(resid(lmm.invS.16s)) # not bad

hist(resid(lmm.invS.16s)) # seems normal, a little skewness

Anova(lmm.invS.16s, type = 3, test.statistic = "F")


summ.lmm.invS.16s = summary(lmm.invS.16s)
summ.lmm.invS.16s
write.csv(summ.lmm.invS.16s$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/coef_lmm_invS_16s_fieldasrandom.csv")

############
# soil types
############
lmm.invS.16s.stype = lmer(mean ~ 1+soil_depth*period_numeric*soil_compname + (1+period_numeric|field), data = invS.16s)

plot(lmm.invS.16s.stype)
qqnorm(resid(lmm.invS.16s.stype)) # normality not bad
qqline(resid(lmm.invS.16s.stype))

hist(resid(lmm.invS.16s.stype)) # looks almost normal

(summ.lmm.invS.16s.stype = summary(lmm.invS.16s.stype))
summ.lmm.invS.16s.stype
write.csv(summ.lmm.invS.16s.stype$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_invS_16s_soiltype.csv")
a = Anova(lmm.invS.16s.stype, type = 3, test.statistic = "F")
# emmeans(lmm.invS.16s.stype, pairwise ~ soil_compname)
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_invS_16s_soiltypes.csv")

##########
# potato intensity
##########


lmm.invS.16s.pota = lmer(mean ~ 1+soil_depth*period_numeric*percPota + (1+period_numeric|field), data = invS.16s)

plot(lmm.invS.16s.pota)
qqnorm(resid(lmm.invS.16s.pota)) # normality not bad
qqline(resid(lmm.invS.16s.pota))

hist(resid(lmm.invS.16s.pota)) # looks normal

(summ.lmm.invS.16s.pota = summary(lmm.invS.16s.pota))
summ.lmm.invS.16s.pota
write.csv(summ.lmm.invS.16s.pota$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_invS_16s_potato.csv")
a = Anova(lmm.invS.16s.pota, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_invS_16s_pota.csv")
# emmeans(lmm.obs.16s.pota, pairwise ~ percPota)

######################
# crop diversity
#####################
lmm.invS.16s.crops = lmer(mean ~ 1+soil_depth*period_numeric*totalCrops + (1+period_numeric|field), data = invS.16s)

plot(lmm.invS.16s.crops)
qqnorm(resid(lmm.invS.16s.crops)) # normality not bad
qqline(resid(lmm.invS.16s.crops))

hist(resid(lmm.invS.16s.crops)) # a little skewness to the left

(summ.lmm.invS.16s.crops = summary(lmm.invS.16s.crops))
summ.lmm.invS.16s.crops
write.csv(summ.lmm.invS.16s.crops$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_invS_16s_crops.csv")
a = Anova(lmm.invS.16s.crops, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_invS_16s_crops.csv")
######################
# total years
#####################
lmm.invS.16s.rotyr = lmer(mean ~ 1+soil_depth*period_numeric*totalRotYrs + (1+period_numeric|fplots), data = invS.16s)

plot(lmm.invS.16s.rotyr)
qqnorm(resid(lmm.invS.16s.rotyr)) # normality not bad
qqline(resid(lmm.invS.16s.rotyr))

hist(resid(lmm.invS.16s.rotyr)) # looks normal

summ.lmm.invS.16s.rotyr = summary(lmm.invS.16s.rotyr)
summ.lmm.invS.16s.rotyr
write.csv(summ.lmm.invS.16s.rotyr$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_invS_16s_rotyr.csv")
a = Anova(lmm.invS.16s.rotyr, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_invS_16s_rotyrs.csv")




#######################
# soil pH
######################
lmm.invS.16s.pH = lmer(mean ~ 1+soil_depth*period_numeric*soil_ph_value + (1+period_numeric|field), data = invS.16s)

library(optimx)
lmm.invS.16s.pH = lmer(mean ~ 1+soil_depth*period_numeric*soil_ph_value + (1+period_numeric|fplots), data = invS.16s, REML = T, control = lmerControl(optimizer = "optimx", optCtrl = list(method="L-BFGS-B")))

plot(lmm.invS.16s.pH) # not bad
qqnorm(resid(lmm.invS.16s.pH)) # normality not bad
qqline(resid(lmm.invS.16s.pH))
hist(resid(lmm.invS.16s.pH)) # looks normal

(summ.lmm.invS.16s.pH = summary(lmm.invS.16s.pH))

write.csv(summ.lmm.invS.16s.pH$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_invS_16s_pH.csv")
a = Anova(lmm.invS.16s.pH, type = 3, test.statistic = "F")
a
write.csv(a, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_invS_16s_pH.csv")









#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# fungal
library(mgcv)
library(splines)
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

obs.its =  subset(alphaITS, measure == "Observed richness")
shan.its = subset(alphaITS, measure == "e^Shannon")
invS.its = subset(alphaITS, measure == "inverse Simpson")


plot(mean ~ period,     data = obs.its)
plot(mean ~ soil_depth, data = obs.its)
plot(mean ~ cul_type,   data = obs.its)

plot(mean ~ period,     data = shan.its)
plot(mean ~ soil_depth, data = shan.its)
plot(mean ~ cul_type,   data = shan.its)

plot(mean ~ period,     data = invS.its)
plot(mean ~ soil_depth, data = invS.its)
plot(mean ~ cul_type,   data = invS.its)



# lmm.obs.its = lmer(mean ~ 1+soil_depth*period_numeric*cul_type + (bs(period_numeric)|fplots), data = obs.its)
lmm.obs.its = lmer(mean ~ 1+soil_depth*period_numeric*cul_type + (1+period_numeric|field), data = obs.its)

plot(lmm.obs.its) # better with period numeric as a random intercept
qqnorm(resid(lmm.obs.its)) # a little skewed
qqline(resid(lmm.obs.its))
hist(resid(lmm.obs.its)) # a little skewed

summ.lmm.obs.its = summary(lmm.obs.its)
summ.lmm.obs.its
write.csv(summ.lmm.obs.its$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/coef_lmm_obs_its_fieldasrandom.csv")


Anova(lmm.obs.its, type = 3, test.statistic = "F")

emmeans(lmm.obs.its, pairwise ~ soil_depth)



############
# soil types
############
lmm.obs.its.stype = lmer(mean ~ 1+soil_depth*period_numeric*soil_compname + (1+period_numeric|field), data = obs.its)

plot(lmm.obs.its.stype)
qqnorm(resid(lmm.obs.its.stype)) # a little skewed
qqline(resid(lmm.obs.its.stype))

hist(resid(lmm.obs.its.stype)) # a little skewed

summ.lmm.obs.its.stype = summary(lmm.obs.its.stype)
summ.lmm.obs.its.stype
write.csv(summ.lmm.obs.its.stype$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_obs_its_soiltype.csv")
b = Anova(lmm.obs.its.stype, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_obs_its_soiltypes.csv")

emmeans(lmm.obs.its.stype, pairwise ~ soil_compname)


##########
# potato intensity
##########


lmm.obs.its.pota = lmer(mean ~ 1+soil_depth*period_numeric*percPota + (1+period_numeric|field), data = obs.its)

plot(lmm.obs.its.pota)
qqnorm(resid(lmm.obs.its.pota)) # normality not bad
qqline(resid(lmm.obs.its.pota))

hist(resid(lmm.obs.its.pota)) # looks almost normal

summ.lmm.obs.its.pota = summary(lmm.obs.its.pota)
summ.lmm.obs.its.pota
write.csv(summ.lmm.obs.its.pota$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_obs_its_potato.csv")
b = Anova(lmm.obs.its.pota, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_obs_its_pota.csv")
# emmeans(lmm.obs.16s.pota, pairwise ~ percPota)

######################
# crop diversity
#####################
lmm.obs.its.crops = lmer(mean ~ 1+soil_depth*period_numeric*totalCrops + (1+period_numeric|field), data = obs.its)

plot(lmm.obs.its.crops)
qqnorm(resid(lmm.obs.its.crops)) # normality not bad
qqline(resid(lmm.obs.its.crops)) # a little right skewed

hist(resid(lmm.obs.its.crops)) # looks normal

summ.lmm.obs.its.crops = summary(lmm.obs.its.crops)
summ.lmm.obs.its.crops
write.csv(summ.lmm.obs.its.crops$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_obs_its_crops.csv")
b = Anova(lmm.obs.its.crops, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_obs_its_crops.csv")
######################
# total years
#####################
lmm.obs.its.rotyr = lmer(mean ~ 1+soil_depth*period_numeric*totalRotYrs + (1+period_numeric|fplots), data = obs.its)

plot(lmm.obs.its.rotyr)
qqnorm(resid(lmm.obs.its.rotyr)) # normality not bad
qqline(resid(lmm.obs.its.rotyr))

hist(resid(lmm.obs.its.rotyr)) # looks normal

summ.lmm.obs.its.rotyr = summary(lmm.obs.its.rotyr)
summ.lmm.obs.its.rotyr
write.csv(summ.lmm.obs.its.rotyr$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_obs_its_rotyr.csv")
b = Anova(lmm.obs.its.rotyr, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_obs_its_rotyrs.csv")


######################
# pH
#####################
lmm.obs.its.pH = lmer(mean ~ 1+soil_depth*period_numeric*soil_ph_value+ (1+period_numeric|field), data = obs.its)

plot(lmm.obs.its.pH)
qqnorm(resid(lmm.obs.its.pH)) # normality not bad
qqline(resid(lmm.obs.its.pH))

hist(resid(lmm.obs.its.pH)) # looks normal

summ.lmm.obs.its.pH = summary(lmm.obs.its.pH)
summ.lmm.obs.its.pH
write.csv(summ.lmm.obs.its.pH$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_obs_its_pH.csv")
b = Anova(lmm.obs.its.pH, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_obs_its_pH.csv")



#||||||||||||||||||||||||||||||||||

# shannon
lmm.shan.its = lmer(mean ~ 1+soil_depth*period_numeric*cul_type + (1+period_numeric|field), data = shan.its)
plot(lmm.shan.its) # not bad
qqnorm(resid(lmm.shan.its)) # skewness, not quite normal
qqline(resid(lmm.shan.its))
hist(resid(lmm.shan.its))

summ.lmm.shan.its = summary(lmm.shan.its)
summ.lmm.shan.its
write.csv(summ.lmm.shan.its$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/coef_lmm_shan_its_fieldasrandom.csv")

Anova(lmm.shan.its, type = 3, test.statistic = "F")

emmeans(lmm.shan.its, pairwise ~ cul_type)


############
# soil types
############
lmm.shan.its.stype = lmer(mean ~ 1+soil_depth*period_numeric*soil_compname + (1+period_numeric|field), data = shan.its)

plot(lmm.shan.its.stype)
qqnorm(resid(lmm.shan.its.stype)) # skewed
qqline(resid(lmm.shan.its.stype))

hist(resid(lmm.shan.its.stype)) # a little skewed

summ.lmm.shan.its.stype = summary(lmm.shan.its.stype)
summ.lmm.shan.its.stype
write.csv(summ.lmm.shan.its.stype$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_shan_its_soiltype.csv")
b = Anova(lmm.shan.its.stype, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_shan_its_soiltypes.csv")
emmeans(lmm.shan.its.stype, pairwise ~ soil_compname)


##########
# potato intensity
##########


lmm.shan.its.pota = lmer(mean ~ 1+soil_depth*period_numeric*percPota + (1+period_numeric|field), data = shan.its)

plot(lmm.shan.its.pota) # not bad
qqnorm(resid(lmm.shan.its.pota)) # skewed
qqline(resid(lmm.shan.its.pota))

hist(resid(lmm.shan.its.pota)) # a little skewed

summ.lmm.shan.its.pota = summary(lmm.shan.its.pota)
summ.lmm.shan.its.pota
write.csv(summ.lmm.shan.its.pota$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_shan_its_potato.csv")
b = Anova(lmm.shan.its.pota, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_shan_its_pota.csv")
# emmeans(lmm.obs.16s.pota, pairwise ~ percPota)

######################
# crop diversity
#####################
lmm.shan.its.crops = lmer(mean ~ 1+soil_depth*period_numeric*totalCrops + (1+period_numeric|field), data = shan.its) # failed to converge with default Nelder-Mead


lmm.shan.its.crops = lmer(mean ~ 1+soil_depth*period_numeric*totalCrops + (1+period_numeric|fplots), data = shan.its, REML = T, control = lmerControl(optimizer = "optimx", optCtrl = list(method="L-BFGS-B")))


plot(lmm.shan.its.crops)
qqnorm(resid(lmm.shan.its.crops)) #
qqline(resid(lmm.shan.its.crops)) # a little right skewed

hist(resid(lmm.shan.its.crops)) # a little skewed

summ.lmm.shan.its.crops = summary(lmm.shan.its.crops)
summ.lmm.shan.its.crops
write.csv(summ.lmm.shan.its.crops$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_shan_its_crops.csv")
b = Anova(lmm.shan.its.crops, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_shan_its_crops.csv")

######################
# total years
#####################
lmm.shan.its.rotyr = lmer(mean ~ 1+soil_depth*period_numeric*totalRotYrs + (1+period_numeric|fplots), data = shan.its)

plot(lmm.shan.its.rotyr)
qqnorm(resid(lmm.shan.its.rotyr)) # skewed
qqline(resid(lmm.shan.its.rotyr))

hist(resid(lmm.shan.its.rotyr)) # a little skewed

summ.lmm.shan.its.rotyr = summary(lmm.shan.its.rotyr)
summ.lmm.shan.its.rotyr
write.csv(summ.lmm.shan.its.rotyr$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_shan_its_rotyr.csv")
b = Anova(lmm.shan.its.rotyr, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_shan_its_rotyrs.csv")

############
# soil pH
###########

lmm.shan.its.pH = lmer(mean ~ 1+soil_depth*period_numeric*soil_ph_value + (1+period_numeric|field), data = shan.its)

plot(lmm.shan.its.pH)
qqnorm(resid(lmm.shan.its.pH)) # skewed
qqline(resid(lmm.shan.its.pH))

hist(resid(lmm.shan.its.pH)) # a little skewed

summ.lmm.shan.its.pH = summary(lmm.shan.its.pH)
summ.lmm.shan.its.pH
write.csv(summ.lmm.shan.its.pH$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_shan_its_pH.csv")
b = Anova(lmm.shan.its.pH, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_shan_its_pH.csv")



#|||||||||||||||||||||||||||||


# invSimpson
# lmm.invS.its = lmer(mean ~ 1+soil_depth*period_numeric*cul_type + (bs(period_numeric)|fplots), data = invS.its)
lmm.invS.its = lmer(mean ~ 1+soil_depth*period_numeric*cul_type + (1+period_numeric|field), data = invS.its)
plot(lmm.invS.its) # not homogeneous, a trend observed
qqnorm(resid(lmm.invS.its)) # normality seems fine
qqline(resid(lmm.invS.its))
hist(resid(lmm.invS.its))



summ.lmm.invS.its = summary(lmm.invS.its)
summ.lmm.invS.its
write.csv(summ.lmm.invS.its$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/coef_lmm_invS_its_fieldasrandom.csv")


Anova(lmm.invS.its, type = 3, test.statistic = "F")

emmeans(lmm.invS.its, pairwise ~ cul_type)

############
# soil types
############
lmm.invS.its.stype = lmer(mean ~ 1+soil_depth*period_numeric*soil_compname + (1+period_numeric|field), data = invS.its)

plot(lmm.invS.its.stype) # seems not homogeneous
qqnorm(resid(lmm.invS.its.stype)) # a little skewed
qqline(resid(lmm.invS.its.stype))

hist(resid(lmm.invS.its.stype)) # a little skewed

summ.lmm.invS.its.stype = summary(lmm.invS.its.stype)
summ.lmm.invS.its.stype
write.csv(summ.lmm.invS.its.stype$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_invS_its_soiltype.csv")
b = Anova(lmm.invS.its.stype, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_invS_its_soiltypes.csv")
emmeans(lmm.invS.its.stype, pairwise ~ soil_compname)


##########
# potato intensity
##########


lmm.invS.its.pota = lmer(mean ~ 1+soil_depth*period_numeric*percPota + (1+period_numeric|field), data = invS.its)

plot(lmm.invS.its.pota) # not bad
qqnorm(resid(lmm.invS.its.pota)) # not bad
qqline(resid(lmm.invS.its.pota))

hist(resid(lmm.invS.its.pota)) # not bad

summ.lmm.invS.its.pota = summary(lmm.invS.its.pota)
summ.lmm.invS.its.pota
write.csv(summ.lmm.invS.its.pota$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_invS_its_potato.csv")
b = Anova(lmm.invS.its.pota, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_invS_its_pota.csv")
# emmeans(lmm.obs.16s.pota, pairwise ~ percPota)

######################
# crop diversity
#####################
lmm.invS.its.crops = lmer(mean ~ 1+soil_depth*period_numeric*totalCrops + (1+period_numeric|field), data = invS.its) 


lmm.invS.its.crops = lmer(mean ~ 1+soil_depth*period_numeric*totalCrops + (1+period_numeric|fplots), data = invS.its, REML = T, control = lmerControl(optimizer = "optimx", optCtrl = list(method="L-BFGS-B")))


plot(lmm.invS.its.crops) # a little heterogeneous
qqnorm(resid(lmm.invS.its.crops)) #
qqline(resid(lmm.invS.its.crops)) # a little skewness

hist(resid(lmm.invS.its.crops)) # not bad

(summ.lmm.invS.its.crops = summary(lmm.invS.its.crops))
summ.lmm.invS.its.crops
write.csv(summ.lmm.invS.its.crops$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_invS_its_crops.csv")
b = Anova(lmm.invS.its.crops, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_invS_its_crops.csv")

######################
# total years
#####################
lmm.invS.its.rotyr = lmer(mean ~ 1+soil_depth*period_numeric*totalRotYrs + (1+period_numeric|fplots), data = invS.its)

plot(lmm.invS.its.rotyr)
qqnorm(resid(lmm.invS.its.rotyr)) # a little skewed
qqline(resid(lmm.invS.its.rotyr))

hist(resid(lmm.invS.its.rotyr)) # not bad

summ.lmm.invS.its.rotyr = summary(lmm.invS.its.rotyr)
summ.lmm.invS.its.rotyr
write.csv(summ.lmm.invS.its.rotyr$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_invS_its_rotyr.csv")
b = Anova(lmm.invS.its.rotyr, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_invS_its_rotyrs.csv")


######################
# soil ph
#####################
lmm.invS.its.pH = lmer(mean ~ 1+soil_depth*period_numeric*soil_ph_value + (1+period_numeric|field), data = invS.its)

plot(lmm.invS.its.pH)
qqnorm(resid(lmm.invS.its.pH)) # a little skewed
qqline(resid(lmm.invS.its.pH))

hist(resid(lmm.invS.its.pH)) # not bad

summ.lmm.invS.its.pH = summary(lmm.invS.its.pH)
summ.lmm.invS.its.pH
write.csv(summ.lmm.invS.its.pH$coefficients, "./microcosm_soil_analysis_output/alpha_diversity/other factors/coef_lmm_invS_its_pH.csv")
b = Anova(lmm.invS.its.pH, type = 3, test.statistic = "F")
b
write.csv(b, "./microcosm_soil_analysis_output/alpha_diversity/other factors/anova_invS_its_pH.csv")



#########################################################################################

# Beta diversity

rm(list=ls())
source("./tools/plotCAP2.R")
dir.create("./microcosm_soil_analysis_output/beta_diversity")
# load pruned data
mydata.16s.pruned = readRDS("./microcosm_soil_analysis_output/data/mydata_16s_pruned.rds")
mydata.its.pruned = readRDS("./microcosm_soil_analysis_output/data/mydata_its_pruned.rds")

mydata16s.norm = mydata.16s.pruned %>%
        transform_sample_counts(function(x)x/sum(x))

mydataits.norm = mydata.its.pruned %>%
        transform_sample_counts(function(x)x/sum(x))



otu.16s = data.frame(t(otu_table(mydata16s.norm)))
otu.its = data.frame(t(otu_table(mydataits.norm)))

# convert soil depth to numeric
env.16s = data.frame(sample_data(mydata16s.norm))
env.its = data.frame(sample_data(mydataits.norm))

env.16s$soil_depth_numeric = as.numeric(gsub("cm", "", env.16s$soil_depth))
env.its$soil_depth_numeric = as.numeric(gsub("cm", "", env.its$soil_depth))
#########################################################################################


######################################
# CAP
###################################

#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# Bacteria
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

colnames(env.16s)


cap.16s = capscale(otu.16s ~cul_type+soil_compname+totalCrops+totalRotYrs+percPota + period_numeric + soil_depth_numeric + soil_ph_value, data = env.16s, distance = "bray", add = T)

vif.cca(cap.16s) # totalRotYrs and cultype > 10, so remove totalRotYrs

cap.16s = capscale(otu.16s ~cul_type+soil_compname+totalCrops+percPota+period_numeric+soil_depth_numeric+soil_ph_value, data = env.16s, distance = "bray", add=T)
vif.cca(cap.16s)

cap.16s

# check this model
set.seed(123)
anova(cap.16s, permutations = 1000) # 0.001
anova(cap.16s, by="terms", permutations = 1000) # all significant
anova(cap.16s, by="axis", permutations = 1000) # all 10

# # get parsimony model
# # by AIC
cap.16s.select = ordistep(capscale(otu.16s ~ 1, data = env.16s, distance = "bray", add=T), scope = formula(cap.16s), direction = "forward", permutations = 1000)
# save file
# saveRDS(cap.16s.select, "./microcosm_soil_analysis_output/beta_diversity/cap16s_AIC_selected.rds")

cap.16s.select = readRDS("./microcosm_soil_analysis_output/beta_diversity/cap16s_AIC_selected.rds")
cap.16s.select
cap.16s.select.r2 = readRDS("./microcosm_soil_analysis_output/beta_diversity/cap16s_R2_selected.rds")
cap.16s.select.r2

# adjust p value
cap.16s.select.padj = cap.16s.select
cap.16s.select.padj$anova$`Pr(>F)` = p.adjust(cap.16s.select$anova$`Pr(>F)`, method = "BH")
cap.16s.select.padj$anova
cap.16s.select.padj

# # by R2 
cap.16s.select.r2 = ordiR2step(capscale(otu.16s~1, data=env.16s, distance = "bray", add=T), scope = formula(cap.16s), direction = "forward", permutations = 1000, R2permutations = 1000)
# cap.16s.select.r2  otu.16s ~ cul_type + percPota + totalCrops + soil_compname + 
# soil_ph_value + period_numeric

# save file
# saveRDS(cap.16s.select.r2, "./microcosm_soil_analysis_output/beta_diversity/cap16s_R2_selected.rds")


cap.16s.select.r2.padj = cap.16s.select.r2
cap.16s.select.r2.padj$anova$`Pr(>F)` = p.adjust(cap.16s.select.r2$anova$`Pr(>F)`, method = "BH")
cap.16s.select.r2.padj$anova

cap.16s.final = capscale(otu.16s ~ cul_type + percPota + totalCrops + soil_compname + soil_ph_value +
                                 period_numeric, data = env.16s, distance="bray", add=T)

# saveRDS(cap.16s.final, "./microcosm_soil_analysis_output/beta_diversity/cap_16s_final.rds")


#------------------------------
### plot
#------------------------------
cap.16s.final = readRDS("./microcosm_soil_analysis_output/beta_diversity/cap_16s_final.rds")

set.seed(104)
anova(cap.16s.final, permutations = 1000, step = 1000)
anova(cap.16s.final, by="terms", permutations = 1000, step = 1000)
anova(cap.16s.final, by="axis", permutations = 1000, step = 1000)

# cap plot

taxcols.16s = as.data.frame(tax_table(mydata16s.norm))
taxcols.16s$ASV = rownames(taxcols.16s)

svg("./microcosm_soil_analysis_output/beta_diversity/cap16s_cultype.svg", width = 10, height = 12)
(graph.cap.16s = plotCAP(cap.16s.final, env = env.16s, tax_cols = taxcols.16s, display = 3, splitDisplay = F, sa.xadj = 0.2, sa.yadj = - 0.5, sa.colorSamples = T, sa.sample_group = "cul_type", sa.showCentroidGroup = "cul_type", sa.biplot = T,sa.ellipse = T, sa.addSp = T, sa.labelsp = T, sa.tax_level = "Genus")+ scale_color_discrete(labels = c("Conventional", "Mixed", "Organic"))+labs(color = "Cropping systems"))
dev.off()

svg("./microcosm_soil_analysis_output/beta_diversity/cap16s_soiltype.svg", width = 10, height = 12)
(graph.cap.16.soiltype = plotCAP(cap.16s.final, env = env.16s, tax_cols = taxcols.16s, display = 1, splitDisplay = F, sa.xadj = 0.2, sa.yadj = -0.5, sa.colorSamples = T, sa.sample_group = "soil_compname", sa.showCentroidGroup = "soil_compname", sa.biplot = T,sa.ellipse = T))
dev.off()



#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# Fungi
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

cap.its = capscale(otu.its ~cul_type+soil_compname+totalCrops+totalRotYrs+percPota + period_numeric + soil_depth_numeric + soil_ph_value, data = env.its, distance = "bray", add = T)
cap.its

vif.cca(cap.its) # colinearity between cul_type and totalRotYrs


cap.its = capscale(otu.its ~cul_type+soil_compname+totalCrops+percPota+period_numeric + soil_depth_numeric + soil_ph_value, data = env.its, distance = "bray", add=T)
cap.its

vif(cap.its)

set.seed(111)
anova(cap.its, permutations = 1000) # sig
anova(cap.its, by="terms", permutations = 1000) # soil depth not significant
anova(cap.its, by="axis", permutations = 1000) # top 9


# select by AIC
cap.its.select = ordistep(capscale(otu.its ~ 1, data=env.its, distance = "bray", add=T), scope = formula(cap.its), direction = "forward", permutations = 1000)
cap.its.select #formula = otu.its ~ cul_type + percPota + soil_compname + period_numeric + totalCrops +
# soil_ph_value + soil_depth_numeric

cap.its.select = readRDS("./microcosm_soil_analysis_output/beta_diversity/capITS_AIC_select.rds")
cap.its.select # otu.its ~ cul_type + percPota + soil_compname + period_numeric +      totalCrops + soil_ph_value + soil_depth_numeric
saveRDS(cap.its.select, "./microcosm_soil_analysis_output/beta_diversity/capITS_AIC_select.rds")
cap.its.select.padj = cap.its.select
cap.its.select.padj$anova$`Pr(>F)` = p.adjust(cap.its.select$anova$`Pr(>F)`, method = "BH")
cap.its.select.padj$anova # soil_compname borderline --- go use R2 selection instead



# select by R2
cap.its.select.r2 = ordiR2step(capscale(otu.its~1, data=env.its, distance = "bray", add=T), scope = formula(cap.its), direction="forward", permutations = 1000, R2permutations = 1000)
cap.its.select.r2

cap.its.select.r2 = readRDS("./microcosm_soil_analysis_output/beta_diversity/capITS_R2_select.rds")
cap.its.select.r2 #  otu.its ~ cul_type + percPota + soil_compname + period_numeric + totalCrops + soil_ph_value + soil_depth_numeric
saveRDS(cap.its.select.r2, "./microcosm_soil_analysis_output/beta_diversity/capITS_R2_select.rds")

cap.its.select.r2.padj = cap.its.select.r2
cap.its.select.r2.padj$anova$`Pr(>F)` = p.adjust(cap.its.select.r2$anova$`Pr(>F)`, method = "BH")
cap.its.select.r2.padj$anova

cap.its.final = capscale(otu.its ~ cul_type + 
                                 percPota + 
                                 soil_compname + 
                                 period_numeric + 
                                 totalCrops +
                                 soil_ph_value + 
                                 soil_depth_numeric, data = env.its, distance = "bray", add=T)
cap.its.final

# saveRDS(cap.its.final, "./microcosm_soil_analysis_output/beta_diversity/cap_ITS_final.rds")

set.seed(10000)
anova(cap.its.final, permutations = 1000)
anova(cap.its.final, by="terms", permutations = 1000) # soil depth not significant
anova(cap.its.final, by="axis", permutations = 1000)


# remove soil depth
cap.its.final = capscale(otu.its ~ cul_type + 
                                 percPota + 
                                 soil_compname + 
                                 period_numeric + 
                                 totalCrops +
                                 soil_ph_value , data = env.its, distance = "bray", add=T)
cap.its.final
set.seed(10001)
anova(cap.its.final, permutations = 1000)
anova(cap.its.final, by="terms", permutations = 1000) # 
anova(cap.its.final, by="axis", permutations = 1000) # all 9 caps significant


#------------------------------
### plot
#------------------------------

cap.its.final = readRDS("./microcosm_soil_analysis_output/beta_diversity/cap_ITS_final.rds")
# cap plot

taxcols.its = as.data.frame(tax_table(mydataits.norm))
taxcols.its$ASV = rownames(taxcols.its)

svg("./microcosm_soil_analysis_output/beta_diversity/capITS_cultype.svg", width = 10, height = 12)
(graph.cap.its = plotCAP(cap.its.final, env = env.its, tax_cols = taxcols.its, display = 3, splitDisplay = F, sa.xadj = 0.3, sa.yadj = -0.5,sa.colorSamples = T, sa.sample_group = "cul_type", sa.showCentroidGroup = "cul_type", sa.biplot = T, sa.nudegy = 0.3, sa.nudgex = 0.3, sa.ellipse = T, sa.addSp = T, sa.labelsp = T, sa.tax_level = "Genus") + scale_color_discrete(labels = c("Conventional", "Mixed", "Organic"))+labs(color="Cropping systems"))
dev.off()



# svg("./microcosm_soil_analysis_output/beta_diversity/capITS_cultype_samples.svg", width = 10, height = 12)
# (graph.cap.its = plotCAP(cap.its.final, env = env.its, tax_cols = taxcols.its, display = 3, splitDisplay = F, sa.xadj = 0.3, sa.yadj = -0.5,sa.colorSamples = T, sa.sample_group = "cul_type", sa.showCentroidGroup = "cul_type", sa.biplot = T, sa.nudegy = 0.3, sa.nudgex = 0.3, sa.ellipse = T, sa.labelsa = T, sa.labelby = "field", sa.addSp = F, sa.labelsp = F, sa.tax_level = "Genus"))
# dev.off()


svg("./microcosm_soil_analysis_output/beta_diversity/capITS_soiltype.svg", width = 13, height = 14)
(graph.cap.it.soiltype = plotCAP(cap.its.final, env = env.its, tax_cols = taxcols.its, display = 1, splitDisplay = F, sa.xadj = 0.2, sa.yadj = -0.5, sa.colorSamples = T, sa.sample_group = "soil_compname", sa.showCentroidGroup = "soil_compname", sa.biplot = T,sa.ellipse = T))
dev.off()


tiff("./microcosm_soil_analysis_output/beta_diversity/cap_all.tif", width = 2100, height = 2200, compression = "lzw")
ggpubr::ggarrange(graph.cap.16s, graph.cap.its, ncol=2, labels = c("A", "B"), align = "hv", common.legend = T, legend = "bottom")
dev.off()




svg("./microcosm_soil_analysis_output/beta_diversity/cap_all.svg", width = 18, height = 10)
ggpubr::ggarrange(graph.cap.16s, graph.cap.its, ncol=2, labels = c("A", "B"), align = "hv", common.legend = T, legend = "bottom")
dev.off()


svg("./microcosm_soil_analysis_output/beta_diversity/cap_all_soiltypes.svg", width = 18, height = 16)
ggpubr::ggarrange(graph.cap.16.soiltype, graph.cap.it.soiltype, ncol=2, labels = c("A", "B"), align = "hv", common.legend = T, legend = "bottom")
dev.off()


################################

# diversity distance vs time

################################



bray.16s = vegdist(otu.16s)

str(bray.16s)

names_bray_16s = attr(bray.16s, 'Labels')
names(bray.16s) = names_bray_16s

bray16s.mtx = as.matrix(bray.16s)

bray16s.mtx[1:5, 1:5]

anno_col.16s = data.frame(Managements=env.16s$cul_type, Soil_depth = env.16s$soil_depth, Time = env.16s$period)
base::rownames(anno_col.16s) = colnames(bray16s.mtx)


pheatmap::pheatmap(bray16s.mtx, annotation_col = anno_col.16s, cluster_cols =  F)



bray.its = vegdist(otu.its)

str(bray.its)

names_bray_its = attr(bray.its, 'Labels')
names(bray.its) = names_bray_its

brayits.mtx = as.matrix(bray.its)

brayits.mtx[1:5, 1:5]

anno_col.its = data.frame(Managements=env.its$cul_type)
base::rownames(anno_col.its) = colnames(brayits.mtx)

anno_row.its = data.frame(Time = env.its$period)
base::rownames(anno_row.its) = base::rownames(brayits.mtx)

pheatmap::pheatmap(brayits.mtx, annotation_col = anno_col.its, annotation_row = anno_row.its, fontsize_row = 5, fontsize_col = 5)


#########################################################################################################

# PERMANOVA

#########################################################################################################

set.seed(2)
adon.16s = adonis2( otu.16s ~ cul_type + soil_compname + totalCrops + percPota + period_numeric + soil_ph_value, data = env.16s, method = "bray", permutations = 10000)

# adon.16s = adonis2( otu.16s ~ cul_type + soil_compname + totalCrops + percPota + period_numeric + soil_ph_value, data = env.16s, method = "bray", permutations = 10000, by="margin")

adon.16s

# saveRDS(adon.16s, "./microcosm_soil_analysis_output/beta_diversity/adon_16s.rds")



set.seed(2)
disp16s.cul = betadisper(bray.16s, env.16s$cul_type, add = T, bias.adjust = T)
permutest(disp16s.cul, pairwise = T, permutations = 10000) # 0.001

set.seed(2)
disp16s.compname = betadisper(bray.16s, env.16s$soil_compname, add = T, bias.adjust = T)
permutest(disp16s.compname, pairwise = T, permutations = 10000) # 0.001

set.seed(2)
disp16s.period = betadisper(bray.16s, env.16s$period_numeric, add = T, bias.adjust = T)
permutest(disp16s.period, pairwise = T, permutations = 1000) # 0.002
# 
# disp16s.period = betadisper(bray.16s, env.16s$period)
# permutest(disp16s.period, pairwise = T) # 0.001
# 
# disp16s.totalCrops = betadisper(bray.16s, env.16s$totalCrops)
# permutest(disp16s.totalCrops, pairwise = T) # 0.001
# 
# disp16s.percPota = betadisper(bray.16s, env.16s$percPota)
# permutest(disp16s.percPota, pairwise = T) # 0.001

set.seed(3)
adon.its = adonis2( otu.its ~ cul_type + soil_compname + totalCrops + percPota + period_numeric + soil_ph_value, data = env.its, method = "bray", permutations = 10000)

adon.its

saveRDS(adon.its, "./microcosm_soil_analysis_output/beta_diversity/adon_its.rds")


set.seed(3)
dispits.cul = betadisper(bray.its, env.its$cul_type, add = T, bias.adjust = T)
permutest(dispits.cul, pairwise = T, permutations = 10000) # 0.001

set.seed(3)
dispits.compname = betadisper(bray.its, env.its$soil_compname, add = T, bias.adjust = T)
permutest(dispits.compname, pairwise = T, permutations = 10000) # 0.001


set.seed(3)
dispits.period = betadisper(bray.its, env.its$period_numeric, add = T, bias.adjust = T)
permutest(dispits.period, pairwise = T, permutations = 1000) # 0.001


###############################################################################################


# compositional analysis

rm(list= ls())

mydata.16s.pruned = readRDS("./microcosm_soil_analysis_output/data/mydata_16s_pruned.rds")
mydata.its.pruned = readRDS("./microcosm_soil_analysis_output/data/mydata_its_pruned.rds")

dir.create("./microcosm_soil_analysis_output/taxonomy_composition")

source("./tools/plotCompHeatmap.R")
source("./tools/topAbundTaxa.R")
###############################################################################################

chooseCol = rev(brewer.pal(9, "YlGnBu"))

meta16s = data.frame(sample_data(mydata.16s.pruned))

levels(meta16s$cul_type) = c("Conventional", "Mixed", "Organic")
levels(meta16s$cul_type)

hmap16s_phy = plotCompHeatMap(mydata.16s.pruned, 
                              taxa_rank = "Phylum",
                              userMeta = T,
                              metaInput = meta16s,
                              show_n = 15, 
                              hp_colr_alpha = 0.3,
                              hp_colr = chooseCol,
                              anno_list = list(column=c("period", "cul_type")),
                              anno_names = list(cols=c("Time", "Cropping systems")),
                              use_hclust = "bray",
                              hclust_pos = "col",
                              colsplit = "soil_depth",
                              dend_col = F,
                              dend_row = F,
                              name = "Relative abundance (%)",
                              change_colnames = "seq",
                              column_names_gp = gpar(fontsize=8),
                              row_names_gp = gpar(fontsize=10),
                              row_names_side = "left",
                              show_column_names = F,
                              heatmap_legend_param = list(direction = "horizontal")
                        )

saveRDS(hmap16s_phy, "./microcosm_soil_analysis_output/taxonomy_composition/hmap16s_phy.rds")

svg("./microcosm_soil_analysis_output/taxonomy_composition/heatmap_16s_phylum.svg", width = 12, height = 8)
draw(hmap16s_phy, heatmap_legend_side = "bottom")
dev.off()

hmap16s_gen = plotCompHeatMap(mydata.16s.pruned, 
                              taxa_rank = "Genus",
                              userMeta = T,
                              metaInput = meta16s,
                              show_n = 15,
                              hp_colr_alpha = 0.3,
                              hp_colr = chooseCol,
                              anno_list = list(column=c("period", "cul_type")),
                              anno_names = list(cols=c("Time", "Cropping systems")),
                              use_hclust = "bray",
                              hclust_pos = "col",
                              colsplit = "soil_depth",
                              dend_col = F,
                              dend_row = F,
                              name = "Relative abundance (%)",
                              change_colnames = "seq",
                              column_names_gp = gpar(fontsize=8),
                              row_names_gp = gpar(fontsize=10),
                              row_names_side = "left",
                              show_column_names = F,
                              heatmap_legend_param = list(direction = "horizontal"))

saveRDS(hmap16s_gen, "./microcosm_soil_analysis_output/taxonomy_composition/hmap16s_gen.rds")

svg("./microcosm_soil_analysis_output/taxonomy_composition/heatmap_16s_gen.svg", width = 12, height = 8)
draw(hmap16s_gen, heatmap_legend_side = "bottom")
dev.off()

# ITS

metaITS = data.frame(sample_data(mydata.its.pruned))

levels(metaITS$cul_type) = c("Conventional", "Mixed", "Organic")
levels(metaITS$cul_type)



hmapits_phy = plotCompHeatMap(mydata.its.pruned, 
                              taxa_rank = "Phylum", 
                              userMeta = T,
                              metaInput = metaITS,
                              show_n = 10,
                              hp_colr_alpha = 0.3,
                              hp_colr = chooseCol,
                              anno_list = list(column=c("period", "cul_type")),
                              anno_names = list(cols=c("Time", "Cropping systems")),
                              use_hclust = "bray",
                              hclust_pos = "col",
                              colsplit = "soil_depth",
                              dend_col = F,
                              dend_row = F,
                              name = "Relative abundance (%)",
                              change_colnames = "seq",
                              column_names_gp = gpar(fontsize=8),
                              row_names_gp = gpar(fontsize=10),
                              row_names_side = "left",
                              show_column_names = F,
                              heatmap_legend_param = list(direction = "horizontal"))

saveRDS(hmapits_phy, "./microcosm_soil_analysis_output/taxonomy_composition/hmapits_phy.rds")

svg("./microcosm_soil_analysis_output/taxonomy_composition/heatmap_its_phylum.svg", width = 12, height = 8)
draw(hmapits_phy, heatmap_legend_side = "bottom")
dev.off()

hmapits_gen = plotCompHeatMap(mydata.its.pruned, 
                              taxa_rank = "Genus",
                              userMeta = T,
                              metaInput = metaITS,
                              show_n = 20, 
                              hp_colr_alpha = 0.5,
                              hp_colr = chooseCol,
                              anno_list = list(column=c("period", "cul_type")),
                              anno_names = list(cols=c("Time", "Cropping systems")),
                              use_hclust = "bray",
                              hclust_pos = "col",
                              colsplit = "soil_depth",
                              dend_col = F,
                              dend_row = F,
                              name = "Relative abundance (%)",
                              change_colnames = "seq",
                              column_names_gp = gpar(fontsize=8),
                              row_names_gp = gpar(fontsize=10),
                              row_names_side = "left",
                              show_column_names = F,
                              heatmap_legend_param = list(direction = "horizontal"))

saveRDS(hmapits_gen, "./microcosm_soil_analysis_output/taxonomy_composition/hmapits_gen.rds")

svg("./microcosm_soil_analysis_output/taxonomy_composition/heatmap_its_gen.svg", width = 12, height = 8)
draw(hmapits_gen, heatmap_legend_side = "bottom")
dev.off()
####################################################

# top phyla in each time in each cultural practices

####################################################




#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# top 10 bacteria
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

showTopN(mydata.16s.pruned, n = 10, taxlevel = "Phylum")

# Phylum           avg_rel_abd
# <fct>                  <dbl>
#         1 Proteobacteria       0.340  
# 2 Actinobacteria       0.225  
# 3 Acidobacteria        0.137  
# 4 Chloroflexi          0.0761 
# 5 Firmicutes           0.0634 
# 6 Gemmatimonadetes     0.0533 
# 7 Bacteroidetes        0.0453 
# 8 Verrucomicrobia      0.0182 
# 9 Planctomycetes       0.0109 
# 10 Nitrospirae          0.00865


showTopN(mydata.16s.pruned, n = 10, taxlevel = "Genus")


# Selecting by avg_rel_abd
# # A tibble: 10 x 2
# Genus             avg_rel_abd
# <fct>                   <dbl>
#         1 Pseudarthrobacter      0.0442
# 2 c__Subgroup 6          0.0396
# 3 Bacillus               0.0308
# 4 Sphingomonas           0.0281
# 5 Pseudomonas            0.0199
# 6 RB41                   0.0196
# 7 Nocardioides           0.0168
# 8 Streptomyces           0.0130
# 9 c__KD4-96              0.0114
# 10 Massilia               0.0106


# conventional

my16s.C.t0 = subset_samples(mydata.16s.pruned, cul_type == "Conventional" & period == "0wk")
my16s.C.t1 = subset_samples(mydata.16s.pruned, cul_type == "Conventional" & period == "1wk")
my16s.C.t3 = subset_samples(mydata.16s.pruned, cul_type == "Conventional" & period == "3wk")
my16s.C.t6 = subset_samples(mydata.16s.pruned, cul_type == "Conventional" & period == "6wk")

my16s.O.t0 = subset_samples(mydata.16s.pruned, cul_type == "Organic" & period == "0wk")
my16s.O.t1 = subset_samples(mydata.16s.pruned, cul_type == "Organic" & period == "1wk")
my16s.O.t3 = subset_samples(mydata.16s.pruned, cul_type == "Organic" & period == "3wk")
my16s.O.t6 = subset_samples(mydata.16s.pruned, cul_type == "Organic" & period == "6wk")


my16s.G.t0 = subset_samples(mydata.16s.pruned, cul_type == "Green" & period == "0wk")
my16s.G.t1 = subset_samples(mydata.16s.pruned, cul_type == "Green" & period == "1wk")
my16s.G.t3 = subset_samples(mydata.16s.pruned, cul_type == "Green" & period == "3wk")
my16s.G.t6 = subset_samples(mydata.16s.pruned, cul_type == "Green" & period == "6wk")


phy_Ct0 = showTopN(my16s.C.t0, 10, "Phylum")
phy_Ct1 = showTopN(my16s.C.t1, 10, "Phylum")
phy_Ct3 = showTopN(my16s.C.t3, 10, "Phylum")
phy_Ct6 = showTopN(my16s.C.t6, 10, "Phylum")

phy_Ot0 = showTopN(my16s.O.t0, 10, "Phylum")
phy_Ot1 = showTopN(my16s.O.t1, 10, "Phylum")
phy_Ot3 = showTopN(my16s.O.t3, 10, "Phylum")
phy_Ot6 = showTopN(my16s.O.t6, 10, "Phylum")


phy_Gt0 = showTopN(my16s.G.t0, 10, "Phylum")
phy_Gt1 = showTopN(my16s.G.t1, 10, "Phylum")
phy_Gt3 = showTopN(my16s.G.t3, 10, "Phylum")
phy_Gt6 = showTopN(my16s.G.t6, 10, "Phylum")


phyC_top10 = data.frame(Type = rep("Conventional", 
                                        nrow(phy_Ct0)+
                                        nrow(phy_Ct1)+
                                        nrow(phy_Ct3)+
                                        nrow(phy_Ct6)),
                        Time = c(rep("T0", nrow(phy_Ct0)), 
                                 rep("T1", nrow(phy_Ct1)), 
                                 rep("T3", nrow(phy_Ct3)), 
                                 rep("T6", nrow(phy_Ct6))))

top10_phyC = cbind(phyC_top10, rbind(phy_Ct0, 
                                     phy_Ct1, 
                                     phy_Ct3, 
                                     phy_Ct6))

# organic
phyO_top10 = data.frame(Type = rep("Organic", nrow(phy_Ot0)+
                                                   nrow(phy_Ot1)+
                                                   nrow(phy_Ot3)+
                                                   nrow(phy_Ot6)),
                        Time = c(rep("T0", nrow(phy_Ot0)), 
                                 rep("T1", nrow(phy_Ot1)), 
                                 rep("T3", nrow(phy_Ot3)), 
                                 rep("T6", nrow(phy_Ot6))))

top10_phyO = cbind(phyO_top10, rbind(phy_Ot0, 
                                     phy_Ot1, 
                                     phy_Ot3, 
                                     phy_Ot6))


# green
phyG_top10 = data.frame(Type = rep("Green", 
                                           nrow(phy_Gt0)+
                                           nrow(phy_Gt1)+
                                           nrow(phy_Gt3)+
                                           nrow(phy_Gt6)),
                        Time = c(rep("T0", nrow(phy_Gt0)), 
                                 rep("T1", nrow(phy_Gt1)), 
                                 rep("T3", nrow(phy_Gt3)), 
                                 rep("T6", nrow(phy_Gt6))))

top10_phyG = cbind(phyG_top10, rbind(phy_Gt0, 
                                     phy_Gt1, 
                                     phy_Gt3, 
                                     phy_Gt6))


# top10 16s in Phylum along time

top10_bytime_16s = rbind(top10_phyC, top10_phyO, top10_phyG)
head(top10_bytime_16s)


# genus level

gen_Ct0 = showTopN(my16s.C.t0, 10, "Genus")
gen_Ct1 = showTopN(my16s.C.t1, 10, "Genus")
gen_Ct3 = showTopN(my16s.C.t3, 10, "Genus")
gen_Ct6 = showTopN(my16s.C.t6, 10, "Genus")

gen_Ot0 = showTopN(my16s.O.t0, 10, "Genus")
gen_Ot1 = showTopN(my16s.O.t1, 10, "Genus")
gen_Ot3 = showTopN(my16s.O.t3, 10, "Genus")
gen_Ot6 = showTopN(my16s.O.t6, 10, "Genus")


gen_Gt0 = showTopN(my16s.G.t0, 10, "Genus")
gen_Gt1 = showTopN(my16s.G.t1, 10, "Genus")
gen_Gt3 = showTopN(my16s.G.t3, 10, "Genus")
gen_Gt6 = showTopN(my16s.G.t6, 10, "Genus")


# conventional
genC_top10 = data.frame(Type = rep("Conventional", 
                                           nrow(gen_Ct0)+
                                           nrow(gen_Ct1)+
                                           nrow(gen_Ct3)+
                                           nrow(gen_Ct6)),
                        Time = c(rep("T0", nrow(gen_Ct0)), 
                                 rep("T1", nrow(gen_Ct1)), 
                                 rep("T3", nrow(gen_Ct3)), 
                                 rep("T6", nrow(gen_Ct6))))

top10_genC = cbind(genC_top10, rbind(gen_Ct0, 
                                     gen_Ct1, 
                                     gen_Ct3, 
                                     gen_Ct6))

# organic
genO_top10 = data.frame(Type = rep("Organic", 
                                           nrow(gen_Ot0)+
                                           nrow(gen_Ot1)+
                                           nrow(gen_Ot3)+
                                           nrow(gen_Ot6)),
                        Time = c(rep("T0", nrow(gen_Ot0)), 
                                 rep("T1", nrow(gen_Ot1)), 
                                 rep("T3", nrow(gen_Ot3)), 
                                 rep("T6", nrow(gen_Ot6))))

top10_genO = cbind(genO_top10, rbind(gen_Ot0, 
                                     gen_Ot1, 
                                     gen_Ot3, 
                                     gen_Ot6))


# green
genG_top10 = data.frame(Type = rep("Green", 
                                           nrow(gen_Gt0)+
                                           nrow(gen_Gt1)+
                                           nrow(gen_Gt3)+
                                           nrow(gen_Gt6)),
                        Time = c(rep("T0", nrow(gen_Gt0)), 
                                 rep("T1", nrow(gen_Gt1)), 
                                 rep("T3", nrow(gen_Gt3)), 
                                 rep("T6", nrow(gen_Gt6))))

top10_genG = cbind(genG_top10, rbind(gen_Gt0, 
                                     gen_Gt1, 
                                     gen_Gt3, 
                                     gen_Gt6))


top10_bytime_16s_genus = rbind(top10_genC, top10_genO, top10_genG)

saveRDS(top10_bytime_16s, "./microcosm_soil_analysis_output/taxonomy_composition/top10_bytimie_phylum_16s.rds")
saveRDS(top10_bytime_16s_genus, "./microcosm_soil_analysis_output/taxonomy_composition/top10_bytimie_genus_16s.rds")

head(top10_bytime_16s_genus)

write.csv(top10_bytime_16s, "./microcosm_soil_analysis_output/taxonomy_composition/top10_bytime_phyla_16s.csv", quote = F, row.names = F)
write.csv(top10_bytime_16s_genus, "./microcosm_soil_analysis_output/taxonomy_composition/top10_bytime_genus_16s.csv", quote = F, row.names = F)


head(top10_bytime_16s)
tail(top10_bytime_16s)


types = c(`Conventional`="Conventional",
          `Green`="Mixed",
          `Organic`="Organic")

svg("./microcosm_soil_analysis_output/taxonomy_composition/top10_phyla_16s.svg", width = 12, height = 8)
(graph_top10_phy.16s = ggplot(top10_bytime_16s, aes(x = Time, y = 100*avg_rel_abd, color = Phylum, group = Phylum)) +
        geom_point(size=3) +
        geom_line() + 
        facet_grid( ~ Type, labeller = as_labeller(types)) +
        scale_color_manual(values = colorRampPalette(glasbey())(40)) +
        theme_bw() +
        labs(y="Mean relative abundance (%)", x = "") +
        theme(legend.position = "bottom")
)        
dev.off()

svg("./microcosm_soil_analysis_output/taxonomy_composition/top10_genus_16s.svg", width = 12, height = 8)
(graph_top10_gen.16s = ggplot(top10_bytime_16s_genus, aes(x = Time, y = 100*avg_rel_abd, color = Genus, group = Genus)) +
                geom_point(size=3) +
                geom_line() + 
                facet_grid( ~ Type, labeller = as_labeller(types)) +
                scale_color_manual(values = colorRampPalette(glasbey())(30)) +
                theme_bw() +
                labs(y="Mean relative abundance (%)", x = "Time(week)") + 
                theme(legend.position = "bottom")
)        
dev.off()
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# top 10 fungi
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

showTopN(mydata.its.pruned, n = 10, taxlevel = "Phylum")

# Selecting by avg_rel_abd
# # A tibble: 10 x 2
# Phylum                 avg_rel_abd
# <fct>                        <dbl>
#         1 Ascomycota               0.795    
# 2 Basidiomycota            0.0781   
# 3 Mortierellomycota        0.0772   
# 4 Chytridiomycota          0.0268   
# 5 unidentified             0.0216   
# 6 Glomeromycota            0.000305 
# 7 Basidiobolomycota        0.000238 
# 8 Mucoromycota             0.000173 
# 9 Olpidiomycota            0.0000900
# 10 Calcarisporiellomycota   0.0000763

showTopN(mydata.its.pruned, n = 10, taxlevel = "Genus")

# Genus            avg_rel_abd
# <fct>                  <dbl>
#         1 Pseudogymnoascus      0.0773
# 2 Mortierella           0.0751
# 3 Gibberella            0.0592
# 4 Fusarium              0.0563
# 5 Gibellulopsis         0.0514
# 6 Alternaria            0.0390
# 7 Plectosphaerella      0.0357
# 8 Solicoccozyma         0.0274
# 9 f__Chaetomiaceae      0.0233
# 10 Acremonium            0.0222

myits.C.t0 = subset_samples(mydata.its.pruned, cul_type == "Conventional" & period == "0wk")
myits.C.t1 = subset_samples(mydata.its.pruned, cul_type == "Conventional" & period == "1wk")
myits.C.t3 = subset_samples(mydata.its.pruned, cul_type == "Conventional" & period == "3wk")
myits.C.t6 = subset_samples(mydata.its.pruned, cul_type == "Conventional" & period == "6wk")
myits.O.t0 = subset_samples(mydata.its.pruned, cul_type == "Organic" & period == "0wk")
myits.O.t1 = subset_samples(mydata.its.pruned, cul_type == "Organic" & period == "1wk")
myits.O.t3 = subset_samples(mydata.its.pruned, cul_type == "Organic" & period == "3wk")
myits.O.t6 = subset_samples(mydata.its.pruned, cul_type == "Organic" & period == "6wk")
myits.G.t0 = subset_samples(mydata.its.pruned, cul_type == "Green" & period == "0wk")
myits.G.t1 = subset_samples(mydata.its.pruned, cul_type == "Green" & period == "1wk")
myits.G.t3 = subset_samples(mydata.its.pruned, cul_type == "Green" & period == "3wk")
myits.G.t6 = subset_samples(mydata.its.pruned, cul_type == "Green" & period == "6wk")


phy_Ct0.its = showTopN(myits.C.t0, 10, "Phylum")
phy_Ct1.its = showTopN(myits.C.t1, 10, "Phylum")
phy_Ct3.its = showTopN(myits.C.t3, 10, "Phylum")
phy_Ct6.its = showTopN(myits.C.t6, 10, "Phylum")
phy_Ot0.its = showTopN(myits.O.t0, 10, "Phylum")
phy_Ot1.its = showTopN(myits.O.t1, 10, "Phylum")
phy_Ot3.its = showTopN(myits.O.t3, 10, "Phylum")
phy_Ot6.its = showTopN(myits.O.t6, 10, "Phylum")
phy_Gt0.its = showTopN(myits.G.t0, 10, "Phylum")
phy_Gt1.its = showTopN(myits.G.t1, 10, "Phylum")
phy_Gt3.its = showTopN(myits.G.t3, 10, "Phylum")
phy_Gt6.its = showTopN(myits.G.t6, 10, "Phylum")


phyC_top10.its = data.frame(Type = rep("Conventional", 
                                           nrow(phy_Ct0.its)+
                                           nrow(phy_Ct1.its)+
                                           nrow(phy_Ct3.its)+
                                           nrow(phy_Ct6.its)),
                        Time = c(rep("T0", nrow(phy_Ct0.its)), 
                                 rep("T1", nrow(phy_Ct1.its)), 
                                 rep("T3", nrow(phy_Ct3.its)), 
                                 rep("T6", nrow(phy_Ct6.its))))

top10_phyC.its = cbind(phyC_top10.its, rbind(phy_Ct0.its, 
                                         phy_Ct1.its, 
                                         phy_Ct3.its, 
                                         phy_Ct6.its))

# organic
phyO_top10.its = data.frame(Type = rep("Organic", 
                                           nrow(phy_Ot0.its)+
                                           nrow(phy_Ot1.its)+
                                           nrow(phy_Ot3.its)+
                                           nrow(phy_Ot6.its)),
                        Time = c(rep("T0", nrow(phy_Ot0.its)), 
                                 rep("T1", nrow(phy_Ot1.its)), 
                                 rep("T3", nrow(phy_Ot3.its)), 
                                 rep("T6", nrow(phy_Ot6.its))))

top10_phyO.its = cbind(phyO_top10.its, rbind(phy_Ot0.its, 
                                             phy_Ot1.its, 
                                             phy_Ot3.its, 
                                             phy_Ot6.its))


# green
phyG_top10.its = data.frame(Type = rep("Green", 
                                           nrow(phy_Gt0.its)+
                                           nrow(phy_Gt1.its)+
                                           nrow(phy_Gt3.its)+
                                           nrow(phy_Gt6.its)),
                        Time = c(rep("T0", nrow(phy_Gt0.its)), 
                                 rep("T1", nrow(phy_Gt1.its)), 
                                 rep("T3", nrow(phy_Gt3.its)), 
                                 rep("T6", nrow(phy_Gt6.its))))

top10_phyG.its = cbind(phyG_top10.its, rbind(phy_Gt0.its, 
                                             phy_Gt1.its, 
                                             phy_Gt3.its, 
                                             phy_Gt6.its))


# top10 16s in Phylum along time

top10_bytime_its = rbind(top10_phyC.its, top10_phyO.its, top10_phyG.its)
head(top10_bytime_its)


# genus level

gen_Ct0.its = showTopN(myits.C.t0, 10, "Genus")
gen_Ct1.its = showTopN(myits.C.t1, 10, "Genus")
gen_Ct3.its = showTopN(myits.C.t3, 10, "Genus")
gen_Ct6.its = showTopN(myits.C.t6, 10, "Genus")
gen_Ot0.its = showTopN(myits.O.t0, 10, "Genus")
gen_Ot1.its = showTopN(myits.O.t1, 10, "Genus")
gen_Ot3.its = showTopN(myits.O.t3, 10, "Genus")
gen_Ot6.its = showTopN(myits.O.t6, 10, "Genus")
gen_Gt0.its = showTopN(myits.G.t0, 10, "Genus")
gen_Gt1.its = showTopN(myits.G.t1, 10, "Genus")
gen_Gt3.its = showTopN(myits.G.t3, 10, "Genus")
gen_Gt6.its = showTopN(myits.G.t6, 10, "Genus")


# conventional
genC_top10.its = data.frame(Type = rep("Conventional", 
                                           nrow(gen_Ct0.its)+
                                           nrow(gen_Ct1.its)+
                                           nrow(gen_Ct3.its)+
                                           nrow(gen_Ct6.its)),
                        Time = c(rep("T0", nrow(gen_Ct0.its)), 
                                 rep("T1", nrow(gen_Ct1.its)), 
                                 rep("T3", nrow(gen_Ct3.its)), 
                                 rep("T6", nrow(gen_Ct6.its))))

top10_genC.its = cbind(genC_top10.its, rbind(gen_Ct0.its, 
                                             gen_Ct1.its, 
                                             gen_Ct3.its, 
                                             gen_Ct6.its))

# organic
genO_top10.its = data.frame(Type = rep("Organic", 
                                           nrow(gen_Ot0.its)+
                                           nrow(gen_Ot1.its)+
                                           nrow(gen_Ot3.its)+
                                           nrow(gen_Ot6.its)),
                        Time = c(rep("T0", nrow(gen_Ot0.its)), 
                                 rep("T1", nrow(gen_Ot1.its)), 
                                 rep("T3", nrow(gen_Ot3.its)), 
                                 rep("T6", nrow(gen_Ot6.its))))

top10_genO.its = cbind(genO_top10.its, rbind(gen_Ot0.its, 
                                             gen_Ot1.its, 
                                             gen_Ot3.its, 
                                             gen_Ot6.its))


# green
genG_top10.its = data.frame(Type = rep("Green", 
                                           nrow(gen_Gt0.its)+
                                           nrow(gen_Gt1.its)+
                                           nrow(gen_Gt3.its)+
                                           nrow(gen_Gt6.its)),
                        Time = c(rep("T0", nrow(gen_Gt0.its)), 
                                 rep("T1", nrow(gen_Gt1.its)), 
                                 rep("T3", nrow(gen_Gt3.its)), 
                                 rep("T6", nrow(gen_Gt6.its))))

top10_genG.its = cbind(genG_top10.its, rbind(gen_Gt0.its, 
                                             gen_Gt1.its, 
                                             gen_Gt3.its, 
                                             gen_Gt6.its))


top10_bytime_its_genus = rbind(top10_genC.its, top10_genO.its, top10_genG.its)
head(top10_bytime_its_genus)


saveRDS(top10_bytime_16s, "./microcosm_soil_analysis_output/taxonomy_composition/top10_bytime_phylum_16s.rds")
saveRDS(top10_bytime_its, "./microcosm_soil_analysis_output/taxonomy_composition/top10_bytime_phylum_its.rds")
saveRDS(top10_bytime_16s_genus, "./microcosm_soil_analysis_output/taxonomy_composition/top10_bytime_genus_16s.rds")
saveRDS(top10_bytime_its_genus, "./microcosm_soil_analysis_output/taxonomy_composition/top10_bytime_genus_its.rds")



write.csv(top10_bytime_its, "./microcosm_soil_analysis_output/taxonomy_composition/top10_bytime_phyla_its.csv", quote = F, row.names = F)
write.csv(top10_bytime_its_genus, "./microcosm_soil_analysis_output/taxonomy_composition/top10_bytime_genus_its.csv", quote = F, row.names = F)


head(top10_bytime_its)
tail(top10_bytime_its)

top10_bytime_its$Phylum = as.character(top10_bytime_its$Phylum)


svg("./microcosm_soil_analysis_output/taxonomy_composition/top10_phyla_its.svg", width = 12, height = 8)
(graph_top10_phy.its = top10_bytime_its %>%
                ggplot(aes(x = Time, y = 100*avg_rel_abd, color = Phylum, group = Phylum)) +
                geom_point(size=3) +
                geom_line() + 
                facet_grid( ~ Type, labeller = as_labeller(types)) +
                scale_color_manual(values = colorRampPalette(glasbey())(30)) +
                theme_bw() +
                labs(y="Mean relative abundance (%)", x = "Time(week)") +
                theme(legend.position = "bottom")
)        
dev.off()

svg("./microcosm_soil_analysis_output/taxonomy_composition/top10_genus_its.svg", width = 12, height = 8)
(graph_top10_gen.its = 
                top10_bytime_its_genus %>%
                ggplot(aes(x = Time, y = 100*avg_rel_abd, color = Genus, group = Genus)) +
                geom_point(size=3) +
                geom_line() + 
                facet_grid( ~ Type, labeller = as_labeller(types)) +
                scale_color_manual(values = colorRampPalette(glasbey())(40)) +
                theme_bw() +
                labs(y="Mean relative abundance (%)", x = "Time(week)") + 
                theme(legend.position = "bottom")
)     
dev.off()

svg("./microcosm_soil_analysis_output/taxonomy_composition/top10.svg", width = 12, height = 8)
ggarrange(graph_top10_phy.16s, graph_top10_phy.its, labels=c("C", "D"), ncol = 1, legend = "right", align = "hv")
dev.off()


######################################################################################

# DESeq2
rm(list=ls())

source("./tools/plotDeseq2.R")
dir.create("./microcosm_soil_analysis_output/DESeq2")
dir.create("./microcosm_soil_analysis_output/DESeq2/16s")
dir.create("./microcosm_soil_analysis_output/DESeq2/its")

mydata.16s.pruned = readRDS("./microcosm_soil_analysis_output/data/mydata_16s_pruned.rds")
mydata.its.pruned = readRDS("./microcosm_soil_analysis_output/data/mydata_its_pruned.rds")

env.16s = as.data.frame(sample_data(mydata.16s.pruned))
env.its = as.data.frame(sample_data(mydata.its.pruned))

# we will use contrast_group for contrast comparison
env.16s$contrast_group = paste(env.16s$cul_type, env.16s$period, sep="_")
env.its$contrast_group = paste(env.its$cul_type, env.its$period, sep="_")

sample_data(mydata.16s.pruned) = env.16s
sample_data(mydata.its.pruned) = env.its

m16s.dds = phyloseq_to_deseq2(mydata.16s.pruned, ~ contrast_group)
mits.dds = phyloseq_to_deseq2(mydata.its.pruned, ~ contrast_group)


# geomean
gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

######################################################################################




#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# Bacterial communities
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

geoMean.16s = apply(counts(m16s.dds), 1, gm_mean)
m16s.dds.mean = estimateSizeFactors(m16s.dds, geoMean = geoMean.16s)
m16s.dds.mean = DESeq(m16s.dds.mean, fitType = "local")

png("./microcosm_soil_analysis_output/DESeq2/disper_geomean_16s.png")
plotDispEsts(m16s.dds.mean)
dev.off()


# saveRDS(m16s.dds.mean, "./microcosm_soil_analysis_output/DESeq2/m16s_dds_Cultype_Time.rds")

(g16s_C_1wk_vs_0wk = plotDeseq2(m16s.dds.mean, contrasts = c("contrast_group", "Conventional_1wk", "Conventional_0wk"), mydata.16s.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = "./microcosm_soil_analysis_output/DESeq2/16s/"))



(g16s_C_3wk_vs_0wk = plotDeseq2(m16s.dds.mean, contrasts = c("contrast_group", "Conventional_3wk", "Conventional_0wk"), mydata.16s.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = "./microcosm_soil_analysis_output/DESeq2/16s/"))




(g16s_C_6wk_vs_0wk = plotDeseq2(m16s.dds.mean, contrasts = c("contrast_group", "Conventional_6wk", "Conventional_0wk"), mydata.16s.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = "./microcosm_soil_analysis_output/DESeq2/16s/"))



g16s_O_1wk_vs_0wk = plotDeseq2(m16s.dds.mean, contrasts = c("contrast_group", "Organic_1wk", "Organic_0wk"), mydata.16s.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = "./microcosm_soil_analysis_output/DESeq2/16s/")
g16s_O_3wk_vs_0wk = plotDeseq2(m16s.dds.mean, contrasts = c("contrast_group", "Organic_3wk", "Organic_0wk"), mydata.16s.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = "./microcosm_soil_analysis_output/DESeq2/16s/")
g16s_O_6wk_vs_0wk = plotDeseq2(m16s.dds.mean, contrasts = c("contrast_group", "Organic_6wk", "Organic_0wk"), mydata.16s.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = "./microcosm_soil_analysis_output/DESeq2/16s/")

g16s_G_1wk_vs_0wk = plotDeseq2(m16s.dds.mean, contrasts = c("contrast_group", "Green_1wk", "Green_0wk"), mydata.16s.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = "./microcosm_soil_analysis_output/DESeq2/16s/")
g16s_G_3wk_vs_0wk = plotDeseq2(m16s.dds.mean, contrasts = c("contrast_group", "Green_3wk", "Green_0wk"), mydata.16s.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = "./microcosm_soil_analysis_output/DESeq2/16s/")
g16s_G_6wk_vs_0wk = plotDeseq2(m16s.dds.mean, contrasts = c("contrast_group", "Green_6wk", "Green_0wk"), mydata.16s.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = "./microcosm_soil_analysis_output/DESeq2/16s/")

# save the objects because it is time consuming to produce them
g16s_comparisons = list(
        g16s_C_1wk_vs_0wk = g16s_C_1wk_vs_0wk,
        g16s_C_3wk_vs_0wk = g16s_C_3wk_vs_0wk,
        g16s_C_6wk_vs_0wk = g16s_C_6wk_vs_0wk,
        g16s_O_1wk_vs_0wk = g16s_O_1wk_vs_0wk,
        g16s_O_3wk_vs_0wk = g16s_O_3wk_vs_0wk,
        g16s_O_6wk_vs_0wk = g16s_O_6wk_vs_0wk,
        g16s_G_1wk_vs_0wk = g16s_G_1wk_vs_0wk,
        g16s_G_3wk_vs_0wk = g16s_G_3wk_vs_0wk,
        g16s_G_6wk_vs_0wk = g16s_G_6wk_vs_0wk
        
)

saveRDS(g16s_comparisons, "./microcosm_soil_analysis_output/DESeq2/g16s_comparisons.rds")

svg("./microcosm_soil_analysis_output/DESeq2/16s/Conventional_1wk_vs_Conventional_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
g16s_C_1wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/16s/Conventional_3wk_vs_Conventional_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
g16s_C_3wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/16s/Conventional_6wk_vs_Conventional_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
g16s_C_6wk_vs_0wk$graph
dev.off()


svg("./microcosm_soil_analysis_output/DESeq2/16s/Org_1wk_vs_Org_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
g16s_O_1wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/16s/Org_3wk_vs_Org_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
g16s_O_3wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/16s/Org_6wk_vs_Org_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
g16s_O_6wk_vs_0wk$graph
dev.off()


svg("./microcosm_soil_analysis_output/DESeq2/16s/Grn_1wk_vs_Grn_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
g16s_G_1wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/16s/Grn_3wk_vs_Grn_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
g16s_G_3wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/16s/Grn_6wk_vs_Grn_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
g16s_G_6wk_vs_0wk$graph
dev.off()




######################################################
#                                                    #
#                    ITS                             #
######################################################



geoMean.its = apply(counts(mits.dds), 1, gm_mean)
mits.dds.mean = estimateSizeFactors(mits.dds, geoMean = geoMean.its)
mits.dds.mean = DESeq(mits.dds.mean, fitType = "local")
saveRDS(mits.dds.mean, "./microcosm_soil_analysis_output/DESeq2/mits_dds_Cultype_Time.rds")



tiff("./microcosm_soil_analysis_output/DESeq2/disper_geomean_its.tiff")
plotDispEsts(mits.dds.mean)
dev.off()

outdds_its = file.path("./microcosm_soil_analysis_output/DESeq2/its/")

gits_C_1wk_vs_0wk = plotDeseq2(mits.dds.mean, contrasts = c("contrast_group", "Conventional_1wk", "Conventional_0wk"), mydata.its.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = outdds_its)

gits_C_3wk_vs_0wk = plotDeseq2(mits.dds.mean, contrasts = c("contrast_group", "Conventional_3wk", "Conventional_0wk"), mydata.its.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = outdds_its)

gits_C_6wk_vs_0wk = plotDeseq2(mits.dds.mean, contrasts = c("contrast_group", "Conventional_6wk", "Conventional_0wk"), mydata.its.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = outdds_its)



gits_O_1wk_vs_0wk = plotDeseq2(mits.dds.mean, contrasts = c("contrast_group", "Organic_1wk", "Organic_0wk"), mydata.its.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = outdds_its)

gits_O_3wk_vs_0wk = plotDeseq2(mits.dds.mean, contrasts = c("contrast_group", "Organic_3wk", "Organic_0wk"), mydata.its.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = outdds_its)

gits_O_6wk_vs_0wk = plotDeseq2(mits.dds.mean, contrasts = c("contrast_group", "Organic_6wk", "Organic_0wk"), mydata.its.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = outdds_its)


gits_G_1wk_vs_0wk = plotDeseq2(mits.dds.mean, contrasts = c("contrast_group", "Green_1wk", "Green_0wk"), mydata.its.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = outdds_its)

gits_G_3wk_vs_0wk = plotDeseq2(mits.dds.mean, contrasts = c("contrast_group", "Green_3wk", "Green_0wk"), mydata.its.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = outdds_its)

gits_G_6wk_vs_0wk = plotDeseq2(mits.dds.mean, contrasts = c("contrast_group", "Green_6wk", "Green_0wk"), mydata.its.pruned, tax_level = "Genus", alpha = 0.01, plottype = "scatter", log2FC = 2, legendPos = "bottom", outputpath = outdds_its)



# save the objects because it is time consuming to produce them
gits_comparisons = list(
        gits_C_1wk_vs_0wk = gits_C_1wk_vs_0wk,
        gits_C_3wk_vs_0wk = gits_C_3wk_vs_0wk,
        gits_C_6wk_vs_0wk = gits_C_6wk_vs_0wk,
        gits_O_1wk_vs_0wk = gits_O_1wk_vs_0wk,
        gits_O_3wk_vs_0wk = gits_O_3wk_vs_0wk,
        gits_O_6wk_vs_0wk = gits_O_6wk_vs_0wk,
        gits_G_1wk_vs_0wk = gits_G_1wk_vs_0wk,
        gits_G_3wk_vs_0wk = gits_G_3wk_vs_0wk,
        gits_G_6wk_vs_0wk = gits_G_6wk_vs_0wk
        
)

saveRDS(gits_comparisons, "./microcosm_soil_analysis_output/DESeq2/gits_comparisons.rds")





svg("./microcosm_soil_analysis_output/DESeq2/its/Conventional_1wk_vs_Conventional_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
gits_C_1wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/its/Conventional_3wk_vs_Conventional_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
gits_C_3wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/its/Conventional_6wk_vs_Conventional_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
gits_C_6wk_vs_0wk$graph
dev.off()


svg("./microcosm_soil_analysis_output/DESeq2/its/Org_1wk_vs_Org_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
gits_O_1wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/its/Org_3wk_vs_Org_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
gits_O_3wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/its/Org_6wk_vs_Org_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
gits_O_6wk_vs_0wk$graph
dev.off()


svg("./microcosm_soil_analysis_output/DESeq2/its/Grn_1wk_vs_Grn_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
gits_G_1wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/its/Grn_3wk_vs_Grn_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
gits_G_3wk_vs_0wk$graph
dev.off()
svg("./microcosm_soil_analysis_output/DESeq2/its/Grn_6wk_vs_Grn_0wk_all_significantOTU_level0.01.svg", width = 12, height = 10)
gits_G_6wk_vs_0wk$graph
dev.off()


######################################
# combine all
######################################

g16s_df = do.call(rbind, list( g16s_C_1wk_vs_0wk$df, 
                               g16s_C_3wk_vs_0wk$df,
                               g16s_C_6wk_vs_0wk$df,
                               g16s_O_1wk_vs_0wk$df,
                               g16s_O_3wk_vs_0wk$df,
                               g16s_O_6wk_vs_0wk$df,
                               g16s_G_1wk_vs_0wk$df,
                               g16s_G_3wk_vs_0wk$df,
                               g16s_G_6wk_vs_0wk$df))

g16s_df$community = "Bacteria"

# saveRDS(g16s_df, "./microcosm_soil_analysis_output/DESeq2/g16s_df.rds")


gits_df = do.call(rbind, list( gits_C_1wk_vs_0wk$df, 
                               gits_C_3wk_vs_0wk$df,
                               gits_C_6wk_vs_0wk$df,
                               gits_O_1wk_vs_0wk$df,
                               gits_O_3wk_vs_0wk$df,
                               gits_O_6wk_vs_0wk$df,
                               gits_G_1wk_vs_0wk$df,
                               gits_G_3wk_vs_0wk$df,
                               gits_G_6wk_vs_0wk$df))

gits_df$community = "Fungi"

# saveRDS(gits_df, "./microcosm_soil_analysis_output/DESeq2/gits_df.rds")

x = tapply(g16s_df$log2FoldChange, g16s_df$Phylum, function(x)max(x))
x = sort(x, T)

x

g16s_df$Phylum = factor(as.character(g16s_df$Phylum), levels=names(x))

x = tapply(g16s_df$log2FoldChange, g16s_df$Genus, function(x)max(x))
x = sort(x, T)

g16s_df$Genus = factor(as.character(g16s_df$Genus), levels=names(x))




(graph16s = ggplot(g16s_df, aes(x = Genus, y = log2FoldChange, color = Phylum)) + 
        facet_grid( cultypes ~ contrasts, scales = "free") +
        geom_jitter(size=3, width = 0.0, height = 0.2, alpha=0.4) +
        geom_hline(aes(yintercept = 0), color = "grey") + 
        scale_color_manual(values = glasbey()) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5, size=10), 
              legend.position = "bottom"))


x = tapply(gits_df$log2FoldChange, gits_df$Phylum, function(x)max(x))
x = sort(x, T)



gits_df$Phylum = factor(as.character(gits_df$Phylum), levels=names(x))

x = tapply(gits_df$log2FoldChange, gits_df$Genus, function(x)max(x))
x = sort(x, T)

gits_df$Genus = factor(as.character(gits_df$Genus), levels=names(x))


(graphITS = ggplot(gits_df, aes(x = Genus, y = log2FoldChange, color = Phylum)) + 
                facet_grid( cultypes ~ contrasts, scales = "free") +
                geom_jitter(size=3, width = 0.0, height = 0.2, alpha=0.4) +
                geom_hline(aes(yintercept = 0), color = "grey") + 
                theme_bw() +
                scale_color_manual(values=brewer.dark2(5)) +
                theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5, size=10), 
                      legend.position = "bottom"))

graphITS

outdds = "./microcosm_soil_analysis_output/DESeq2/"
svg(file.path(outdds, "scatter_16s_deseq2.svg"), width = 18, height = 12, pointsize = 15)
graph16s
dev.off()


svg(file.path(outdds, "scatter_its_deseq2.svg"), width = 22, height = 14)
graphITS
dev.off()



########################################################
# make tree plot for DESeq2 ASVs detected significant (0.01)


rm(list = ls())
source("./tools/PlotTree.R")
mydata.16s.pruned = readRDS("./microcosm_soil_analysis_output/data/mydata_16s_pruned.rds")
mydata.its.pruned = readRDS("./microcosm_soil_analysis_output/data/mydata_its_pruned.rds")
#########################################################
# analyzing tree 16s

# loading signficant ASVs

conv_all_sig.16s = c()
org_all_sig.16s = c()
grn_all_sig.16s = c()
for(i in list.files("./microcosm_soil_analysis_output/DESeq2/16s/")){
        if(startsWith(i, "Conventional") & endsWith(i, "level0.01.csv")){
                conv_all_sig.16s = c(conv_all_sig.16s, i)
        } else if(startsWith(i, "Organic") & endsWith(i, "level0.01.csv")){
                org_all_sig.16s = c(org_all_sig.16s, i)
        } else if(startsWith(i, "Green") & endsWith(i, "level0.01.csv")){
                grn_all_sig.16s = c(grn_all_sig.16s, i)
        }
}



conv_all_sig.16s
conv_df.16s = lapply(conv_all_sig.16s, addTime, filepath = "./microcosm_soil_analysis_output/DESeq2/16s/")

# add time indicator to each dataframe
conv_df2.16s = do.call(rbind, conv_df.16s)
conv_df2.16s[1:5, 1:5]

org_all_sig.16s
org_df.16s = lapply(org_all_sig.16s, addTime, filepath = "./microcosm_soil_analysis_output/DESeq2/16s/")
org_df2.16s = do.call(rbind, org_df.16s)

grn_df.16s = lapply(grn_all_sig.16s, addTime, filepath = "./microcosm_soil_analysis_output/DESeq2/16s/")
grn_df2.16s = do.call(rbind, grn_df.16s)


ASVs_C.16s = as.character(conv_df2.16s$X)
ASVs_O.16s = as.character(org_df2.16s$X)
ASVs_G.16s = as.character(grn_df2.16s$X)


C16s = prune_taxa(ASVs_C.16s, mydata.16s.pruned)
O16s = prune_taxa(ASVs_O.16s, mydata.16s.pruned)
G16s = prune_taxa(ASVs_G.16s, mydata.16s.pruned)


# plot tree

tax_table(C16s)[, "Genus"] = gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "A-N-P-Rhizobium", tax_table(C16s)[, "Genus"])

tax_table(G16s)[tax_table(G16s)[, "Genus"] == "uncultured alpha proteobacterium"]

tax_table(O16s)[, "Genus"] = gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "A-N-P-Rhizobium", tax_table(O16s)[, "Genus"])

tax_table(G16s)[, "Genus"] = gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "A-N-P-Rhizobium", tax_table(G16s)[, "Genus"])

UseColors = colorRampPalette(brewer.pal(9, "Set1"))

Phylum16s = c(tax_table(C16s)[, "Phylum"], tax_table(O16s)[, "Phylum"], tax_table(G16s)[, "Phylum"])
assignColors = data.frame(Phylum = unique(Phylum16s), Colors = UseColors(length(unique(Phylum16s))))
assignCol = as.character(assignColors[, "Colors"])
assignCol
names(assignCol) = unique(Phylum16s)

source("./tools/PlotTree.R")
# set same scale also modified the code to add function that transforms matrix
conv_matrix.16S = transformMatrix(conv_df2.16s, colselected = c("1wk", "3wk", "6wk"), otucol = "X")
org_matrix.16S = transformMatrix(org_df2.16s, colselected = c("1wk", "3wk", "6wk"), otucol = "X")
grn_matrix.16S = transformMatrix(grn_df2.16s, colselected = c("1wk", "3wk", "6wk"), otucol = "X")

rngMtrix = range(conv_matrix.16S, org_matrix.16S, grn_matrix.16S)
rngMtrix # -8.121, 11.04



# save graphs

svg("./microcosm_soil_analysis_output/DESeq2/treeplot_16s_C_scaled.svg", width = 12, height = 14)
(C16s_graph = plotggtree(C16s, conv_matrix.16S, scaleColorRange = T, rngMtrix = rngMtrix, breaks_by = 3,colFUN = assignCol, changeColnames = c("1", "3", "6"), addColnames = T, hmap_offset = 3, hmap_width = 0.5, colnames_offset_y = -0.5) + 
                guides(color = guide_legend(order=0), fill = guide_legend(order = 1)))
dev.off()


svg("./microcosm_soil_analysis_output/DESeq2/treeplot_16s_O_scaled.svg", width = 12, height = 14)
(O16s_graph = plotggtree(O16s, org_matrix.16S, colFUN = assignCol,  changeColnames = c("1", "3", "6"), scaleColorRange = T, rngMtrix = rngMtrix, breaks_by = 3, addColnames = T, hmap_offset = 3, hmap_width = 0.5, colnames_offset_y = -0.5) + 
                guides(color = guide_legend(order=0), fill = guide_legend(order = 1)))
dev.off()


svg("./microcosm_soil_analysis_output/DESeq2/treeplot_16s_G_scaled.svg", width = 12, height = 14)
(G16s_graph = out = plotggtree(G16s, grn_matrix.16S, colFUN = assignCol,changeColnames = c("1", "3", "6"), scaleColorRange = T, rngMtrix = rngMtrix, breaks_by = 3, hmap_offset = 3, hmap_width = 0.5, colnames_offset_y = -0.5)+ 
                guides(color = guide_legend(order=0), fill = guide_legend(order = 1)))
dev.off()



#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
#Fungi tree
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

# analyzing tree ITS

# loading signficant ASVs

conv_all_sig.its = c()
org_all_sig.its = c()
grn_all_sig.its = c()
for(i in list.files("./microcosm_soil_analysis_output/DESeq2/its/")){
        if(startsWith(i, "Conventional") & endsWith(i, "level0.01.csv")){
                conv_all_sig.its = c(conv_all_sig.its, i)
        } else if(startsWith(i, "Organic") & endsWith(i, "level0.01.csv")){
                org_all_sig.its = c(org_all_sig.its, i)
        } else if(startsWith(i, "Green") & endsWith(i, "level0.01.csv")){
                grn_all_sig.its = c(grn_all_sig.its, i)
        }
}


# put data frame together

conv_all_sig.its
# add time indicator to each dataframe
conv_df.its = lapply(conv_all_sig.its, addTime, filepath = "./microcosm_soil_analysis_output/DESeq2/its/")
# combine
conv_df2.its = do.call(rbind, conv_df.its)

org_all_sig.its
org_df.its = lapply(org_all_sig.its, addTime, filepath = "./microcosm_soil_analysis_output/DESeq2/its/")
org_df2.its = do.call(rbind, org_df.its)

grn_df.its = lapply(grn_all_sig.its, addTime, filepath = "./microcosm_soil_analysis_output/DESeq2/its/")
grn_df2.its = do.call(rbind, grn_df.its)


ASVs_C.its = as.character(conv_df2.its$X)
ASVs_O.its = as.character(org_df2.its$X)
ASVs_G.its = as.character(grn_df2.its$X)


Cits = prune_taxa(ASVs_C.its, mydata.its.pruned)
Oits = prune_taxa(ASVs_O.its, mydata.its.pruned)
Gits = prune_taxa(ASVs_G.its, mydata.its.pruned)


source("./tools/PlotTree.R")
# set same scale also modified the code to add function that transforms matrix
conv_matrix.ITS = transformMatrix(conv_df2.its, colselected = c("1wk", "3wk", "6wk"), otucol = "X")
org_matrix.ITS = transformMatrix(org_df2.its, colselected = c("1wk", "3wk", "6wk"), otucol = "X")
grn_matrix.ITS = transformMatrix(grn_df2.its, colselected = c("1wk", "3wk", "6wk"), otucol = "X")

rngMtrix = range(conv_matrix.ITS, org_matrix.ITS, grn_matrix.ITS)
rngMtrix #-8.89311 15.16232

# plot tree


# assign save color for phylum for each plot
UseColors = colorRampPalette(brewer.pal(9, "Set1"))

PhylumITS = c(tax_table(Cits)[, "Phylum"], tax_table(Oits)[, "Phylum"], tax_table(Gits)[, "Phylum"])
assignColors = data.frame(Phylum = unique(PhylumITS), Colors = UseColors(length(unique(PhylumITS))))
assignCol = as.character(assignColors[, "Colors"])
assignCol
names(assignCol) = unique(PhylumITS)


svg("./microcosm_soil_analysis_output/DESeq2/treeplot_its_C_scaled.svg", width = 12, height = 16)
(Cits_graph = plotggtree(Cits, conv_matrix.ITS, colFUN = assignCol, changeColnames = c("1", "3", "6"), scaleColorRange = T, rngMtrix = rngMtrix, breaks_by = 3, addColnames = T, hmap_offset = 3, hmap_width = 0.5, colnames_offset_y = -0.5)+
                guides(color = guide_legend(order=0), fill = guide_legend(order = 1)))
dev.off()


svg("./microcosm_soil_analysis_output/DESeq2/treeplot_its_O_scaled.svg", width = 12, height = 16)
(Oits_graph = plotggtree(Oits, org_matrix.ITS, colFUN = assignCol, scaleColorRange = T, rngMtrix = rngMtrix, breaks_by = 3, changeColnames = c("1", "3", "6"), addColnames = T, hmap_offset = 3, hmap_width = 0.5, colnames_offset_y = -0.5) + 
                guides(color = guide_legend(order=0), fill = guide_legend(order = 1)))
dev.off()


svg("./microcosm_soil_analysis_output/DESeq2/treeplot_its_G_scaled.svg", width = 12, height = 16)
(Gits_graph = out = plotggtree(Gits, grn_matrix.ITS, colFUN = assignCol, scaleColorRange = T, rngMtrix = rngMtrix, breaks_by = 3,  changeColnames = c("1", "3", "6"), hmap_offset = 3, hmap_width = 0.5, colnames_offset_y = -0.5)+ guides(color = guide_legend(order=0), fill = guide_legend(order = 1)))
dev.off()



# ####################################################################################################
# 
# # network
# rm(list=ls())
# source("./tools/phyloTools.R")
# source("./tools/networkTools.R")
# 
# 
# bst16s = readRDS("./microcosm_soil_analysis_output/data/mydata_16s_pruned.rds")
# bstITS = readRDS("./microcosm_soil_analysis_output/data/mydata_its_pruned.rds")
# 
# 
# dir.create("./microcosm_soil_analysis_output/network/16s", recursive = T)
# dir.create("./microcosm_soil_analysis_output/network/its", recursive = T)
# 
# 
# out16s = "./microcosm_soil_analysis_output/network/16s"
# outits = "./microcosm_soil_analysis_output/network/its"
# 
# 
# ####################################################################################################
# 
# 
# 
# ############ d10
# 
# bst16s.C.d10t0 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "10cm" & period == "0wk")
# bst16s.C.d10t1 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "10cm" & period == "1wk")
# bst16s.C.d10t3 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "10cm" & period == "3wk")
# bst16s.C.d10t6 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "10cm" & period == "6wk")
# 
# bst16s.O.d10t0 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "10cm" & period == "0wk")
# bst16s.O.d10t1 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "10cm" & period == "1wk")
# bst16s.O.d10t3 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "10cm" & period == "3wk")
# bst16s.O.d10t6 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "10cm" & period == "6wk")
# 
# bst16s.G.d10t0 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "10cm" & period == "0wk")
# bst16s.G.d10t1 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "10cm" & period == "1wk")
# bst16s.G.d10t3 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "10cm" & period == "3wk")
# bst16s.G.d10t6 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "10cm" & period == "6wk")
# 
# ###############d20
# 
# bst16s.C.d20t0 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "20cm" & period == "0wk")
# bst16s.C.d20t1 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "20cm" & period == "1wk")
# bst16s.C.d20t3 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "20cm" & period == "3wk")
# bst16s.C.d20t6 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "20cm" & period == "6wk")
# 
# 
# bst16s.O.d20t0 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "20cm" & period == "0wk")
# bst16s.O.d20t1 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "20cm" & period == "1wk")
# bst16s.O.d20t3 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "20cm" & period == "3wk")
# bst16s.O.d20t6 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "20cm" & period == "6wk")
# 
# bst16s.G.d20t0 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "20cm" & period == "0wk")
# bst16s.G.d20t1 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "20cm" & period == "1wk")
# bst16s.G.d20t3 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "20cm" & period == "3wk")
# bst16s.G.d20t6 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "20cm" & period == "6wk")
# 
# ############### d30
# 
# bst16s.C.d30t0 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "30cm" & period == "0wk")
# bst16s.C.d30t1 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "30cm" & period == "1wk")
# bst16s.C.d30t3 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "30cm" & period == "3wk")
# bst16s.C.d30t6 = subset_samples(bst16s, cul_type=="Conventional" & soil_depth == "30cm" & period == "6wk")
# 
# 
# bst16s.O.d30t0 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "30cm" & period == "0wk")
# bst16s.O.d30t1 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "30cm" & period == "1wk")
# bst16s.O.d30t3 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "30cm" & period == "3wk")
# bst16s.O.d30t6 = subset_samples(bst16s, cul_type=="Organic" & soil_depth == "30cm" & period == "6wk")
# 
# bst16s.G.d30t0 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "30cm" & period == "0wk")
# bst16s.G.d30t1 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "30cm" & period == "1wk")
# bst16s.G.d30t3 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "30cm" & period == "3wk")
# bst16s.G.d30t6 = subset_samples(bst16s, cul_type=="Green" & soil_depth == "30cm" & period == "6wk")
# 
# 
# conv16s = list(d10 = list(bst16s.C.d10t0, bst16s.C.d10t1, bst16s.C.d10t3, bst16s.C.d10t6),
#                d20 = list(bst16s.C.d20t0, bst16s.C.d20t1, bst16s.C.d20t3, bst16s.C.d20t6),
#                d30 = list(bst16s.C.d30t0, bst16s.C.d30t1, bst16s.C.d30t3, bst16s.C.d30t6))
# 
# org16s = list(d10 = list(bst16s.O.d10t0, bst16s.O.d10t1, bst16s.O.d10t3, bst16s.O.d10t6),
#               d20 = list(bst16s.O.d20t0, bst16s.O.d20t1, bst16s.O.d20t3, bst16s.O.d20t6),
#               d30 = list(bst16s.O.d30t0, bst16s.O.d30t1, bst16s.O.d30t3, bst16s.O.d30t6))
# 
# grn16s = list(d10 = list(bst16s.G.d10t0, bst16s.G.d10t1, bst16s.G.d10t3, bst16s.G.d10t6),
#               d20 = list(bst16s.G.d20t0, bst16s.G.d20t1, bst16s.G.d20t3, bst16s.G.d20t6),
#               d30 = list(bst16s.G.d30t0, bst16s.G.d30t1, bst16s.G.d30t3, bst16s.G.d30t6))
# 
# 
# cultypes16s =  list(Conventional = conv16s, Organic=org16s, Green=grn16s)
# 
# 
# periods = c("t0", "t1", "t3", "t6")
# 
# for(i in c("Conventional", "Organic", "Green")){
#         for(j in c("d10", "d20", "d30")){
#                 
#                 
#                 each_depth = cultypes16s[[i]][[j]]
#                 
#                 dir.create(file.path(out16s, paste0("datasource/", i, "/", j)), recursive = T)
#                 
#                 out = file.path(out16s, paste0("datasource/", i, "/", j))
#                 
#                 for(k in 1:length(each_depth)){
#                         
#                         otu_out = paste0("otu_", substr(i,1,1), "_", paste0(j, "t", periods[k]))
#                         tax_out = paste0("tax_", substr(i,1,1), "_", paste0(j, "t", periods[k]))
#                         meta_out = paste0("meta_", substr(i,1,1), "_", paste0(j, "t", periods[k]))
#                         
#                         
#                         each_phylo = each_depth[[k]]
#                         output_phylo(each_phylo,
#                                      otuoutfmt = "tsv",
#                                      outfilenames = c(otu_out, tax_out, meta_out),
#                                      otuout_transpose = T,
#                                      methods = "compositional",
#                                      outputdir = out,
#                                      detection = 0.01/100,
#                                      prevalence = 70/100)
#                 }
#                 
#         }
# }
# 
# 
# #####################################
# #ITS
# #####################################
# 
# # 
# ############ d10
# 
# bstITS.C.d10t0 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "10cm" & period == "0wk")
# bstITS.C.d10t1 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "10cm" & period == "1wk")
# bstITS.C.d10t3 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "10cm" & period == "3wk")
# bstITS.C.d10t6 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "10cm" & period == "6wk")
# bstITS.O.d10t0 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "10cm" & period == "0wk")
# bstITS.O.d10t1 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "10cm" & period == "1wk")
# bstITS.O.d10t3 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "10cm" & period == "3wk")
# bstITS.O.d10t6 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "10cm" & period == "6wk")
# bstITS.G.d10t0 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "10cm" & period == "0wk")
# bstITS.G.d10t1 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "10cm" & period == "1wk")
# bstITS.G.d10t3 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "10cm" & period == "3wk")
# bstITS.G.d10t6 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "10cm" & period == "6wk")
# 
# ###############d20
# 
# bstITS.C.d20t0 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "20cm" & period == "0wk")
# bstITS.C.d20t1 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "20cm" & period == "1wk")
# bstITS.C.d20t3 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "20cm" & period == "3wk")
# bstITS.C.d20t6 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "20cm" & period == "6wk")
# bstITS.O.d20t0 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "20cm" & period == "0wk")
# bstITS.O.d20t1 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "20cm" & period == "1wk")
# bstITS.O.d20t3 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "20cm" & period == "3wk")
# bstITS.O.d20t6 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "20cm" & period == "6wk")
# bstITS.G.d20t0 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "20cm" & period == "0wk")
# bstITS.G.d20t1 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "20cm" & period == "1wk")
# bstITS.G.d20t3 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "20cm" & period == "3wk")
# bstITS.G.d20t6 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "20cm" & period == "6wk")
# 
# ############### d30
# 
# bstITS.C.d30t0 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "30cm" & period == "0wk")
# bstITS.C.d30t1 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "30cm" & period == "1wk")
# bstITS.C.d30t3 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "30cm" & period == "3wk")
# bstITS.C.d30t6 = subset_samples(bstITS, cul_type=="Conventional" & soil_depth == "30cm" & period == "6wk")
# bstITS.O.d30t0 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "30cm" & period == "0wk")
# bstITS.O.d30t1 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "30cm" & period == "1wk")
# bstITS.O.d30t3 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "30cm" & period == "3wk")
# bstITS.O.d30t6 = subset_samples(bstITS, cul_type=="Organic" & soil_depth == "30cm" & period == "6wk")
# bstITS.G.d30t0 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "30cm" & period == "0wk")
# bstITS.G.d30t1 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "30cm" & period == "1wk")
# bstITS.G.d30t3 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "30cm" & period == "3wk")
# bstITS.G.d30t6 = subset_samples(bstITS, cul_type=="Green" & soil_depth == "30cm" & period == "6wk")
# 
# 
# convits = list(d10 = list(bstITS.C.d10t0, bstITS.C.d10t1, bstITS.C.d10t3, bstITS.C.d10t6),
#                d20 = list(bstITS.C.d20t0, bstITS.C.d20t1, bstITS.C.d20t3, bstITS.C.d20t6),
#                d30 = list(bstITS.C.d30t0, bstITS.C.d30t1, bstITS.C.d30t3, bstITS.C.d30t6))
# 
# orgits = list(d10 = list(bstITS.O.d10t0, bstITS.O.d10t1, bstITS.O.d10t3, bstITS.O.d10t6),
#               d20 = list(bstITS.O.d20t0, bstITS.O.d20t1, bstITS.O.d20t3, bstITS.O.d20t6),
#               d30 = list(bstITS.O.d30t0, bstITS.O.d30t1, bstITS.O.d30t3, bstITS.O.d30t6))
# 
# grnits = list(d10 = list(bstITS.G.d10t0, bstITS.G.d10t1, bstITS.G.d10t3, bstITS.G.d10t6),
#               d20 = list(bstITS.G.d20t0, bstITS.G.d20t1, bstITS.G.d20t3, bstITS.G.d20t6),
#               d30 = list(bstITS.G.d30t0, bstITS.G.d30t1, bstITS.G.d30t3, bstITS.G.d30t6))
# 
# cultypesITS =  list(Conventional = convits, Organic=orgits, Green=grnits)
# 
# periods = c("t0", "t1", "t3", "t6")
# 
# for(i in c("Conventional", "Organic", "Green")){
#         for(j in c("d10", "d20", "d30")){
#                 
#                 
#                 each_depth = cultypesITS[[i]][[j]]
#                 
#                 dir.create(file.path(outits, paste0("datasource/", i, "/", j)), recursive = T)
#                 
#                 out = file.path(outits, paste0("datasource/", i, "/", j))
#                 
#                 for(k in 1:length(each_depth)){
#                         
#                         otu_out = paste0("otu_", substr(i,1,1), "_", paste0(j, "t", periods[k]))
#                         tax_out = paste0("tax_", substr(i,1,1), "_", paste0(j, "t", periods[k]))
#                         meta_out = paste0("meta_", substr(i,1,1), "_", paste0(j, "t", periods[k]))
#                         
#                         
#                         each_phylo = each_depth[[k]]
#                         
#                         output_phylo(each_phylo,
#                                      otuoutfmt = "tsv",
#                                      outfilenames = c(otu_out, tax_out, meta_out),
#                                      otuout_transpose = T,
#                                      methods = "compositional",
#                                      outputdir = out,
#                                      detection = 0.01/100,
#                                      prevalence = 70/100)
#                 }
#                 
#         }
# }
# 
# ########################
# # change ASV label with time stamp
# ######################
# 
# 
# 
# b.C.d10t0 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_otu_C_d10tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_meta_C_d10tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_tax_C_d10tt0.tsv","microcosms","t0")
# b.C.d10t1 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_otu_C_d10tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_meta_C_d10tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_tax_C_d10tt1.tsv","microcosms", "t1")
# b.C.d10t3 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_otu_C_d10tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_meta_C_d10tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_tax_C_d10tt3.tsv","microcosms", "t3")
# b.C.d10t6 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_otu_C_d10tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_meta_C_d10tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d10/phylo_out/core/CORE_tax_C_d10tt6.tsv","microcosms", "t6")
# b.C.d20t0 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_otu_C_d20tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_meta_C_d20tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_tax_C_d20tt0.tsv","microcosms", "t0")
# b.C.d20t1 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_otu_C_d20tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_meta_C_d20tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_tax_C_d20tt1.tsv","microcosms", "t1")
# b.C.d20t3 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_otu_C_d20tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_meta_C_d20tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_tax_C_d20tt3.tsv","microcosms", "t3")
# b.C.d20t6 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_otu_C_d20tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_meta_C_d20tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d20/phylo_out/core/CORE_tax_C_d20tt6.tsv","microcosms", "t6")
# b.C.d30t0 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_otu_C_d30tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_meta_C_d30tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_tax_C_d30tt0.tsv","microcosms", "t0")
# b.C.d30t1 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_otu_C_d30tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_meta_C_d30tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_tax_C_d30tt1.tsv","microcosms", "t1")
# b.C.d30t3 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_otu_C_d30tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_meta_C_d30tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_tax_C_d30tt3.tsv","microcosms", "t3")
# b.C.d30t6 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_otu_C_d30tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_meta_C_d30tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Conventional/d30/phylo_out/core/CORE_tax_C_d30tt6.tsv","microcosms", "t6")
# 
# 
# #---------------
# # Organic
# # d10
# b.O.d10t0 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_otu_O_d10tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_meta_O_d10tt0.tsv",
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_tax_O_d10tt0.tsv","microcosms", "t0")
# b.O.d10t1 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_otu_O_d10tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_meta_O_d10tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_tax_O_d10tt1.tsv","microcosms", "t1")
# b.O.d10t3 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_otu_O_d10tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_meta_O_d10tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_tax_O_d10tt3.tsv","microcosms", "t3")
# b.O.d10t6 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_otu_O_d10tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_meta_O_d10tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d10/phylo_out/core/CORE_tax_O_d10tt6.tsv","microcosms", "t6")
# b.O.d20t0 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_otu_O_d20tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_meta_O_d20tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_tax_O_d20tt0.tsv","microcosms", "t0")
# b.O.d20t1 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_otu_O_d20tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_meta_O_d20tt1.tsv",
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_tax_O_d20tt1.tsv","microcosms", "t1")
# b.O.d20t3 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_otu_O_d20tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_meta_O_d20tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_tax_O_d20tt3.tsv","microcosms", "t3")
# b.O.d20t6 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_otu_O_d20tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_meta_O_d20tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d20/phylo_out/core/CORE_tax_O_d20tt6.tsv","microcosms", "t6")
# b.O.d30t0 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_otu_O_d30tt0.tsv",
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_meta_O_d30tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_tax_O_d30tt0.tsv","microcosms", "t0")
# b.O.d30t1 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_otu_O_d30tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_meta_O_d30tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_tax_O_d30tt1.tsv","microcosms", "t1")
# b.O.d30t3 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_otu_O_d30tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_meta_O_d30tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_tax_O_d30tt3.tsv","microcosms", "t3")
# b.O.d30t6 = reAssignOTU("./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_otu_O_d30tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_meta_O_d30tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/16s/datasource/Organic/d30/phylo_out/core/CORE_tax_O_d30tt6.tsv","microcosms", "t6")
# 
# 
# 
# 
# ##############################################
# # ITS
# ##############################################
# 
# # conventional
# # d10
# f.C.d10t0 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_otu_C_d10tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_meta_C_d10tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_tax_C_d10tt0.tsv","microcosms", "t0")
# f.C.d10t1 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_otu_C_d10tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_meta_C_d10tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_tax_C_d10tt1.tsv","microcosms", "t1")
# f.C.d10t3 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_otu_C_d10tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_meta_C_d10tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_tax_C_d10tt3.tsv","microcosms", "t3")
# f.C.d10t6 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_otu_C_d10tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_meta_C_d10tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d10/phylo_out/core/CORE_tax_C_d10tt6.tsv","microcosms", "t6")
# # some typeo in microcosm column causing duplicates, change B506 rows
# f.C.d20t0 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_otu_C_d20tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_meta_C_d20tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_tax_C_d20tt0.tsv","microcosms", "t0")
# f.C.d20t1 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_otu_C_d20tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_meta_C_d20tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_tax_C_d20tt1.tsv","microcosms", "t1")
# f.C.d20t3 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_otu_C_d20tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_meta_C_d20tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_tax_C_d20tt3.tsv","microcosms", "t3")
# f.C.d20t6 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_otu_C_d20tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_meta_C_d20tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d20/phylo_out/core/CORE_tax_C_d20tt6.tsv","microcosms", "t6")
# f.C.d30t0 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_otu_C_d30tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_meta_C_d30tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_tax_C_d30tt0.tsv","microcosms", "t0")
# f.C.d30t1 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_otu_C_d30tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_meta_C_d30tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_tax_C_d30tt1.tsv","microcosms", "t1")
# f.C.d30t3 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_otu_C_d30tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_meta_C_d30tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_tax_C_d30tt3.tsv","microcosms", "t3")
# f.C.d30t6 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_otu_C_d30tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_meta_C_d30tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Conventional/d30/phylo_out/core/CORE_tax_C_d30tt6.tsv","microcosms", "t6")
# 
# 
# #---------------
# # Organic
# # d10
# f.O.d10t0 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_otu_O_d10tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_meta_O_d10tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_tax_O_d10tt0.tsv","microcosms", "t0")
# f.O.d10t1 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_otu_O_d10tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_meta_O_d10tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_tax_O_d10tt1.tsv","microcosms", "t1")
# f.O.d10t3 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_otu_O_d10tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_meta_O_d10tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_tax_O_d10tt3.tsv","microcosms", "t3")
# f.O.d10t6 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_otu_O_d10tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_meta_O_d10tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d10/phylo_out/core/CORE_tax_O_d10tt6.tsv","microcosms", "t6")
# f.O.d20t0 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_otu_O_d20tt0.tsv",
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_meta_O_d20tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_tax_O_d20tt0.tsv","microcosms", "t0")
# f.O.d20t1 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_otu_O_d20tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_meta_O_d20tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_tax_O_d20tt1.tsv", "microcosms", "t1")
# f.O.d20t3 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_otu_O_d20tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_meta_O_d20tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_tax_O_d20tt3.tsv","microcosms", "t3")
# f.O.d20t6 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_otu_O_d20tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_meta_O_d20tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d20/phylo_out/core/CORE_tax_O_d20tt6.tsv","microcosms", "t6")
# f.O.d30t0 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_otu_O_d30tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_meta_O_d30tt0.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_tax_O_d30tt0.tsv","microcosms", "t0")
# f.O.d30t1 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_otu_O_d30tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_meta_O_d30tt1.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_tax_O_d30tt1.tsv","microcosms", "t1")
# f.O.d30t3 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_otu_O_d30tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_meta_O_d30tt3.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_tax_O_d30tt3.tsv", "microcosms", "t3")
# f.O.d30t6 = reAssignOTU("./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_otu_O_d30tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_meta_O_d30tt6.tsv", 
#                         "./microcosm_soil_analysis_output/network/its/datasource/Organic/d30/phylo_out/core/CORE_tax_O_d30tt6.tsv","microcosms", "t6")
# 
# 
# 
# 
# # stack up and make cross domain matrixs
# # conventionl
# cross.C.d10 = stackOTUMatrix(b.C.d10t0, b.C.d10t1, b.C.d10t3, b.C.d10t6,
#                              f.C.d10t0, f.C.d10t1, f.C.d10t3, f.C.d10t6, tax_level = "Phylum")
# 
# 
# cross.C.d20 = stackOTUMatrix(b.C.d20t0, b.C.d20t1, b.C.d20t3, b.C.d20t6,
#                              f.C.d20t0, f.C.d20t1, f.C.d20t3, f.C.d20t6, tax_level = "Phylum")
# 
# cross.C.d30 = stackOTUMatrix(b.C.d30t0, b.C.d30t1, b.C.d30t3, b.C.d30t6,
#                              f.C.d30t0, f.C.d30t1, f.C.d30t3, f.C.d30t6, tax_level = "Phylum")
# 
# # organic
# cross.O.d10 = stackOTUMatrix(b.O.d10t0, b.O.d10t1, b.O.d10t3, b.O.d10t6,
#                              f.O.d10t0, f.O.d10t1, f.O.d10t3, f.O.d10t6, tax_level = "Phylum")
# 
# cross.O.d20 = stackOTUMatrix(b.O.d20t0, b.O.d20t1, b.O.d20t3, b.O.d20t6,
#                              f.O.d20t0, f.O.d20t1, f.O.d20t3, f.O.d20t6, tax_level = "Phylum")
# 
# cross.O.d30 = stackOTUMatrix(b.O.d30t0, b.O.d30t1, b.O.d30t3, b.O.d30t6,
#                              f.O.d30t0, f.O.d30t1, f.O.d30t3, f.O.d30t6, tax_level = "Phylum")
# 
# 
# 
# # C
# cross.C.d10.corr = createCorDF(cross.C.d10, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01)
# cross.C.d20.corr = createCorDF(cross.C.d20, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01)
# cross.C.d30.corr = createCorDF(cross.C.d30, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01)
# 
# # O # invalid times argument
# cross.O.d10.corr = createCorDF(cross.O.d10, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01)
# cross.O.d20.corr = createCorDF(cross.O.d20, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01)
# cross.O.d30.corr = createCorDF(cross.O.d30, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01)
# 
# 
# 
# 
# core_list = list( C_d10 = cross.C.d10.corr,
#                   C_d20 = cross.C.d20.corr,
#                   C_d30 = cross.C.d30.corr,
#                   O_d10 = cross.O.d10.corr,
#                   O_d20 = cross.O.d20.corr,
#                   O_d30 = cross.O.d30.corr)
# 
# 
# multiplot_colors = taxColorScheme(core_list, colorPals = glasbey())$all_taxa_color
# 
# 
# source("./tools/networkTools.R")
# 
# dir.create("./microcosm_soil_analysis_output/network/chord_graph")
# out = "./microcosm_soil_analysis_output/network/chord_graph/"
# 
# # both correlation
# for(i in names(core_list)){
#         corObj = core_list[[i]]
#         cat("\n", i, "\n")
#         svg(file.path(out, paste0("both_chord_", i, ".svg")), width = 5, height = 5, pointsize = 12)
#         cor2chord(corObj, legendCex = 0.7, showRelation = 0, addLegend = F, chord_scale = F, outputLegend = F, combinedLegend = T, taxColorSchemeDF = multiplot_colors)
#         dev.off()
# }
# 
# 
# 
# # for(i in names(core_list)){
# #         corObj = core_list[[i]]
# #         svg(file.path(out, paste0("pos_chord_", i, ".svg")), width = 5, height = 5, pointsize = 12)
# #         cor2chord(corObj, legendCex = 0.7, showRelation = 1, addLegend = F, chord_scale = F, outputLegend = F, combinedLegend = T, taxColorSchemeDF = multiplot_colors)
# #         dev.off()
# #         
# # }
# 
# # plot legends
# svg(file.path(out, "legend_both_chord.svg"), width = 4, height = 8)
# plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# legend("bottomleft", legend = multiplot_colors$taxaName, pch = 15, col = as.character(multiplot_colors$taxaCols), cex = 0.8, bty = "n")
# dev.off()
# 
# 
# 
# #################################
# # igraph
# #################################
# # SUMMARY OF the NETWORK
# ################################
# 
# source("./tools/networkTools.R")
# 
# 
# 
# 
# i_c_d10 = cor2igraph(cross.C.d10.corr)
# i_c_d20 = cor2igraph(cross.C.d20.corr)
# i_c_d30 = cor2igraph(cross.C.d30.corr)
# 
# i_o_d10 = cor2igraph(cross.O.d10.corr)
# i_o_d20 = cor2igraph(cross.O.d20.corr)
# i_o_d30 = cor2igraph(cross.O.d30.corr)
# 
# 
# vcount(i_c_d10) # 37 nodes
# vcount(i_c_d20) # 30 nodes
# vcount(i_c_d30) # 31 nodes
# vcount(i_o_d10) # 13 nodes
# vcount(i_o_d20) # 28 nodes
# vcount(i_o_d30) # 25 nodes
# 
# 
# ecount(i_c_d10) # 208 edges
# ecount(i_c_d20) # 217 edges
# ecount(i_c_d30) # 209 edges
# ecount(i_o_d10) # 9 edges
# ecount(i_o_d20) # 36 edges
# ecount(i_o_d30) # 27 edges
# 
# sum(edge.attributes(i_c_d10)$r > 0) #186
# sum(edge.attributes(i_c_d10)$r < 0) #22
# 
# sum(edge.attributes(i_c_d20)$r > 0) #178
# sum(edge.attributes(i_c_d20)$r < 0) #39
# 
# sum(edge.attributes(i_c_d30)$r > 0) #184
# sum(edge.attributes(i_c_d30)$r < 0) #25
# 
# 
# sum(edge.attributes(i_o_d10)$r > 0) #5
# sum(edge.attributes(i_o_d10)$r < 0) #4
# 
# sum(edge.attributes(i_o_d20)$r > 0) #28
# sum(edge.attributes(i_o_d20)$r < 0) #8
# 
# sum(edge.attributes(i_o_d30)$r > 0) #20
# sum(edge.attributes(i_o_d30)$r < 0) #7
# 
# 
# 
# network_summary_df = data.frame(practices = c(rep("Conventional", 3), rep("Organic", 3)), 
#                                 depth = rep(c("d10", "d20", "d30"), 2), 
#                                 nodes = c(vcount(i_c_d10), 
#                                           vcount(i_c_d20), 
#                                           vcount(i_c_d30), 
#                                           vcount(i_o_d10), 
#                                           vcount(i_o_d20), 
#                                           vcount(i_o_d30)), 
#                                 edges = c(ecount(i_c_d10), 
#                                           ecount(i_c_d20), 
#                                           ecount(i_c_d30), 
#                                           ecount(i_o_d10), 
#                                           ecount(i_o_d20), 
#                                           ecount(i_o_d30)), 
#                                 average_degrees = c(mean(igraph::degree(i_c_d10), na.omit=T), 
#                                                     mean(igraph::degree(i_c_d20),na.omit=T), 
#                                                     mean(igraph::degree(i_c_d30),na.omit=T), 
#                                                     mean(igraph::degree(i_o_d10),na.omit=T), 
#                                                     mean(igraph::degree(i_o_d20),na.omit=T),
#                                                     mean(igraph::degree(i_o_d30),na.omit=T)),
#                                 average_between = c(mean(betweenness(i_c_d10)), 
#                                                     mean(betweenness(i_c_d20)), 
#                                                     mean(betweenness(i_c_d30)), 
#                                                     mean(betweenness(i_o_d10)), 
#                                                     mean(betweenness(i_o_d20)),
#                                                     mean(betweenness(i_o_d30))),
#                                 average_closeness = c(mean(closeness(i_c_d10)), 
#                                                       mean(closeness(i_c_d20)), 
#                                                       mean(closeness(i_c_d30)), 
#                                                       mean(closeness(i_o_d10)), 
#                                                       mean(closeness(i_o_d20)),
#                                                       mean(closeness(i_o_d30))))
# 
# head(network_summary_df)
# 
# t.test(average_degrees ~ practices, data = network_summary_df)
# t.test(average_closeness ~ practices, data = network_summary_df)
# t.test(average_between ~ practices, data = network_summary_df)
# t.test(edges ~ practices, data = network_summary_df)
# t.test(nodes ~ practices, data = network_summary_df)
# 
# 
# write.csv(network_summary_df, file=file.path(out, "network_summary.csv"), quote = F, row.names = F)
# 
# 
# 
# 
# 
# 
# 
# p = par(mfrow = c(2,3))
# hist(igraph::degree(i_c_d10), col = "lightblue", xlab = "Node degrees", ylab = "Frequency", main = "Conventional at 10cm")
# hist(igraph::degree(i_c_d20), col = "lightblue", xlab = "Node degrees", ylab = "Frequency", main = "Conventional at 20cm")
# hist(igraph::degree(i_c_d30), col = "lightblue", xlab = "Node degrees", ylab = "Frequency", main = "Conventional at 30cm")
# hist(igraph::degree(i_o_d10), col = "pink", xlab = "Node degrees", ylab = "Frequency", main = "Organic at 10cm")
# hist(igraph::degree(i_o_d20), col = "pink", xlab = "Node degrees", ylab = "Frequency", main = "Organic at 20cm")
# hist(igraph::degree(i_o_d30), col = "pink", xlab = "Node degrees", ylab = "Frequency", main = "Organic at 30cm")
# par(p)
# 
# # node centrality
# 
# degree_c_d10 = sort(igraph::degree(i_c_d10), decreasing = T)
# degree_c_d20 = sort(igraph::degree(i_c_d20), decreasing = T)
# degree_c_d30 = sort(igraph::degree(i_c_d30), decreasing = T)
# degree_o_d10 = sort(igraph::degree(i_o_d10), decreasing = T)
# degree_o_d20 = sort(igraph::degree(i_o_d20), decreasing = T)
# degree_o_d30 = sort(igraph::degree(i_o_d30), decreasing = T)
# 
# 
# 
# betw_c_d10 = sort(betweenness(i_c_d10), decreasing = T)
# betw_c_d20 = sort(betweenness(i_c_d20), decreasing = T)
# betw_c_d30 = sort(betweenness(i_c_d30), decreasing = T)
# betw_o_d10 = sort(betweenness(i_o_d10), decreasing = T)
# betw_o_d20 = sort(betweenness(i_o_d20), decreasing = T)
# betw_o_d30 = sort(betweenness(i_o_d30), decreasing = T)
# 
# 
# close_c_d10 = sort(closeness(i_c_d10), decreasing = T)
# close_c_d20 = sort(closeness(i_c_d20), decreasing = T)
# close_c_d30 = sort(closeness(i_c_d30), decreasing = T)
# close_o_d10 = sort(closeness(i_o_d10), decreasing = T)
# close_o_d20 = sort(closeness(i_o_d20), decreasing = T)
# close_o_d30 = sort(closeness(i_o_d30), decreasing = T)
# 
# 
# nwk_desc_df = data.frame(
#         list(
#                 
#                 cul_types = c(rep("Conventional", sum(vcount(i_c_d10), vcount(i_c_d20), vcount(i_c_d30))),
#                               rep("Organic", sum(vcount(i_o_d10), vcount(i_o_d20), vcount(i_o_d30)))),
#                 
#                 depths = c(rep("d10", vcount(i_c_d10)), 
#                            rep("d20", vcount(i_c_d20)), 
#                            rep("d30", vcount(i_c_d30)), 
#                            rep("d10", vcount(i_o_d10)), 
#                            rep("d20", vcount(i_o_d20)), 
#                            rep("d30", vcount(i_o_d30))), 
#                 
#                 netwk = c(rep("C_d10", vcount(i_c_d10)), 
#                           rep("C_d20", vcount(i_c_d20)), 
#                           rep("C_d30", vcount(i_c_d30)), 
#                           rep("O_d10", vcount(i_o_d10)), 
#                           rep("O_d20", vcount(i_o_d20)), 
#                           rep("O_d30", vcount(i_o_d30))), 
#                 
#                 degrees_names = c(names(degree_c_d10), 
#                                   names(degree_c_d20),
#                                   names(degree_c_d30),
#                                   names(degree_o_d10),
#                                   names(degree_o_d20),
#                                   names(degree_o_d30)),
#                 
#                 degrees = c(    degree_c_d10, 
#                                 degree_c_d20,
#                                 degree_c_d30,
#                                 degree_o_d10,
#                                 degree_o_d20,
#                                 degree_o_d30),
#                 
#                 
#                 betweenness_names = c(names(betw_c_d10), 
#                                       names(betw_c_d20),
#                                       names(betw_c_d30),
#                                       names(betw_o_d10),
#                                       names(betw_o_d20),
#                                       names(betw_o_d30)),
#                 
#                 betweenness = c(betw_c_d10, 
#                                 betw_c_d20,
#                                 betw_c_d30,
#                                 betw_o_d10,
#                                 betw_o_d20,
#                                 betw_o_d30),
#                 
#                 closeness_names = c(names(close_c_d10), 
#                                     names(close_c_d20),
#                                     names(close_c_d30),
#                                     names(close_o_d10),
#                                     names(close_o_d20),
#                                     names(close_o_d30)),
#                 
#                 closeness =       c(close_c_d10, 
#                                     close_c_d20,
#                                     close_c_d30,
#                                     close_o_d10,
#                                     close_o_d20,
#                                     close_o_d30))
#         
# )
# 
# nwk_desc_df$degree_tp = gsub("^.*_", "", nwk_desc_df$degrees_names)
# nwk_desc_df$between_tp = gsub("^.*_", "", nwk_desc_df$betweenness_names)
# nwk_desc_df$close_tp = gsub("^.*_", "", nwk_desc_df$closeness_names)
# 
# # partitioned by cul_types
# 
# nwk_desc_df %>%
#         group_by(cul_types, degree_tp) %>%
#         summarise(sum(degrees))
# 
# 
# nwk_desc_df %>%
#         group_by(cul_types, between_tp) %>%
#         summarise(sum(betweenness))
# 
# 
# nwk_desc_df %>%
#         group_by(cul_types, close_tp) %>%
#         summarise(sum(closeness))
# 
# 
# 
# 
# # taxa
# 
# head(nwk_desc_df)
# 
# nwk_desc_df %>%
#         select(cul_types, degrees_names, degrees) %>%
#         group_by(cul_types, degrees_names) %>%
#         summarise(sum_degree = sum(degrees)) %>%
#         top_n(5) %>%
#         arrange(-sum_degree)
# 
# 
# 
# library(ggplot2)
# library(dplyr)
# library(tidytext)
# library(stringr)
# # degrees
# 
# dir.create("./microcosm_soil_analysis_output/network/out_summary")
# out_summary = "./microcosm_soil_analysis_output/network/out_summary/"
# 
# svg(file.path(out_summary, "node_degrees_rho08_phylum.svg"), width = 15, height = 15)
# nwk_desc_df %>%
#         dplyr::select(netwk, degrees_names, degrees) %>%
#         group_by(netwk) %>%
#         top_n(20) %>%
#         ungroup() %>%
#         mutate(norder = reorder_within(degrees_names, degrees, netwk)) %>%
#         ggplot(aes(x = norder, y = degrees)) +
#         geom_segment(aes(x=norder, xend = norder, y=0, yend=degrees), color="grey") +
#         geom_point(size = 2, alpha=0.7, fill = "springgreen4", shape = 21) +
#         facet_wrap( ~ netwk, scales = "free") +
#         theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank()) +
#         scale_x_reordered() +
#         labs(y = "Nodes degrees", x = "Nodes", title="Node degrees (top20) for co-occurrence networks") +
#         coord_flip()
# dev.off()        
# 
# 
# 
# svg(file.path(out_summary, "node_betweenness_rho08_phylum.svg"), width = 15, height = 15)
# nwk_desc_df %>%
#         dplyr::select(netwk, betweenness_names, betweenness) %>%
#         group_by(netwk) %>%
#         top_n(20) %>%
#         ungroup() %>%
#         mutate(norder = reorder_within(betweenness_names, betweenness, netwk)) %>%
#         ggplot(aes(x = norder, y = betweenness)) +
#         geom_segment(aes(x=norder, xend = norder, y=0, yend=betweenness), color="grey") +
#         geom_point(size = 2, alpha=0.7, fill = "salmon", shape = 21) +
#         facet_wrap( ~ netwk, scales = "free") +
#         theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank()) +
#         scale_x_reordered() +
#         labs(y = "Nodes betweenness", x = "Nodes", title="Node betweenness (top20) for co-occurrence networks") +
#         coord_flip()
# dev.off()
# 
# 
# 
# svg(file.path(out_summary, "node_closeness_rho08_phylum.svg"), width = 15, height = 15)
# nwk_desc_df %>%
#         dplyr::select(netwk, closeness_names, closeness) %>%
#         group_by(netwk) %>%
#         top_n(20) %>%
#         ungroup() %>%
#         mutate(norder = reorder_within(closeness_names, closeness, netwk)) %>%
#         ggplot(aes(x = norder, y = closeness)) +
#         geom_segment(aes(x=norder, xend = norder, y=0, yend=closeness), color="grey") +
#         geom_point(size = 2, alpha=0.7, fill = "dodgerblue", shape = 21) +
#         facet_wrap( ~ netwk, scales = "free") +
#         theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank()) +
#         scale_x_reordered() +
#         labs(y = "Nodes closeness", x = "Nodes", title="Node closeness (top20) for co-occurrence networks") +
#         coord_flip()
# dev.off()
# 
# 
# 
# 
# 
# #PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
# 
# 
# ###############################################################
# # family level
# ###############################################################
# # conventionl
# 
# cross.C.d10.f = stackOTUMatrix(b.C.d10t0, b.C.d10t1, b.C.d10t3, b.C.d10t6,
#                                f.C.d10t0, f.C.d10t1, f.C.d10t3, f.C.d10t6, tax_level = "Family")
# cross.C.d20.f = stackOTUMatrix(b.C.d20t0, b.C.d20t1, b.C.d20t3, b.C.d20t6,
#                                f.C.d20t0, f.C.d20t1, f.C.d20t3, f.C.d20t6, tax_level = "Family")
# cross.C.d30.f = stackOTUMatrix(b.C.d30t0, b.C.d30t1, b.C.d30t3, b.C.d30t6,
#                                f.C.d30t0, f.C.d30t1, f.C.d30t3, f.C.d30t6, tax_level = "Family")
# 
# 
# # organic
# cross.O.d10.f = stackOTUMatrix(b.O.d10t0, b.O.d10t1, b.O.d10t3, b.O.d10t6,
#                                f.O.d10t0, f.O.d10t1, f.O.d10t3, f.O.d10t6, tax_level = "Family")
# cross.O.d20.f = stackOTUMatrix(b.O.d20t0, b.O.d20t1, b.O.d20t3, b.O.d20t6,
#                                f.O.d20t0, f.O.d20t1, f.O.d20t3, f.O.d20t6, tax_level = "Family")
# cross.O.d30.f = stackOTUMatrix(b.O.d30t0, b.O.d30t1, b.O.d30t3, b.O.d30t6,
#                                f.O.d30t0, f.O.d30t1, f.O.d30t3, f.O.d30t6, tax_level = "Family")
# 
# 
# 
# 
# 
# cross.C.d10.corr.f = createCorDF(cross.C.d10.f, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Family")
# cross.C.d20.corr.f = createCorDF(cross.C.d20.f, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Family")
# cross.C.d30.corr.f = createCorDF(cross.C.d30.f, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Family")
# cross.O.d10.corr.f = createCorDF(cross.O.d10.f, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Family")
# cross.O.d20.corr.f = createCorDF(cross.O.d20.f, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Family")
# cross.O.d30.corr.f = createCorDF(cross.O.d30.f, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Family")
# 
# 
# 
# 
# core_list.f = list( C_d10 = cross.C.d10.corr.f,
#                     C_d20 = cross.C.d20.corr.f,
#                     C_d30 = cross.C.d30.corr.f,
#                     O_d10 = cross.O.d10.corr.f,
#                     O_d20 = cross.O.d20.corr.f,
#                     O_d30 = cross.O.d30.corr.f)
# 
# 
# 
# multiplot_colors.f = taxColorScheme(core_list.f)
# sector_colors.f = multiplot_colors.f$all_sector_colors
# 
# 
# # plot legend
# 
# dir.create("./microcosm_soil_analysis_output/network/chord_graph/fam")
# 
# out = "./microcosm_soil_analysis_output/network/chord_graph/fam/"
# 
# # svg(file.path(out, "legend_both_chord_fam.svg"), width = 4, height = 20)
# # plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# # legend("bottomleft", legend = sector_colors.f$legendNames, pch = 15, col = as.character(sector_colors.f$legendColors), cex = 0.8, bty = "n", ncol = 3)
# # dev.off()
# 
# 
# 
# 
# # # legend.list = list()
# # # positive
# 
# 
# for(i in names(core_list.f)){
#         corObj.f = core_list.f[[i]]
#         svg(file.path(out, paste0("fam_pos_chord_", i, ".svg")), width = 5, height = 5, pointsize = 12)
#         L= cor2chord(corObj.f, legendCex = 0.7, showRelation = 1, addLegend = F, chord_scale = T, outputLegend = T, combinedLegend = F, plot_title = i)
#         dev.off()
# 
#         svg(file.path(out, paste("legend_fam_pos_chord_", i, ".svg")), width = 12, height = 15, pointsize = 13)
#         plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#         legend("bottomleft", legend = L$legendNames, pch = 15, col = as.character(L$legendColors), cex = 0.8, bty = "n", ncol = 3)
#         dev.off()
# 
# }
# # 
# # # # negative
# for(i in names(core_list.f)){
#         corObj.f = core_list.f[[i]]
#         svg(file.path(out, paste0("fam_neg_chord_", i, ".svg")), width = 5, height = 5, pointsize = 12)
#         L = cor2chord(corObj.f, legendCex = 0.7, showRelation = -1, addLegend = F, chord_scale = T, outputLegend = T, combinedLegend = F, plot_title = i)
#         dev.off()
# 
#         svg(file.path(out, paste("legend_fam_neg_chord_", i, ".svg")), width = 12, height = 15, pointsize = 13)
#         plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#         legend("bottomleft", legend = L$legendNames, pch = 15, col = as.character(L$legendColors), cex = 0.8, bty = "n", ncol = 3)
#         dev.off()
# 
# }
# 
# # both correlation
# 
# # legend.list = list()
# for(i in names(core_list.f)){
#         corObj.f = core_list.f[[i]]
#         svg(file.path(out, paste0("fam_both_chord_", i, ".svg")), width = 5, height = 5, pointsize = 12)
#         L = cor2chord(corObj.f, legendCex = 0.7, showRelation = 0, addLegend = F, chord_scale = T, outputLegend = T, combinedLegend = F, plot_title = i)
#         dev.off()
# 
#         svg(file.path(out, paste("both_legend_fam_chord_", i, ".svg")), width = 12, height = 15, pointsize = 13)
#         plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#         legend("bottomleft", legend = L$legendNames, pch = 15, col = as.character(L$legendColors), cex = 0.8, bty = "n", ncol = 3)
#         dev.off()
# }
# 
# 
# 
# #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 
# 
# 
# i_c_d10.f = cor2igraph(cross.C.d10.corr.f)
# i_c_d20.f = cor2igraph(cross.C.d20.corr.f)
# i_c_d30.f = cor2igraph(cross.C.d30.corr.f)
# i_o_d10.f = cor2igraph(cross.O.d10.corr.f)
# i_o_d20.f = cor2igraph(cross.O.d20.corr.f)
# i_o_d30.f = cor2igraph(cross.O.d30.corr.f)
# 
# 
# 
# 
# vcount(i_c_d10.f) # 132 nodes
# vcount(i_c_d20.f) # 121 nodes
# vcount(i_c_d30.f) # 126 nodes
# vcount(i_o_d10.f) # 15 nodes
# vcount(i_o_d20.f) # 65 nodes
# vcount(i_o_d30.f) # 44 nodes
# ecount(i_c_d10.f) # 208 edges
# ecount(i_c_d20.f) # 217 edges
# ecount(i_c_d30.f) # 209 edges
# ecount(i_o_d10.f) # 9 edges
# ecount(i_o_d20.f) # 36 edges
# ecount(i_o_d30.f) # 27 edges
# 
# sum(edge.attributes(i_c_d10.f)$r > 0) #186
# sum(edge.attributes(i_c_d10.f)$r < 0) #22
# sum(edge.attributes(i_c_d20.f)$r > 0) #178
# sum(edge.attributes(i_c_d20.f)$r < 0) #39
# sum(edge.attributes(i_c_d30.f)$r > 0) #184
# sum(edge.attributes(i_c_d30.f)$r < 0) #25
# 
# sum(edge.attributes(i_o_d10.f)$r > 0) #5
# sum(edge.attributes(i_o_d10.f)$r < 0) #4
# sum(edge.attributes(i_o_d20.f)$r > 0) #28
# sum(edge.attributes(i_o_d20.f)$r < 0) #8
# sum(edge.attributes(i_o_d30.f)$r > 0) #20
# sum(edge.attributes(i_o_d30.f)$r < 0) #7
# 
# 
# 
# network_summary_df.f = data.frame(practices = c(rep("Conventional", 3), rep("Organic", 3)), 
#                                   depth = rep(c("d10", "d20", "d30"), 2), 
#                                   nodes = c(vcount(i_c_d10.f), 
#                                             vcount(i_c_d20.f), 
#                                             vcount(i_c_d30.f), 
#                                             vcount(i_o_d10.f), 
#                                             vcount(i_o_d20.f), 
#                                             vcount(i_o_d30.f)), 
#                                   
#                                   edges = c(ecount(i_c_d10.f),
#                                             ecount(i_c_d20.f), 
#                                             ecount(i_c_d30.f), 
#                                             ecount(i_o_d10.f), 
#                                             ecount(i_o_d20.f), 
#                                             ecount(i_o_d30.f)), 
#                                   
#                                   average_degrees = c(mean(igraph::degree(i_c_d10.f), na.rm = T), 
#                                                       mean(igraph::degree(i_c_d20.f), na.rm = T), 
#                                                       mean(igraph::degree(i_c_d30.f), na.rm = T), 
#                                                       mean(igraph::degree(i_o_d10.f), na.rm = T), 
#                                                       mean(igraph::degree(i_o_d20.f), na.rm = T),
#                                                       mean(igraph::degree(i_o_d30.f), na.rm = T)),
#                                   
#                                   average_between = c(mean(igraph::betweenness(i_c_d10.f), na.rm = T), 
#                                                       mean(igraph::betweenness(i_c_d20.f), na.rm = T), 
#                                                       mean(igraph::betweenness(i_c_d30.f), na.rm = T), 
#                                                       mean(igraph::betweenness(i_o_d10.f), na.rm = T), 
#                                                       mean(igraph::betweenness(i_o_d20.f), na.rm = T),
#                                                       mean(igraph::betweenness(i_o_d30.f), na.rm = T)),
#                                   
#                                   average_closeness = c(mean(igraph::closeness(i_c_d10.f), na.rm = T), 
#                                                         mean(igraph::closeness(i_c_d20.f), na.rm = T), 
#                                                         mean(igraph::closeness(i_c_d30.f), na.rm = T), 
#                                                         mean(igraph::closeness(i_o_d10.f), na.rm = T), 
#                                                         mean(igraph::closeness(i_o_d20.f), na.rm = T),
#                                                         mean(igraph::closeness(i_o_d30.f), na.rm = T)))
# 
# 
# 
# network_summary_df.f
# 
# 
# 
# 
# t.test(average_degrees ~ practices, data = network_summary_df.f)
# t.test(average_closeness ~ practices, data = network_summary_df.f)
# t.test(average_between ~ practices, data = network_summary_df.f)
# t.test(edges ~ practices, data = network_summary_df.f)
# t.test(nodes ~ practices, data = network_summary_df.f)
# 
# 
# write.csv(network_summary_df.f, file=file.path(out, "network_summary_spearman08_fam.csv"), quote = F, row.names = F)
# 
# 
# 
# p = par(mfrow = c(2,3))
# hist(igraph::degree(i_c_d10.f), col = "lightblue", xlab = "Node degrees", ylab = "Frequency", main = "Conventional at 10cm")
# hist(igraph::degree(i_c_d20.f), col = "lightblue", xlab = "Node degrees", ylab = "Frequency", main = "Conventional at 20cm")
# hist(igraph::degree(i_c_d30.f), col = "lightblue", xlab = "Node degrees", ylab = "Frequency", main = "Conventional at 30cm")
# hist(igraph::degree(i_o_d10.f), col = "pink", xlab = "Node degrees", ylab = "Frequency", main = "Organic at 10cm")
# hist(igraph::degree(i_o_d20.f), col = "pink", xlab = "Node degrees", ylab = "Frequency", main = "Organic at 20cm")
# hist(igraph::degree(i_o_d30.f), col = "pink", xlab = "Node degrees", ylab = "Frequency", main = "Organic at 30cm")
# par(p)
# 
# # node centrality
# 
# degree_c_d10.f = sort(igraph::degree(i_c_d10.f), decreasing = T)
# degree_c_d20.f = sort(igraph::degree(i_c_d20.f), decreasing = T)
# degree_c_d30.f = sort(igraph::degree(i_c_d30.f), decreasing = T)
# degree_o_d10.f = sort(igraph::degree(i_o_d10.f), decreasing = T)
# degree_o_d20.f = sort(igraph::degree(i_o_d20.f), decreasing = T)
# degree_o_d30.f = sort(igraph::degree(i_o_d30.f), decreasing = T)
# 
# 
# 
# betw_c_d10.f = sort(betweenness(i_c_d10.f), decreasing = T)
# betw_c_d20.f = sort(betweenness(i_c_d20.f), decreasing = T)
# betw_c_d30.f = sort(betweenness(i_c_d30.f), decreasing = T)
# betw_o_d10.f = sort(betweenness(i_o_d10.f), decreasing = T)
# betw_o_d20.f = sort(betweenness(i_o_d20.f), decreasing = T)
# betw_o_d30.f = sort(betweenness(i_o_d30.f), decreasing = T)
# 
# 
# close_c_d10.f = sort(closeness(i_c_d10.f), decreasing = T)
# close_c_d20.f = sort(closeness(i_c_d20.f), decreasing = T)
# close_c_d30.f = sort(closeness(i_c_d30.f), decreasing = T)
# close_o_d10.f = sort(closeness(i_o_d10.f), decreasing = T)
# close_o_d20.f = sort(closeness(i_o_d20.f), decreasing = T)
# close_o_d30.f = sort(closeness(i_o_d30.f), decreasing = T)
# 
# 
# nwk_desc_df.f = data.frame(
#         list(
#                 
#                 cul_types = c(rep("Conventional", 
#                                   sum(vcount(i_c_d10.f), 
#                                       vcount(i_c_d20.f), 
#                                       vcount(i_c_d30.f))),
#                               rep("Organic", 
#                                   sum(vcount(i_o_d10.f), 
#                                       vcount(i_o_d20.f), 
#                                       vcount(i_o_d30.f)))),
#                 
#                 depths = c(rep("d10", vcount(i_c_d10.f)), 
#                            rep("d20", vcount(i_c_d20.f)), 
#                            rep("d30", vcount(i_c_d30.f)), 
#                            rep("d10", vcount(i_o_d10.f)), 
#                            rep("d20", vcount(i_o_d20.f)), 
#                            rep("d30", vcount(i_o_d30.f))), 
#                 
#                 netwk = c(rep("C_d10", vcount(i_c_d10.f)), 
#                           rep("C_d20", vcount(i_c_d20.f)), 
#                           rep("C_d30", vcount(i_c_d30.f)), 
#                           rep("O_d10", vcount(i_o_d10.f)), 
#                           rep("O_d20", vcount(i_o_d20.f)), 
#                           rep("O_d30", vcount(i_o_d30.f))), 
#                 
#                 degrees_names = c(names(degree_c_d10.f), 
#                                   names(degree_c_d20.f),
#                                   names(degree_c_d30.f),
#                                   names(degree_o_d10.f),
#                                   names(degree_o_d20.f),
#                                   names(degree_o_d30.f)),
#                 
#                 degrees = c(    degree_c_d10.f, 
#                                 degree_c_d20.f,
#                                 degree_c_d30.f,
#                                 degree_o_d10.f,
#                                 degree_o_d20.f,
#                                 degree_o_d30.f),
#                 
#                 
#                 betweenness_names = c(names(betw_c_d10.f), 
#                                       names(betw_c_d20.f),
#                                       names(betw_c_d30.f),
#                                       names(betw_o_d10.f),
#                                       names(betw_o_d20.f),
#                                       names(betw_o_d30.f)),
#                 
#                 betweenness = c(betw_c_d10.f, 
#                                 betw_c_d20.f,
#                                 betw_c_d30.f,
#                                 betw_o_d10.f,
#                                 betw_o_d20.f,
#                                 betw_o_d30.f),
#                 
#                 closeness_names = c(names(close_c_d10.f), 
#                                     names(close_c_d20.f),
#                                     names(close_c_d30.f),
#                                     names(close_o_d10.f),
#                                     names(close_o_d20.f),
#                                     names(close_o_d30.f)),
#                 
#                 closeness =       c(close_c_d10.f, 
#                                     close_c_d20.f,
#                                     close_c_d30.f,
#                                     close_o_d10.f,
#                                     close_o_d20.f,
#                                     close_o_d30.f))
#         
# )
# 
# # determine which time point
# 
# head(nwk_desc_df.f)
# 
# nwk_desc_df.f$degree_tp = gsub("^.*_", "", nwk_desc_df.f$degrees_names)
# nwk_desc_df.f$between_tp = gsub("^.*_", "", nwk_desc_df.f$betweenness_names)
# nwk_desc_df.f$close_tp = gsub("^.*_", "", nwk_desc_df.f$closeness_names)
# 
# 
# # partitioned by cul_types
# 
# nwk_desc_df.f %>%
#         group_by(cul_types, degree_tp) %>%
#         summarise(sum(degrees))
# 
# 
# nwk_desc_df.f %>%
#         group_by(cul_types, between_tp) %>%
#         summarise(sum(betweenness))
# 
# 
# nwk_desc_df.f %>%
#         group_by(cul_types, close_tp) %>%
#         summarise(sum(closeness))
# 
# 
# # taxa
# 
# nwk_desc_df.f %>%
#         select(cul_types, degrees_names, degrees) %>%
#         group_by(cul_types, degrees_names) %>%
#         summarise(sum_degree = sum(degrees)) %>%
#         top_n(5) %>%
#         arrange(-sum_degree)
# 
# 
# svg(file.path(out_summary, "node_degrees_rho08_fam.svg"), width = 15, height = 15)
# nwk_desc_df.f %>%
#         dplyr::select(netwk, degrees_names, degrees) %>%
#         group_by(netwk) %>%
#         top_n(20) %>%
#         ungroup() %>%
#         mutate(norder = reorder_within(degrees_names, degrees, netwk)) %>%
#         ggplot(aes(x = norder, y = degrees)) +
#         geom_segment(aes(x=norder, xend = norder, y=0, yend=degrees), color="grey") +
#         geom_point(size = 2, alpha=0.7, fill = "springgreen4", shape = 21) +
#         facet_wrap( ~ netwk, scales = "free") +
#         theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank()) +
#         scale_x_reordered() +
#         labs(y = "Nodes degrees", x = "Nodes", title="Node degrees (top 20) for co-occurrence networks") +
#         coord_flip()
# 
# dev.off()
# 
# 
# 
# 
# 
# svg(file.path(out_summary, "node_betweenness_rho08_fam.svg"), width = 15, height = 15)
# nwk_desc_df.f %>%
#         dplyr::select(netwk, betweenness_names, betweenness) %>%
#         group_by(netwk) %>%
#         top_n(20) %>%
#         ungroup() %>%
#         mutate(norder = reorder_within(betweenness_names, betweenness, netwk)) %>%
#         ggplot(aes(x = norder, y = betweenness)) +
#         geom_segment(aes(x=norder, xend = norder, y=0, yend=betweenness), color="grey") +
#         geom_point(size = 2, alpha=0.7, fill = "salmon", shape = 21) +
#         facet_wrap( ~ netwk, scales = "free") +
#         theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank()) +
#         scale_x_reordered() +
#         labs(y = "Nodes betweenness", x = "Nodes", title="Node betweenness (top 20) for co-occurrence networks") +
#         coord_flip()
# dev.off()
# 
# 
# 
# svg(file.path(out_summary, "node_closeness_rho08_fam.svg"), width = 15, height = 15)
# nwk_desc_df.f %>%
#         dplyr::select(netwk, closeness_names, closeness) %>%
#         group_by(netwk) %>%
#         top_n(20) %>%
#         ungroup() %>%
#         mutate(norder = reorder_within(closeness_names, closeness, netwk)) %>%
#         ggplot(aes(x = norder, y = closeness)) +
#         geom_segment(aes(x=norder, xend = norder, y=0, yend=closeness), color="grey") +
#         geom_point(size = 2, alpha=0.7, fill = "dodgerblue", shape = 21) +
#         facet_wrap( ~ netwk, scales = "free") +
#         theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank()) +
#         scale_x_reordered() +
#         labs(y = "Nodes closeness", x = "Nodes", title="Node closeness (top 20) for co-occurrence networks") +
#         coord_flip()
# 
# dev.off()
# 
# 
# 
# ######################################################################################################################
# # genus level
# ######################################################################################################################
# 
# # conventionl
# 
# cross.C.d10.g = stackOTUMatrix(b.C.d10t0, b.C.d10t1, b.C.d10t3, b.C.d10t6,
#                                f.C.d10t0, f.C.d10t1, f.C.d10t3, f.C.d10t6, tax_level = "Genus")
# cross.C.d20.g = stackOTUMatrix(b.C.d20t0, b.C.d20t1, b.C.d20t3, b.C.d20t6,
#                                f.C.d20t0, f.C.d20t1, f.C.d20t3, f.C.d20t6, tax_level = "Genus")
# cross.C.d30.g = stackOTUMatrix(b.C.d30t0, b.C.d30t1, b.C.d30t3, b.C.d30t6,
#                                f.C.d30t0, f.C.d30t1, f.C.d30t3, f.C.d30t6, tax_level = "Genus")
# 
# 
# # organic
# cross.O.d10.g = stackOTUMatrix(b.O.d10t0, b.O.d10t1, b.O.d10t3, b.O.d10t6,
#                                f.O.d10t0, f.O.d10t1, f.O.d10t3, f.O.d10t6, tax_level = "Genus")
# cross.O.d20.g = stackOTUMatrix(b.O.d20t0, b.O.d20t1, b.O.d20t3, b.O.d20t6,
#                                f.O.d20t0, f.O.d20t1, f.O.d20t3, f.O.d20t6, tax_level = "Genus")
# cross.O.d30.g = stackOTUMatrix(b.O.d30t0, b.O.d30t1, b.O.d30t3, b.O.d30t6,
#                                f.O.d30t0, f.O.d30t1, f.O.d30t3, f.O.d30t6, tax_level = "Genus")
# 
# 
# 
# 
# 
# cross.C.d10.corr.g = createCorDF(cross.C.d10.g, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Genus")
# cross.C.d20.corr.g = createCorDF(cross.C.d20.g, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Genus")
# cross.C.d30.corr.g = createCorDF(cross.C.d30.g, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Genus")
# 
# cross.O.d10.corr.g = createCorDF(cross.O.d10.g, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Genus")
# cross.O.d20.corr.g = createCorDF(cross.O.d20.g, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Genus")
# cross.O.d30.corr.g = createCorDF(cross.O.d30.g, corTest = "rcorr", padj_method = "BH", cor_cutoff = 0.8, padj_cutoff = 0.01, taxlevel = "Genus")
# 
# 
# 
# 
# core_list.g = list( C_d10 = cross.C.d10.corr.g,
#                     C_d20 = cross.C.d20.corr.g,
#                     C_d30 = cross.C.d30.corr.g,
#                     O_d10 = cross.O.d10.corr.g,
#                     O_d20 = cross.O.d20.corr.g,
#                     O_d30 = cross.O.d30.corr.g)
# 
# 
# 
# multiplot_colors.g = taxColorScheme(core_list.g)
# sector_colors.g = multiplot_colors.g$all_sector_colors
# 
# 
# # plot legend
# 
# dir.create("./microcosm_soil_analysis_output/network/chord_graph/gen")
# out = "./microcosm_soil_analysis_output/network/chord_graph/gen"
# 
# svg(file.path(out, "legend_both_chord_genus.svg"), width = 4, height = 20)
# plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# legend("bottomleft", legend = sector_colors.g$legendNames, pch = 15, col = as.character(sector_colors.g$legendColors), cex = 0.8, bty = "n", ncol = 3)
# dev.off()
# 
# 
# 
# # # legend.list = list()
# # # positive
# 
# # 
# # for(i in names(core_list.g)){
# #         corObj.g = core_list.g[[i]]
# #         svg(file.path(out, paste0("genus_pos_chord_", i, ".svg")), width = 5, height = 5, pointsize = 12)
# #         L= cor2chord(corObj.g, legendCex = 0.7, showRelation = 1, addLegend = F, chord_scale = F, outputLegend = T, combinedLegend = F, plot_title = i)
# #         dev.off()
# # 
# #         svg(file.path(out, paste("legend_genus_pos_chord_", i, ".svg")), width = 12, height = 15, pointsize = 13)
# #         plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# #         legend("bottomleft", legend = L$legendNames, pch = 15, col = as.character(L$legendColors), cex = 0.8, bty = "n", ncol = 3)
# #         dev.off()
# # 
# # }
# # 
# # # # negative
# # for(i in names(core_list.g)){
# #         corObj.g = core_list.g[[i]]
# #         svg(file.path(out, paste0("genus_neg_chord_", i, ".svg")), width = 5, height = 5, pointsize = 12)
# #         L = cor2chord(corObj.g, legendCex = 0.7, showRelation = -1, addLegend = F, chord_scale = F, outputLegend = T, combinedLegend = F, plot_title = i)
# #         dev.off()
# # 
# #         svg(file.path(out, paste("legend_genus_neg_chord_", i, ".svg")), width = 12, height = 15, pointsize = 13)
# #         plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# #         legend("bottomleft", legend = L$legendNames, pch = 15, col = as.character(L$legendColors), cex = 0.8, bty = "n", ncol = 3)
# #         dev.off()
# # 
# # }
# 
# # # both correlation
# 
# # legend.list = list()
# for(i in names(core_list.g)){
#         corObj.g = core_list.g[[i]]
#         svg(file.path(out, paste0("genus_both_chord_", i, ".svg")), width = 10, height = 10)
#         L = cor2chord(corObj.g, legendCex = 0.7, showRelation = 0, addLegend = F, chord_scale = F, outputLegend = T, combinedLegend = F, plot_title = i)
#         dev.off()
# 
#         svg(file.path(out, paste("legend_genus_chord_", i, ".svg")), width = 10, height = 15)
#         plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#         legend("bottomleft", legend = L$legendNames, pch = 15, col = as.character(L$legendColors), cex = 0.8, bty = "n", ncol = 3)
#         dev.off()
# }
# 
# svg(file.path(out, paste("legend_genus_chord_", "test", ".svg")), width = 10, height = 15)
# plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# legend("bottomleft", legend = L$legendNames, pch = 15, col = as.character(L$legendColors), cex = 0.8, bty = "n", ncol = 3, text.width = c(0.2,0.2,0.2))
# dev.off()
# 
# #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 
# 
# 
# i_c_d10.g = cor2igraph(cross.C.d10.corr.g)
# i_c_d20.g = cor2igraph(cross.C.d20.corr.g)
# i_c_d30.g = cor2igraph(cross.C.d30.corr.g)
# 
# i_o_d10.g = cor2igraph(cross.O.d10.corr.g)
# i_o_d20.g = cor2igraph(cross.O.d20.corr.g)
# i_o_d30.g = cor2igraph(cross.O.d30.corr.g)
# 
# 
# 
# 
# vcount(i_c_d10.g) # 148 nodes
# vcount(i_c_d20.g) # 134 nodes
# vcount(i_c_d30.g) # 141 nodes
# vcount(i_o_d10.g) # 16 nodes
# vcount(i_o_d20.g) # 63 nodes
# vcount(i_o_d30.g) # 46 nodes
# 
# 
# ecount(i_c_d10.g) # 208 edges
# ecount(i_c_d20.g) # 217 edges
# ecount(i_c_d30.g) # 209 edges
# ecount(i_o_d10.g) # 9 edges
# ecount(i_o_d20.g) # 36 edges
# ecount(i_o_d30.g) # 27 edges
# 
# sum(edge.attributes(i_c_d10.g)$r > 0) #186
# sum(edge.attributes(i_c_d10.g)$r < 0) #22
# 
# sum(edge.attributes(i_c_d20.g)$r > 0) #178
# sum(edge.attributes(i_c_d20.g)$r < 0) #39
# 
# sum(edge.attributes(i_c_d30.g)$r > 0) #184
# sum(edge.attributes(i_c_d30.g)$r < 0) #25
# 
# 
# sum(edge.attributes(i_o_d10.g)$r > 0) #5
# sum(edge.attributes(i_o_d10.g)$r < 0) #4
# 
# sum(edge.attributes(i_o_d20.g)$r > 0) #28
# sum(edge.attributes(i_o_d20.g)$r < 0) #8
# 
# sum(edge.attributes(i_o_d30.g)$r > 0) #20
# sum(edge.attributes(i_o_d30.g)$r < 0) #7
# 
# 
# 
# network_summary_df.g = data.frame(practices = c(rep("Conventional", 3), rep("Organic", 3)),
#                                   depth = rep(c("d10", "d20", "d30"), 2),
#                                   nodes = c(vcount(i_c_d10.g),
#                                             vcount(i_c_d20.g),
#                                             vcount(i_c_d30.g),
#                                             vcount(i_o_d10.g),
#                                             vcount(i_o_d20.g),
#                                             vcount(i_o_d30.g)),
# 
#                                   edges = c(ecount(i_c_d10.g),
#                                             ecount(i_c_d20.g),
#                                             ecount(i_c_d30.g),
#                                             ecount(i_o_d10.g),
#                                             ecount(i_o_d20.g),
#                                             ecount(i_o_d30.g)),
# 
#                                   average_degrees = c(mean(igraph::degree(i_c_d10.g), na.rm = T),
#                                                       mean(igraph::degree(i_c_d20.g), na.rm = T),
#                                                       mean(igraph::degree(i_c_d30.g), na.rm = T),
#                                                       mean(igraph::degree(i_o_d10.g), na.rm = T),
#                                                       mean(igraph::degree(i_o_d20.g), na.rm = T),
#                                                       mean(igraph::degree(i_o_d30.g), na.rm = T)),
# 
#                                   average_between = c(mean(igraph::betweenness(i_c_d10.g), na.rm = T),
#                                                       mean(igraph::betweenness(i_c_d20.g), na.rm = T),
#                                                       mean(igraph::betweenness(i_c_d30.g), na.rm = T),
#                                                       mean(igraph::betweenness(i_o_d10.g), na.rm = T),
#                                                       mean(igraph::betweenness(i_o_d20.g), na.rm = T),
#                                                       mean(igraph::betweenness(i_o_d30.g), na.rm = T)),
# 
#                                   average_closeness = c(mean(igraph::closeness(i_c_d10.g), na.rm = T),
#                                                         mean(igraph::closeness(i_c_d20.g), na.rm = T),
#                                                         mean(igraph::closeness(i_c_d30.g), na.rm = T),
#                                                         mean(igraph::closeness(i_o_d10.g), na.rm = T),
#                                                         mean(igraph::closeness(i_o_d20.g), na.rm = T),
#                                                         mean(igraph::closeness(i_o_d30.g), na.rm = T)))
# 
# 
# 
# network_summary_df.g
# 
# 
# 
# 
# summary(aov(average_degrees ~ practices, data = network_summary_df.g))
# summary(aov(average_closeness ~ practices, data = network_summary_df.g))
# summary(aov(average_between ~ practices, data = network_summary_df.g))
# summary(aov(edges ~ practices, data = network_summary_df.g))
# summary(aov(nodes ~ practices, data = network_summary_df.g))
# 
# 
# write.csv(network_summary_df.g, file=file.path(out, "network_summary_rho08_genus.csv"), quote = F, row.names = F)
# 
# 
# 
# 
# 
# p = par(mfrow = c(2,3))
# hist(igraph::degree(i_c_d10.g), col = "lightblue", xlab = "Node degrees", ylab = "Frequency", main = "Conventional at 10cm")
# hist(igraph::degree(i_c_d20.g), col = "lightblue", xlab = "Node degrees", ylab = "Frequency", main = "Conventional at 20cm")
# hist(igraph::degree(i_c_d30.g), col = "lightblue", xlab = "Node degrees", ylab = "Frequency", main = "Conventional at 30cm")
# hist(igraph::degree(i_o_d10.g), col = "pink", xlab = "Node degrees", ylab = "Frequency", main = "Organic at 10cm")
# hist(igraph::degree(i_o_d20.g), col = "pink", xlab = "Node degrees", ylab = "Frequency", main = "Organic at 20cm")
# hist(igraph::degree(i_o_d30.g), col = "pink", xlab = "Node degrees", ylab = "Frequency", main = "Organic at 30cm")
# par(p)
# 
# # node centrality
# 
# degree_c_d10.g = sort(igraph::degree(i_c_d10.g), decreasing = T)
# degree_c_d20.g = sort(igraph::degree(i_c_d20.g), decreasing = T)
# degree_c_d30.g = sort(igraph::degree(i_c_d30.g), decreasing = T)
# degree_o_d10.g = sort(igraph::degree(i_o_d10.g), decreasing = T)
# degree_o_d20.g = sort(igraph::degree(i_o_d20.g), decreasing = T)
# degree_o_d30.g = sort(igraph::degree(i_o_d30.g), decreasing = T)
# 
# 
# 
# betw_c_d10.g = sort(betweenness(i_c_d10.g), decreasing = T)
# betw_c_d20.g = sort(betweenness(i_c_d20.g), decreasing = T)
# betw_c_d30.g = sort(betweenness(i_c_d30.g), decreasing = T)
# betw_o_d10.g = sort(betweenness(i_o_d10.g), decreasing = T)
# betw_o_d20.g = sort(betweenness(i_o_d20.g), decreasing = T)
# betw_o_d30.g = sort(betweenness(i_o_d30.g), decreasing = T)
# 
# 
# close_c_d10.g = sort(closeness(i_c_d10.g), decreasing = T)
# close_c_d20.g = sort(closeness(i_c_d20.g), decreasing = T)
# close_c_d30.g = sort(closeness(i_c_d30.g), decreasing = T)
# close_o_d10.g = sort(closeness(i_o_d10.g), decreasing = T)
# close_o_d20.g = sort(closeness(i_o_d20.g), decreasing = T)
# close_o_d30.g = sort(closeness(i_o_d30.g), decreasing = T)
# 
# 
# nwk_desc_df.g = data.frame(
#         list(
# 
#                 cul_types = c(rep("Conventional",
#                                   sum(vcount(i_c_d10.g),
#                                       vcount(i_c_d20.g),
#                                       vcount(i_c_d30.g))),
#                               rep("Organic",
#                                   sum(vcount(i_o_d10.g),
#                                       vcount(i_o_d20.g),
#                                       vcount(i_o_d30.g)))),
# 
#                 depths = c(rep("d10", vcount(i_c_d10.g)),
#                            rep("d20", vcount(i_c_d20.g)),
#                            rep("d30", vcount(i_c_d30.g)),
#                            rep("d10", vcount(i_o_d10.g)),
#                            rep("d20", vcount(i_o_d20.g)),
#                            rep("d30", vcount(i_o_d30.g))),
# 
#                 netwk = c(rep("C_d10", vcount(i_c_d10.g)),
#                           rep("C_d20", vcount(i_c_d20.g)),
#                           rep("C_d30", vcount(i_c_d30.g)),
#                           rep("O_d10", vcount(i_o_d10.g)),
#                           rep("O_d20", vcount(i_o_d20.g)),
#                           rep("O_d30", vcount(i_o_d30.g))),
# 
#                 degrees_names = c(names(degree_c_d10.g),
#                                   names(degree_c_d20.g),
#                                   names(degree_c_d30.g),
#                                   names(degree_o_d10.g),
#                                   names(degree_o_d20.g),
#                                   names(degree_o_d30.g)),
# 
#                 degrees = c(    degree_c_d10.g,
#                                 degree_c_d20.g,
#                                 degree_c_d30.g,
#                                 degree_o_d10.g,
#                                 degree_o_d20.g,
#                                 degree_o_d30.g),
# 
# 
#                 betweenness_names = c(names(betw_c_d10.g),
#                                       names(betw_c_d20.g),
#                                       names(betw_c_d30.g),
#                                       names(betw_o_d10.g),
#                                       names(betw_o_d20.g),
#                                       names(betw_o_d30.g)),
# 
#                 betweenness = c(betw_c_d10.g,
#                                 betw_c_d20.g,
#                                 betw_c_d30.g,
#                                 betw_o_d10.g,
#                                 betw_o_d20.g,
#                                 betw_o_d30.g),
# 
#                 closeness_names = c(names(close_c_d10.g),
#                                     names(close_c_d20.g),
#                                     names(close_c_d30.g),
#                                     names(close_o_d10.g),
#                                     names(close_o_d20.g),
#                                     names(close_o_d30.g)),
# 
#                 closeness =       c(close_c_d10.g,
#                                     close_c_d20.g,
#                                     close_c_d30.g,
#                                     close_o_d10.g,
#                                     close_o_d20.g,
#                                     close_o_d30.g))
# 
# )
# 
# 
# 
# library(ggplot2)
# library(dplyr)
# library(tidytext)
# library(stringr)
# # degrees
# 
# tail(nwk_desc_df.g)
# 
# 
# svg(file.path(out_summary, "node_degrees_rho08_genus.svg"), width = 20, height = 15, pointsize = 20)
# nwk_desc_df.g %>%
#         dplyr::select(netwk, degrees_names, degrees) %>%
#         filter(degrees > 5) %>%
#         mutate(norder = reorder_within(degrees_names, degrees, netwk)) %>%
#         ggplot(aes(x = norder, y = degrees)) +
#         geom_segment(aes(x=norder, xend = norder, y=0, yend=degrees), color="grey") +
#         geom_point(size = 2, alpha=0.7, fill = "springgreen4", shape = 21) +
#         facet_wrap( ~ netwk, scales = "free") +
#         theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank()) +
#         scale_x_reordered() +
#         labs(y = "Nodes degrees", x = "Nodes", title="Node degrees (> 5) for co-occurrence networks") +
#         coord_flip()
# 
# dev.off()
# 
# 
# 
# 
# 
# svg(file.path(out_summary, "node_betweenness_rho08_genus.svg"), width = 20, height = 15, pointsize = 20)
# nwk_desc_df.g %>%
#         dplyr::select(netwk, betweenness_names, betweenness) %>%
#         filter(betweenness > 50) %>%
#         mutate(norder = reorder_within(betweenness_names, betweenness, netwk)) %>%
#         ggplot(aes(x = norder, y = betweenness)) +
#         geom_segment(aes(x=norder, xend = norder, y=0, yend=betweenness), color="grey") +
#         geom_point(size = 2, alpha=0.7, fill = "salmon", shape = 21) +
#         facet_wrap( ~ netwk, scales = "free") +
#         theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank()) +
#         scale_x_reordered() +
#         labs(y = "Nodes betweenness", x = "Nodes", title="Node betweenness (> 50) for co-occurrence networks") +
#         coord_flip()
# dev.off()
# 
# 
# 
# svg(file.path(out_summary, "node_closeness_rho08_genus.svg"), width = 20, height = 15, pointsize = 20)
# nwk_desc_df.g %>%
#         dplyr::select(netwk, closeness_names, closeness) %>%
#         filter(closeness > 1.949e-05) %>%
#         group_by(netwk) %>%
#         top_n(10) %>%
#         ungroup() %>%
#         mutate(norder = reorder_within(closeness_names, closeness, netwk)) %>%
#         ggplot(aes(x = norder, y = closeness)) +
#         geom_segment(aes(x=norder, xend = norder, y=0, yend=closeness), color="grey") +
#         geom_point(size = 2, alpha=0.7, fill = "dodgerblue", shape = 21) +
#         facet_wrap( ~ netwk, scales = "free") +
#         theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank()) +
#         scale_x_reordered() +
#         labs(y = "Nodes closeness", x = "Nodes", title="Node closeness (top 10) for co-occurrence networks") +
#         coord_flip()
# 
# dev.off()
# 
# 
# 
# cross.C.g = rbind(cross.C.d10.corr.g$df, cross.C.d20.corr.g$df, cross.C.d30.corr.g$df)
# cross.O.g = rbind(cross.O.d10.corr.g$df, cross.O.d20.corr.g$df, cross.O.d30.corr.g$df)
# 
# cross.C.g$from = as.character(cross.C.g$from)
# cross.C.g$to = as.character(cross.C.g$to)
# 
# cross.O.g$from = as.character(cross.O.g$from)
# cross.O.g$to = as.character(cross.O.g$to)
# 
# # how many bacterial interactions
# sum(startsWith(cross.C.g$from, "b") & startsWith(cross.C.g$to, "b")) # 479
# sum(startsWith(cross.O.g$from, "b") & startsWith(cross.O.g$to, "b")) # 43
# 
# bbC = cross.C.g[startsWith(cross.C.g$from, "b") & startsWith(cross.C.g$to, "b"), ]
# bbO = cross.O.g[startsWith(cross.O.g$from, "b") & startsWith(cross.O.g$to, "b"), ]
# 
# bbC %>%
#         arrange(-r) %>%
#         top_n(5)
# 
# 
# bbC %>%
#         filter(r < 0) %>%
#         arrange(-r) %>%
#         top_n(5)
# 
# 
# bbO %>%
#         filter(r == 1)
# 
# bbO %>%
#         arrange(-r) %>%
#         top_n(5)
# 
# 
# bbO %>%
#         filter(r < 0) %>%
#         arrange(r) %>%
#         top_n(5)
# 
# 
# bbC %>%
#         mutate(pair = paste(from, to, sep="+++")) %>%
#         group_by(pair) %>%
#         tally() %>%
#         mutate(perc = 100* n/479) %>%
#         arrange(-perc)
# 
# bbO %>%
#         mutate(pair = paste(from, to, sep="+++")) %>%
#         group_by(pair) %>%
#         tally() %>%
#         mutate(perc = 100* n/43) %>%
#         arrange(-perc)
# 
# # how many fungal interactions
# sum(startsWith(cross.C.g$from, "f") & startsWith(cross.C.g$to, "f")) # 85
# sum(startsWith(cross.O.g$from, "f") & startsWith(cross.O.g$to, "f")) # 10
# 
# ffC = cross.C.g[startsWith(cross.C.g$from, "f") & startsWith(cross.C.g$to, "f"), ]
# ffO = cross.O.g[startsWith(cross.O.g$from, "f") & startsWith(cross.O.g$to, "f"), ]
# 
# ffC %>%
#         arrange(-r) %>%
#         top_n(5)
# 
# ffC %>%
#         filter(r < 0) %>%
#         arrange(r) %>%
#         top_n(5)
# 
# 
# ffC %>%
#         mutate(pair = paste(from, to, sep="+++")) %>%
#         group_by(pair) %>%
#         tally() %>%
#         mutate(perc = 100* n/85) %>%
#         arrange(-perc)
# 
# ffO %>%
#         mutate(pair = paste(from, to, sep="+++")) %>%
#         group_by(pair) %>%
#         tally() %>%
#         mutate(perc = 100* n/10) %>%
#         arrange(-perc)
# 
# ffO %>%
#         arrange(-r)
# 
# 
# ffO %>%
#         filter(r < 0) %>%
#         arrange(-r)
# 
# # how many bacteria and fungal interactions
# sum(startsWith(cross.C.g$from, "f") & startsWith(cross.C.g$to, "b")) # 1
# sum(startsWith(cross.C.g$from, "b") & startsWith(cross.C.g$to, "f")) # 69
# 
# 
# sum(startsWith(cross.O.g$from, "f") & startsWith(cross.O.g$to, "b")) # 0
# sum(startsWith(cross.O.g$from, "b") & startsWith(cross.O.g$to, "f")) # 19
# 
# 
# bfC = cross.C.g[(startsWith(cross.C.g$from, "f") & startsWith(cross.C.g$to, "b")) | (startsWith(cross.C.g$from, "b") & startsWith(cross.C.g$to, "f")), ]
# bfO = cross.O.g[(startsWith(cross.O.g$from, "f") & startsWith(cross.O.g$to, "b")) | (startsWith(cross.O.g$from, "b") & startsWith(cross.O.g$to, "f")), ]
# 
# bfC %>%
#         arrange(-r) %>%
#         top_n(5)
# 
# bfC %>%
#         filter(r < 0) %>%
#         arrange(r) %>%
#         top_n(5)
# 
# 
# bfC %>%
#         mutate(pair = paste(from, to, sep="+++")) %>%
#         group_by(pair) %>%
#         tally() %>%
#         mutate(perc = 100* n/70) %>%
#         arrange(-perc)
# 
# 
# bfO %>%
#         mutate(pair = paste(from, to, sep="+++")) %>%
#         group_by(pair) %>%
#         tally() %>%
#         mutate(perc = 100* n/19) %>%
#         arrange(-perc)
# 
# 
# bfO %>%
#         arrange(-r)
# 
# bfC %>%
#         filter(r < 0) %>%
#         arrange(r) %>%
#         top_n(5)
