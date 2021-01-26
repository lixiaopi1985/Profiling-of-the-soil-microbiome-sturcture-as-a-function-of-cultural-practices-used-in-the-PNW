showTopN = function(phylo, n, taxlevel, writeout = F, writeoutpath = NULL, writeoutname = NULL,  plotOUT=F, ...){
        require(microbiome)
        require(microbiomeutilities)
        require(dplyr)
        
        aggTaxa = microbiome::aggregate_taxa(phylo, taxlevel)
        cat("\nTranform to relative abundance\n")
        aggTaxa_ldf = microbiomeutilities::phy_to_ldf(aggTaxa, transform.counts = "compositional")
        
        sym_taxlevel = sym(taxlevel)
        
        topN = aggTaxa_ldf %>% 
                select(!!sym_taxlevel, Abundance) %>%
                group_by(!!sym_taxlevel) %>%
                summarise(avg_rel_abd = mean(Abundance)) %>%
                arrange(-avg_rel_abd) %>%
                top_n(n)
        
        if(writeout){
                
                if(is.null(writeoutpath)){
                        getpath = getwd()
                } else {
                        getpath = writeoutpath
                }
                
                if(is.null(writeoutname)){
                        getname = paste0("Top_", n, "_", taxlevel, ".csv")
                } else {
                        getname = writeoutname
                }
                
                write.csv(topN, file = file.path(getpath, getname), quote = F)
                 
        }
        
        
        if(plotOUT){
                aggTaxa_ldf = aggTaxa_ldf %>%
                        filter(!!sym_taxlevel != "Other")
                
                ggstripchart(aggTaxa_ldf, ...)
                             # x = "period", xlab = "Time", y = "Abundance", ylab = "Relative abundance (%)", color = "period", facet.by = c("soil_depth", "Phylum"), add=c("boxplot", "mean_sd"), add.params = list(fill="grey90"), ggtheme = theme_bw())  + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        }
        
        return(topN)
        
}