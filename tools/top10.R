top_taxa = function(phylo, taxlevel, groupby = NULL, perc = T, aggFun = mean, topn=10){
        
        require(rlang)
        
        # agglomerate taxa
        
        
        if(perc){
                phylo.norm = transform_sample_counts(phylo, function(x)100*x/sum(x))
        } else {
                phylo.norm = transform_sample_counts(phylo, function(x)x/sum(x))
        }
        
        glom = tax_glom(phylo.norm, taxrank = taxlevel)

        # psmelt
        dat = psmelt(glom) 
        
        # group by tax level than calculate mean
        dat[[taxlevel]] = as.character(dat[[taxlevel]])
        
        if( (!is.null(groupby)) && (groupby %in% colnames(dat))){
                dat_top = dat %>%
                        group_by(!!sym(groupby), !!sym(taxlevel)) %>%
                        summarise(mean_abd = aggFun(Abundance)) %>%
                        top_n(n=topn, wt=mean_abd) %>%
                        arrange(!!sym(groupby), -mean_abd) %>%
                        print(n=Inf)
        } else {
                dat_top = dat %>%
                        group_by(!!sym(taxlevel)) %>%
                        summarise(mean_abd = aggFun(Abundance)) %>%
                        arrange(desc(mean_abd)) %>%
                        print(n=topn)
        }


        return(dat)
}
