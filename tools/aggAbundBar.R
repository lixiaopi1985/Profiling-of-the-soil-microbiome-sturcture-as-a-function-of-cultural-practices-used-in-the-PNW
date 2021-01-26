aggAbund_barplot = function(phylo, taxlevel, colors, x = "Sample", y = "Abundance", facet_ = "cul_type", xLab=NULL, yLab="Relative abundance (%)", perc = T, abund = 1/100, NArm = T, legpos = "bottom", ...){
        
        
        require(ggplot2)
        require(pals)
        require(tidyverse)
        require(rlang)
        
        colPatte = colorRampPalette(colors)
        
        
        if(perc){
                # agglomerate taxa
                glom = tax_glom(phylo, taxrank = taxlevel, NArm = NArm) %>%
                        transform_sample_counts(function(x)100*x/sum(x))
                abund = 100*abund
        } else {
                glom = tax_glom(phylo, taxrank = taxlevel, NArm = NArm) %>%
                        transform_sample_counts(function(x)x/sum(x))
        }
        

        dat = psmelt(glom) 
        
        # group by tax level than calculate mean
        dat[[taxlevel]] = as.character(dat[[taxlevel]])
        
        # filter
        if(perc){
                dat[[taxlevel]] = ifelse(dat$Abundance < abund, paste("<", abund, "%"),  dat[[taxlevel]])
        } else {
                dat[[taxlevel]] = ifelse(dat$Abundance < abund, paste("<", abund),  dat[[taxlevel]])
        }
        
        
        
        facet_arg = paste0("~", facet_)
        
        n_tax = length(unique(dat[[taxlevel]]))
        
        outplot = 
                dat %>%
                ggplot(aes_string(x=x, y=y, fill=taxlevel)) +
                geom_bar(stat = "identity", position = "fill") +
                facet_grid(facet_arg, scales = "free_x", space="free_x", ...) +
                scale_fill_manual(values = colPatte(n_tax), na.value = "grey45") +
                labs(x=xLab, y=yLab) +
                theme_bw()+
                theme(axis.text.x = element_text(angle = 90, hjust=1), legend.position = legpos)

        print(outplot)
        return(dat)

}
