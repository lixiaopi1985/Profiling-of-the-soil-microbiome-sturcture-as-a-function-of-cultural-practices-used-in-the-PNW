IndicTaxaSumm = function(indic_attr, grp, taxlevel = "Phylum"){
        
        require(tidyverse)
        require(rlang)
        
        # what groups in the data
        cat("\nSummary of indic result:\n\n")
        cat("\n++++++++++++++++++++++++++++++\n")
        cat("Unique groups:\n")
        uniq_grp = unique(indic_attr$indicgroup)
        cat(uniq_grp)
        
        # how many ASVs in the grp
        cat(paste("\n\nHow many ASVs in the group:", grp), ":\n")
        ASVsum = sum(indic_attr$indicgroup == grp)
        cat(ASVsum)
        
        cat(paste("\n\nHow many (Unique) ASVs in the group:", grp), ":\n")
        ASV_unique = length(unique(indic_attr$node[indic_attr$indicgroup == grp]))
        cat(ASV_unique)
        
        
        # how many different taxa in the grp
        cat(paste("\n\nHow many unique", taxlevel, "in the", grp), ":\n")
        taxaCount = length(unique(indic_attr[[taxlevel]][(indic_attr$indicgroup == grp) & (!is.na(indic_attr[[taxlevel]]))]))
        cat(taxaCount)
        
        # percentage
        cat(paste("\n\nPercentage of the", taxlevel), ":\n\n")
        percTaxa = indic_attr %>%
                filter(indicgroup == grp & !is.na(!!sym(taxlevel))) %>%
                select(!!sym(taxlevel)) %>%
                group_by(!!sym(taxlevel)) %>%
                tally() %>%
                mutate(perc = 100*round(n/sum(n), 3)) %>%
                arrange(-perc)
        print(percTaxa)
        cat("\n+++++++++++++++++++++++++++++++++++++\n")
        
        return(list(uniq_group = uniq_grp, sumASVs = ASVsum, taxaCounts = taxaCount, taxaPerc = percTaxa))
}



indic_pieChart = function(splist, legend="Phylum", barwidth = 1, min_prop = 1, x_baselabel = 1.6, n=0.01, inc = 0.05, title=NULL){
        require(ggplot2)
        
        splist_2 = droplevels(splist)
        df = data.frame(sort(table(splist_2), decreasing = T))
        # get frequency
        colnames(df) = c("taxa", "counts")
        
        total_count = sum(df$counts)
        
        df$prop = round(100*(df$counts/total_count),1)
        
        sort_prop = sort(df$prop)
        df$yloc = sort(cumsum(sort_prop) - sort_prop/2, decreasing = T)
        
        
        xlabel = c()
        for(i in 1:length(df$prop)){
                p = df$prop[i]
                if(p < min_prop){
                        xlabel = c(xlabel, x_baselabel+n)
                        n = n+inc
                } else {
                        xlabel = c(xlabel, x_baselabel)
                }
                
        }
        
        
        df$xlabel = xlabel
        
        # segment position
        xseg = 1 + barwidth/2
        xseg_end = ifelse(df$prop < min_prop, xseg + (df$xlabel - xseg)/2, xseg)
        
        df$xseg = xseg
        df$xseg_end = xseg_end
        

        pie_graph = ggplot(df, aes(x="", y=prop, fill = taxa)) +
                        geom_bar(width=barwidth, stat="identity", color = "white", alpha=0.6) +
                        geom_text(aes(x=xlabel, y = yloc, label=paste(prop, "%"))) +
                        geom_segment(aes(x=xseg, xend=xseg_end, yend=yloc), position = position_stack(vjust=0.5)) +
                        scale_fill_manual(values = glasbey(nrow(df))) +

                        theme_void() +
                        labs(title=title, fill=legend) +
                        coord_polar("y", start=0)
        
        return(pie_graph)
        
}