makeScatterPlot = function(sigtab, contrasts, sec_tax, fst_tax = "Phylum", title = NULL, point_size=4, show_n=NULL, log2FC = NULL, legendPos = "right"){
        # sort table
        baselevel = contrasts[3]
        numerator = contrasts[2]
        
        
        sigtab2 = subset(sigtab, !is.na(sec_tax))
        
        x = tapply(sigtab2$log2FoldChange, sigtab2[[fst_tax]], function(x)max(x))
        x = sort(x, T)
        
        sigtab2[fst_tax] = factor(as.character(sigtab2[[fst_tax]]), levels=names(x))
        
        x = tapply(sigtab2$log2FoldChange, sigtab2[[sec_tax]], function(x)max(x))
        x = sort(x, T)
        
        sigtab2[sec_tax] = factor(as.character(sigtab2[[sec_tax]]), levels=names(x))
        
        if(is.null(show_n)){
                sigtab_n = sigtab2
        } else {
                
                if(nrow(sigtab2) < show_n){
                        stop("show_n is larger than the number detected")
                }
                
                sigtab_n = sigtab2[1:show_n, ]
        }
        
        

        g = ggplot(sigtab_n, aes_string(x= sec_tax, y="log2FoldChange", color=fst_tax)) +
                geom_jitter(size=point_size, width = 0, height = 0.1, alpha=0.4) +
                geom_hline(aes(yintercept = 0), color = "grey") + 
                theme_classic() + 
                ggtitle(title, subtitle = paste("Top significant ASVs with |lod2FC| >=", log2FC, ":", numerator, "vs", baselevel)) + 
                theme(plot.subtitle = element_text(size=8, color="grey24", face = "italic"),
                      axis.text.x = element_text(angle=-90, hjust=0, vjust=0), 
                      legend.position = legendPos,
                      text = element_text(size=8))
        
        
        
        return(g)
        
}


makeBarplot = function(sigtab, contrasts, sec_tax, fst_tax = "Phylum", show_n=5, title=NULL, log2FC = NULL, legendPos = "right", textsize = 15){
        # show_n show top and bottom
        

        baselevel = contrasts[3]
        numerator = contrasts[2]

        # x_up = sigtab[order(sigtab$log2FoldChange, decreasing = T), ]
        
        x_up = sigtab %>% tibble::rownames_to_column("ASV_IDs") %>%  dplyr::filter(log2FoldChange > 0) %>% dplyr::arrange(-log2FoldChange) %>% mutate(change = "Enriched")
        
        # x_down = sigtab[order(sigtab$log2FoldChange, decreasing = F), ]
        x_down = sigtab %>% tibble::rownames_to_column("ASV_IDs") %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::arrange(log2FoldChange) %>% mutate(change = "Depleted")
        
        
        if(show_n > nrow(x_up) ){
                xup_show = x_up[1:nrow(x_up), ]
        } else if (is.null(show_n)){
                xup_show = x_up[1:nrow(x_up), ]
        } else{
                xup_show = x_up[1:show_n, ]
        }
        
        if(show_n > nrow(x_down)){
                xdown_show = x_down[1:nrow(x_down), ]
        } else if(is.null(show_n)) {
                xdown_show = x_down[1:nrow(x_down), ]
        } else {
                xdown_show = x_down[1:show_n, ]
        }
        x = rbind(xup_show, xdown_show)
        
        x = as.data.frame(x)
        
        ylabels = paste(x[["ASV_IDs"]], x[[fst_tax]], x[[sec_tax]], sep="_")
        
        x$yDisp = as.character(ylabels)
        
        x$yDisp = gsub("NA", 'unassigned', x$yDisp)
        

        g = ggplot(x, aes(x = log2FoldChange, y = reorder(yDisp, log2FoldChange), fill=change)) +
                geom_bar(stat = "identity", orientation = "y", alpha=0.7) +
                scale_fill_manual(values = c("cyan3", "salmon")) +
                theme_classic() +
                ggtitle(title, subtitle = paste("Top", show_n, "significant ASVs with |lod2FC| >=", log2FC, ":", numerator, "vs", baselevel, "(baselevel)")) + 
                theme(plot.subtitle = element_text(size=textsize, color="grey24", face = "italic"), legend.position = legendPos, plot.margin = unit(c(1,1,1,1), "mm"), axis.text.x = element_text(size=textsize), axis.text.y = element_text(size=textsize)) +
                labs(y = paste("ASV", fst_tax, sec_tax, sep = "_"), fill = "log2 fold change") 
        
        
        return(g)
        
}



makeVolcano = function(sigtab, contrasts, taxTable, x = "log2FoldChange", y = "padj", vlab=NULL, title=NULL, subtitle = T,subtitleLabSize = 8, ...){
        
        require(EnhancedVolcano)
        
        baselevel = contrasts[3]
        numerator = contrasts[2]
        
     
        if(subtitle){
                subt = paste(numerator, "vs", baselevel)       
        } else {
                subt = NULL
        }

        # label
        plab = ""
        

        if(!is.null(vlab)){
                if(toupper(vlab) %in% c("ASV", "OTU")){
                        plab = rownames(sigtab)
                } else if (toupper(vlab) %in% toupper(colnames(taxTable))) {
                        plab = sigtab[[vlab]]
                }
        } else {
                plab = rownames(sigtab)
        }
        

        
        g = EnhancedVolcano(sigtab, 
                            lab = plab, 
                            x = x, 
                            y = y, 
                            title=title, 
                            subtitle = subt, 
                            subtitleLabSize = subtitleLabSize, 
                            ...)
        
        return(g)
}


plotDeseq2 = function(deseq_res, phyloseq_obj, contrasts, lfc_Shrink = T, alpha = 0.05, tax_level = NULL,  log2FC = 2, show_n = NULL, plottype="scatter", title=NULL, outputpath =NULL, summar = F, summout = NULL, legendPos = "right", vlab = NULL, subtitle = T,subtitleLabSize = 8, ...){
        # plottype: scatter
        # plottype: hbar
        require(phyloseq)
        require(DESeq2)
        require(ggplot2)
        require(dplyr)
        
        
        # reorder
        if(class(deseq_res) != "DESeqDataSet"){
                stop("deseq_res must be a deseq object")
        }
        
        if(class(phyloseq_obj) != "phyloseq"){
                stop("phyloseq_obj must be a phyloseq object")
        }
        
        res0 = deseq_res
        

        cat(paste("Contrasts: ", resultsNames(res0)), "\n\n")
        res = results(res0, contrast = contrasts)
        
        if(lfc_Shrink){
                cat("Running LFC shrinkage to get accurate log2 fold change...\n")
                res = lfcShrink(res0, contrast = contrasts, res=res)
        }
        

        
        phylo = phyloseq_obj
        taxTable = tax_table(phylo)
        
        baselevel = contrasts[3]
        numerator = contrasts[2]

        
        res = res[order(res$padj, na.last=NA),]
        sigtab = res[(res$padj < alpha),]
        
        # summary
        
        cat("\nBasic summary:\n\n")
        cat(paste("Parameters: alpha for padj", alpha, "log2FoldChange cutoff:",log2FC), "\n")
        cat(paste("Alpha level: ", alpha, ". Total number of significant ASV detected: ", dim(sigtab)[1]), "\n")
        cat(paste("Number of ASVs with log2 FC >= 1", sum(sigtab$log2FoldChange >=1), ". Number of ASVs with log2 FC >=2", sum(sigtab$log2FoldChange >= 2)), "\n")
        cat(paste("Number of ASVs with log2 FC <= -1", sum(sigtab$log2FoldChange <= -1), ". Number of ASVs with log2 FC <= -2", sum(sigtab$log2FoldChange <= -2)), "\n")
        

        sigtab = cbind(as(sigtab, "data.frame"), as(taxTable[rownames(sigtab), ], "matrix"))
        sigtab_foldchange = sigtab[(abs(sigtab$log2FoldChange) >= log2FC), ]
        sigtab_foldchange$cultypes = gsub("_.*$", "", numerator)
        sigtab_foldchange$contrasts = paste0(gsub("^.*_", "", numerator), "_vs_", gsub("^.*_", "", baselevel))
        
        if(summar){
                # summarise in taxa
                cat("\n++++++++++++++++++++++++++++++++")
                cat("\n\nTaxonomy summary\n\n")
                cat("\n1. ---- log2FC >= 1 at phylum level:\n")
                
                
                
                summ1 = sigtab %>% 
                        dplyr::filter(abs(log2FoldChange) >= 1) %>%
                        dplyr::select(Phylum) %>%
                        dplyr::group_by(Phylum) %>%
                        dplyr::tally() %>%
                        dplyr::mutate(tax_perc = 100*n/sum(n)) %>%
                        dplyr::arrange(-tax_perc)
                
                print(summ1)
                
                cat("\n2. ---- log2FC >= 1 at genus level:\n")
                
                summ2 = sigtab %>% 
                        dplyr::filter(abs(log2FoldChange) >= 1) %>%
                        dplyr::select(Genus) %>%
                        dplyr::filter(!is.na(Genus)) %>%
                        dplyr::group_by(Genus) %>%
                        dplyr::tally() %>%
                        dplyr::mutate(tax_perc = 100*n/sum(n)) %>%
                        dplyr::arrange(-tax_perc)
                
                print(summ2)
                
                cat("\n3. ---- log2FC >= 2 at phylum level:\n")
                
                summ3 = sigtab %>% 
                        dplyr::filter(abs(log2FoldChange) >= 2) %>%
                        dplyr::select(Phylum) %>%
                        dplyr::group_by(Phylum) %>%
                        dplyr::tally() %>%
                        dplyr::mutate(tax_perc = 100*n/sum(n)) %>%
                        dplyr::arrange(-tax_perc)
                
                print(summ3)
                
                cat("\n4. ---- log2FC >= 2 at genus level:\n")
                
                summ4 = sigtab %>% 
                        dplyr::filter(abs(log2FoldChange) >= 2) %>%
                        dplyr::select(Genus) %>%
                        dplyr::filter(!is.na(Genus)) %>%
                        dplyr::group_by(Genus) %>%
                        dplyr::tally() %>%
                        dplyr::mutate(tax_perc = 100*n/sum(n)) %>%
                        dplyr::arrange(-tax_perc)
                
                print(summ4)
        }

        
        
        if(!is.null(summout)){
                
                colnames(summ1) = c("taxa", "counts", "percentage(%)")
                colnames(summ2) = c("taxa", "counts", "percentage(%)")
                summ_all = rbind(summ1, summ2)
                write.csv(summ_all, file = file.path(summout, paste0("DE_summary_", numerator, "_vs_", baselevel, ".csv")), quote = F, row.names = F)
        }
        
        cat("\n++++++++++++++++++++++++++++++++++++\n\n")
                
        # output file
        if(!is.null(outputpath)){
                f1 = file.path(outputpath, paste0(numerator, "_vs_", baselevel, "_", "all_significantOTU_level", alpha, ".csv"))
                write.csv(sigtab, f1)
                
                f2 = file.path(outputpath, paste0(numerator, "_vs_", baselevel, "_", "significantOTU_withfoldChange_lgt", log2FC, ".csv"))
                write.csv(sigtab_foldchange, f2)
        }
        
        # output graph
        
        if(!is.null(tax_level)){
                
                if(length(tax_level) == 1){
                        
                        if(!tax_level %in% rank_names(phylo)){
                                stop("tax_level is not in the rank names")
                        }
                        cat("tax_level is one, default phylum used as the first taxa\n")
                        
                        
                        if(plottype=="scatter"){
                                cat("\n\nMaking scatter plot...\n\n")
                                g = makeScatterPlot(sigtab_foldchange, contrasts = contrasts, sec_tax = tax_level, title = title,show_n = show_n, log2FC = log2FC, legendPos = legendPos)    
                        } else if(plottype == "hbar"){
                                cat("\n\nMaking horizontal bar plot...\n\n")
                                g = makeBarplot(sigtab_foldchange, contrasts = contrasts, sec_tax = tax_level, show_n = show_n, title = title, log2FC = log2FC, legendPos = legendPos)
                        } else if(plottype == "volcano"){
                                cat("\n\nMaking volcano plot...\n\n")
                                g = makeVolcano(sigtab, contrasts = contrasts, taxTable = taxTable, vlab=vlab, title = title, subtitle = subtitle, subtitleLabSize = subtitleLabSize, ...)
                        }
                        
                        
                        
                        
                }
                
                if(length(tax_level) == 2){
                        if(sum(tax_level %in% rank_names(phylo)) != 2){
                                stop("tax_level must be length of 2 and be of the rank names in phyloseq")
                        }
                        
                        if(plottype == "scatter"){
                                g = makeScatterPlot(sigtab_foldchange, contrasts = contrasts, tax_level[1], tax_level[2], title = title, show_n = show_n, log2FC = log2FC)
                        } else if(plottype == "hbar") {
                                g = makeBarplot(sigtab_foldchange, contrasts = contrasts, sec_tax = tax_level, show_n = show_n, title = title, log2FC = log2FC)
                        } else if(plottype == "volcano"){
                                g = makeVolcano(sigtab, contrasts = contrasts, taxTable = taxTable, vlab = vlab, title = title, subtitle = subtitle, subtitleLabSize = subtitleLabSize, ...)
                        }
                }
                

        }
        
        return(list(graph = g, df = sigtab_foldchange))
        
}


summ_deseq2 = function(ddsout, taxLevel = "Phylum"){
        
        require(tidyverse)
        require(rlang)
        
        df = ddsout$df
        
        # how many significant
        cat("Total ASVs:", nrow(df), "\n")
        
        # how many up
        cat("Enriched: ", nrow(df[df$log2FoldChange > 0,]), "\n")
        
        # how many down
        cat("Depleted: ", nrow(df[df$log2FoldChange < 0,]), "\n")
        
        # summarize taxonomy level
        cat("-------------------------------------\n")
        cat("Taxonomy composition in all significant:\n")
        
        
        taxlevel = sym(taxLevel)
        
        
        df.0 =df %>%
                dplyr::group_by(!!taxlevel) %>%
                tally() %>%
                mutate(total_counts = sum(n), comp_perc = 100*(n/total_counts)) %>%
                arrange(-comp_perc)
        print(df.0)
        
        cat("\nTaxonomy in enriched:\n")
        df.up = df %>%
                dplyr::filter(log2FoldChange > 0) %>%
                dplyr::group_by(!!taxlevel) %>%
                tally() %>%
                mutate(total_counts = sum(n), comp_perc = 100*(n/total_counts)) %>%
                arrange(-comp_perc)
        print(df.up)
        
        cat("\nTaxonomy in depleted:\n")
        df.down = df %>%
                dplyr::filter(log2FoldChange < 0) %>%
                dplyr::group_by(!!taxlevel) %>%
                tally() %>%
                mutate(total_counts = sum(n), comp_perc = 100*(n/total_counts)) %>%
                arrange(-comp_perc)
        print(df.down)
        
        cat("------------------------------------------------------------\n")
        
        return(list(df_all = df.0, df_up = df.up, df_down = df.down))
}