addTime = function(x, filepath, patterns="[0-9]wk" ){
        require(stringr)
        df = read.csv(file.path(filepath, x))
        timelabel = str_extract(x, patterns)
        df$timestamp = timelabel
        return(df)
}

transformMatrix = function(metadata, colselected, otucol){
        # create matrix
        clean_meta = metadata %>%
                select(X, timestamp, log2FoldChange) %>%
                spread(timestamp, log2FoldChange)
        
        # NA sets to 0
        clean_meta[is.na(clean_meta)] = 0
        # make matrix  
        clean_meta_mtrx = as.matrix(clean_meta[, colselected])
        # chane row names
        rownames(clean_meta_mtrx) = clean_meta[[otucol]]
        
        return(clean_meta_mtrx)
}

plotggtree = function(phylo_object, clean_meta_mtrx, scaleColorRange = T, rngMtrix = NULL, low = "steelblue", high = "red", mid = "grey94", midPoint=0, breaks_by = 2, newscalename = "Log2 fold changes", colFUN = pals::glasbey(), changeColnames = NULL, addColnames = T, colPos = "bottom", hmap_offset = 1, hmap_width = 0.4, tip_align = T, xlimTree = 10, ...){
        require(ggtree)
        require(ggstance)
        

        if(!is.null(changeColnames)){
                colnames(clean_meta_mtrx) = changeColnames
        }

        
        n = length(unique(tax_table(phylo_object)[, "Phylum"]))
        
        # plot tree
        if(is.null(xlimTree)) {
                gtree = ggtree(phylo_object, ladderize = F)+ geom_tiplab(aes(label=Genus), align = tip_align, offset = 0.1) + geom_tippoint(aes(color=Phylum), size = 3)
        } else {
                gtree = ggtree(phylo_object, ladderize = F)+ geom_tiplab(aes(label=Genus), align = tip_align, offset = 0.1) + geom_tippoint(aes(color=Phylum), size = 3) + xlim_tree(xlimTree)
        }
        
        # plot heatmap
        if(scaleColorRange){
                if(is.null(rngMtrix)){
                        stop("Need to supply range matrices")
                }
                
                maxLimit = ceiling(rngMtrix[2])
                minLimit = floor(rngMtrix[1])
                newLimit = c(minLimit, maxLimit)
                breakrange = seq(minLimit, maxLimit, by = breaks_by)
                
                ghmap = gheatmap(gtree, clean_meta_mtrx, offset = hmap_offset, width = hmap_width, legend_title = "Log2 Fold Change", colnames = addColnames, colnames_position = colPos, ...) +
                        scale_color_manual(values = colFUN) + 
                        scale_fill_gradient2(low = low, mid = mid, high = high, midpoint = midPoint, breaks = breakrange, limits = newLimit, name = newscalename) +
                        theme(legend.position = "left")
                
                
        } else {
                ghmap = gheatmap(gtree, clean_meta_mtrx, offset = hmap_offset, width = hmap_width, legend_title = "Log2 Fold Change", colnames = addColnames, colnames_position = colPos, ...) +
                scale_color_manual(values = colFUN) + 
                theme(legend.position = "left")
        }
        
        return(ghmap)
}
