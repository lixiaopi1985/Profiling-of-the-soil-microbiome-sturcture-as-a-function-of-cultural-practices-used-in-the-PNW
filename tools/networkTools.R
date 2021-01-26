cordf = function(r.mat, p.mat){
        # construct rcorr output to a data frame
        lwer.tri = lower.tri(r.mat)
        df = data.frame(from=colnames(r.mat)[col(r.mat)[lwer.tri]], 
                        to=rownames(r.mat)[row(r.mat)[lwer.tri]],
                        r=r.mat[lwer.tri],
                        p = p.mat[lwer.tri])
        
        return(df)
}

normFun = function(mtrx, norm_method, norm_scale, taxaARErows = F){
        
        require(phyloseq)
        require(microbiome)
        
        if(taxaARErows){
                input_matrix = mtrx
        } else {
                input_matrix = t(mtrx)
        }
        
        
        phylo = phyloseq(otu_table(input_matrix, taxa_are_rows = T))
        norm_phylo = microbiome::transform(phylo, transform = norm_method, scale=norm_scale)
        normed.mat = as(t(otu_table(norm_phylo, taxa_are_rows = T)), "matrix")
        # taxa are columns and samples are rows
        
        return(normed.mat)
        
} 


permCor = function(mtrx, cor_method, norm_method, norm_scale, taxaARErows = F, nperm=1000, norm = T){
        # modified on Mehdi Layeghifard et al. Microbiome analysis
        require(progress)
        require(vegan)
        
        Ntaxa = dim(mtrx)[2]
        MB.mat = array(0, dim = c(Ntaxa, Ntaxa, nperm+1))
        
        cat("permCor: implementing permutation\n")
        
        # get permutation
        MBperm = vegan::permatswap(mtrx, times=nperm)
        
        # get normalization
        if(norm){
                cat("Normalizing matrix - method ", norm_method, "... \n")
                MB.norm = normFun(mtrx, norm_method = norm_method, norm_scale = norm_scale, taxaARErows = taxaARErows)
        } else {
                cat("No normalizing...\n")
                MB.norm = mtrx
        }
        
        MB.mat[,, 1] = as.matrix(cor(MB.norm, method = cor_method))
        
        
        pb = progress_bar$new(
                format = "permCor Calculating correlation :spin [:bar] :percent in :elapsed",
                clear = F,
                total=length(MBperm$perm))
        
        cat("Normalizing perm matrix - method ", norm_method, "... \n")
        for(i in 1:length(MBperm$perm)){
                pb$tick()
                permed.MB = MBperm$perm[[i]]
                if(norm){

                        permed.MB.norm = normFun(permed.MB, norm_method = norm_method, norm_scale = norm_scale, taxaARErows = taxaARErows)
                } else {
                        permed.MB.norm = permed.MB
                }
                MB.mat[,, i+1] = as.matrix(cor(permed.MB.norm, method = cor_method))
        }
        
        # evaluate p value for each pair of original cor matrix with the rest, get the total for each pair / number of perm
        
        
        pb2 = progress_bar$new(
                format = "permCor Calculating p-values :spin [:bar] :percent in :elapsed",
                clear = F,
                total=Ntaxa^2)
        pv = sapply(1:Ntaxa, function(i)sapply(1:Ntaxa, function(j){
                pb2$tick()
                sum(MB.mat[i,j,1]>MB.mat[i,j,2:nperm])
        }))
        
        pv = pv / nperm # p is a matrix
        cor_matrix = MB.mat[,,1]
        
        # add rownames and colnames
        
        rownames(pv) = colnames(mtrx)
        colnames(pv) = colnames(mtrx)
        rownames(cor_matrix) = colnames(mtrx)
        colnames(cor_matrix) = colnames(mtrx)
        
        return(list(r = cor_matrix, P=pv))
        
}



createCorDF = function(input_list, corTest, corrtype="spearman", nperm = 1000, taxlevel = "Phylum", taxaARErows = F, normalize = T, norm_method = "clr", norm_scale = 1, padj_method = "fdr", cor_cutoff = 0.6, padj_cutoff = 0.01){
        # input is a list containing matrix and tax mapping
        # normalize OTU matrix, col: OTUS, rows samples
        # calculate correlation
        # select by pvalue and rho
        # output adjacency matrix
        
        require(Hmisc)
        require(phyloseq)
        
        
        m_origin = input_list$m
        tax_df = input_list$taxa
        
        if(normalize){
                cat("Normalizing otu matrix ...\n")
                normed.mat = normFun(m_origin, norm_method = norm_method, norm_scale = norm_scale, taxaARErows = taxaARErows)
        } else {
                normed.mat = m_origin
        }
        
        
        # get correlation
        cat("Calculating correlation, method:",corrtype, " \n")
        if(corTest == "perm"){
                cat("Using permutation to get p value. Permutation times: ", nperm, " ...\n")
                cor_obj = permCor(m_origin, cor_method=corrtype, norm_method = norm_method, norm_scale = norm_scale, nperm = nperm, norm = normalize, taxaARErows = taxaARErows)
                
        } else if(corTest == "rcorr"){
                
                cat("Using rcorr in Hmisc package to get p value ...\n")
                cor_obj = Hmisc::rcorr(normed.mat, type=corrtype)
                
        } else if(corTest == "basecor"){
                
                cor.mtrx = cor(normed.mat, method = corrtype)
                cor.p = outer(colnames(normed.mat), colnames(normed.mat), Vectorize(function(i, j)cor.test(normed.mat[, i], normed.mat[, j], method = corrtype)$p.value))
                rownames(cor.p) = colnames(cor.mtrx)
                colnames(cor.p) = colnames(cor.mtrx)
                cor_obj = list(r=cor.mtrx, P=cor.p)
        }
        

        cor_df = cordf(cor_obj$r, cor_obj$P)
        cor_df$padj = p.adjust(cor_df$p, method = padj_method)
        
        # filter by rho and p
        
        cor_df.filtered = cor_df[which(abs(cor_df$r) > cor_cutoff), ]
        cor_df.filtered = cor_df.filtered[which(cor_df.filtered$padj < padj_cutoff), ]
        
        
        # cor_df.filtered: from, to, r, p, padj
        rownames(tax_df) = tax_df$OTUID
        
        
        tax_df[[taxlevel]] = paste0(substr(tax_df$OTUID, 1, 1), ":", tax_df[[taxlevel]])
        
        
        new_df = data.frame(from = tax_df[cor_df.filtered$from, taxlevel],
                            to = tax_df[cor_df.filtered$to, taxlevel],
                            r = cor_df.filtered$r)
        
        return(list(r=cor_obj$r, p=cor_obj$P, df = new_df, taxdf = tax_df))
        
        
}


cor2chord = function(input_list, showRelation = 0, plot_title="", titleCex=0.8,  chord_scale = F, chord_transparency = 0.4, colors=pals::alphabet2(), addLegend = T, legendPos = "bottomleft", legendCex = 0.6, outputLegend = F, combinedLegend = F, taxColorSchemeDF = NULL, directional = 0){
        
        require(RColorBrewer)
        require(circlize)
        
        if(combinedLegend){
                if(is.null(taxColorSchemeDF)){
                        stop("When use combinedLgend, you need to input a color scheme for all unique element, you can use taxColorScheme()")
                }
                new_df = input_list$df
                tax_df = input_list$taxdf
                tax_level = colnames(tax_df)[which(colnames(tax_df) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))]
                circos.clear()
                Leg = createChord(new_df=new_df, tax_df=tax_df, chord_scale = chord_scale, transparency = chord_transparency, colorPals = colors, taxlevel=tax_level, plot_title = plot_title, addLegend = addLegend, showRelation = showRelation, combineLegned = T, taxColorSchemeDF = taxColorSchemeDF, directional = directional)
                circos.clear()

                if(outputLegend){
                        return(Leg)
                }
                
        } else {
                new_df = input_list$df
                tax_df = input_list$taxdf
                tax_level = colnames(tax_df)[which(colnames(tax_df) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))]
                circos.clear()
                Leg = createChord(new_df=new_df, tax_df=tax_df, chord_scale = chord_scale, transparency = chord_transparency, colorPals = colors, taxlevel=tax_level, plot_title = plot_title, addLegend = addLegend, showRelation = showRelation)
                circos.clear()
                
                
                
                if(outputLegend){
                        return(Leg)
                }
        }

        return(Leg)
        
        

}


cor2igraph = function(corrlist){
        
        require(igraph)
        
        new_df = corrlist$df
        # load into igraph
        cor_graph = graph_from_data_frame(new_df, directed = F)


        # add node attributes
        #  time point, taxonomy, abundance?
# 
#         V(cor_graph)$time = ifelse(endsWith(V(cor_graph)$name, "t0"), "pre",
#                                    ifelse(endsWith(V(cor_graph)$name, "t1"), "week_1", ifelse(endsWith(V(cor_graph)$name, "t3"), "week_3", ifelse(endsWith(V(cor_graph)$name, "t6"), "week_6", "NA"))))
# 
#         V(cor_graph)$domain_org = ifelse(startsWith(V(cor_graph)$name, "b"), "Bacteria", "Fungi")
        
        return(cor_graph)
}


createChord = function(new_df, tax_df, taxlevel, colorPals, chord_scale = F, transparency = 0.4, plot_title = "", titleCex = 0.8,  addLegend=T, legendPos = "bottomleft", legendNcol = 1, legendCex = 0.6, showRelation = 0, combineLegned = F, taxColorSchemeDF = NULL, directional = 0){
        
        if(combineLegned){
                
                if(is.null(taxColorSchemeDF)){
                        stop("When use combinedLgend, you need to input a color scheme for all unique element, you can use taxColorScheme()")
                }
                
                cat("\nCombined legend colors used\n")
                taxColorMap = taxColorSchemeDF 
                
        } else {
                uniq_taxa = sort(unique(gsub("_t.?", "", tax_df[, taxlevel]))) # to make color match
                cols = colorRampPalette(colorPals)(length(uniq_taxa))
                taxColorMap = data.frame(taxaName = uniq_taxa, taxaCols = cols, stringsAsFactors = F)
        }

        
        
        
        sectors = union(new_df$from, new_df$to)
        
        sector_colors = sapply(sectors, function(x){
                x1 = gsub("_t.?", "", x)
                taxColorMap$taxaCols[taxColorMap$taxaName == x1]
        })
        names(sector_colors) = sectors
        
        
        # create groups by time
        time_groups = structure(gsub(".+_", "", sectors), names = sectors)
        time_groups_factors = factor(time_groups, levels = c("t0", "t1", "t3", "t6"))
        uniq_time_groups = sort(unique(time_groups_factors))
        
        if(showRelation == 1){
                cat("Showing positive links ...\n")
                linkCol = ifelse(new_df$r > 0, "tomato3", "#00000000") # this sets full transparency

        } else if(showRelation == -1){
                cat("Showing negative links ...\n")
                linkCol = ifelse(new_df$r < 0, "darkturquoise", "#00000000")
 
        } else if(showRelation == 0){
                cat("Showing both positive and negative links ...\n")
                linkCol = ifelse(new_df$r < 0, "darkturquoise", "tomato3")

        } else if(showRelation == 3){
                linkCol = sapply(new_df$from, function(x)sector_colors[x])
        }
        
        circos.initializeWithIdeogram(plotType = NULL)
        chordDiagram(new_df,
                     annotationTrack = "grid",
                     grid.col = sector_colors,
                     col = linkCol,
                     directional = directional,
                     preAllocateTracks = list(track.height = uh(5, "mm"),
                                              track.margin = c(uh(2, "mm"), 0)),
                     transparency = transparency,
                     group = time_groups, # change from time_groups_factors to time_groups
                     scale = chord_scale,
                     big.gap = 2)

        circos.track(track.index=1, ylim = c(0,1), panel.fun = function(x,y){
                
                circos.text(CELL_META$xcenter,
                            CELL_META$ycenter,
                            CELL_META$sector.numeric.index,
                            cex=0.7,
                            col="white",
                            niceFacing = T)
        }, bg.border = NA)
        
        
        
        highLight_label = c("pre", "week 1", "week 3", "week 6")
        highLight_cols = c("#86ABD0", "#86D0AB", "#86D086", "#ABD086")
        hightLightDF = cbind(highLight_label, highLight_cols)
        colnames(hightLightDF) = c("labels", "colors")
        rownames(hightLightDF) = c("t0", "t1", "t3", "t6")
        
        for(i in 1:length(uniq_time_groups)){
                
                highlight.sector(
                                 sectors[which(gsub(".*_", "", sectors) == uniq_time_groups[i])],
                                 CELL_META$xcenter,
                                 CELL_META$ycenter,
                                 track.index=1,
                                 text.col = "black",
                                 text = hightLightDF[i, "labels"],
                                 col = hightLightDF[i, "colors"],
                                 cex = 0.8,
                                 border = NA,
                                 facing = "bending",
                                 niceFacing=T)
        }
        
        title(plot_title, cex=titleCex)
        
        
        # color legend by sectors
        # taxColLegend = structure(taxColorMap$taxaCols, names=taxColorMap$taxaName)
        
        # remove duplicated because time was added on names(sector_colors)
        new_names = gsub("_t[0-9]$", "", names(sector_colors))
        names(sector_colors) = new_names
        sector_colors_nodup = sector_colors[!duplicated(sector_colors)]
        taxColLegend = structure(sector_colors_nodup, names=names(sector_colors_nodup))
        
        if(addLegend){
                legend(legendPos, pch=15, col=taxColLegend, legend = names(taxColLegend), cex=legendCex, box.lty = 0, ncol = legendNcol)
        } else {
                legendDF = data.frame(legendNames = names(taxColLegend), legendColors = taxColLegend)
                
                return(legendDF)
        }
        
        
}



taxColorScheme = function(corr_list, colorPals=pals::alphabet2()){
        # ... all dfs
        
        require(RColorBrewer)
        
        all_tax_df = do.call("rbind", lapply(corr_list, function(x)x$taxdf))
        all_new_df = do.call("rbind", lapply(corr_list, function(x)x$df))
        
        print(head(all_new_df))
                             
        tax_level = colnames(all_tax_df)[which(colnames(all_tax_df) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))]
        
        uniq_taxa = sort(unique(gsub("_t.?", "", all_tax_df[, tax_level]))) # to make color matching associated with taxa 
        cols = colorRampPalette(colorPals)(length(uniq_taxa))
        taxColorMap = data.frame(taxaName = uniq_taxa, taxaCols = cols, stringsAsFactors = F)
        
        sectors = union(all_new_df$from, all_new_df$to)
        
        sector_colors = sapply(sectors, function(x){
                x1 = gsub("_t.?", "", x)
                taxColorMap$taxaCols[taxColorMap$taxaName == x1]
        })
        names(sector_colors) = sectors
        
        # change the sector names to taxa then remove duplicates
        
        new_names = gsub("_t[0-9]$", "", names(sector_colors))
        names(sector_colors) = new_names
        sector_colors_nodup = sector_colors[!duplicated(sector_colors)]
        
        legendDF = data.frame(legendNames = names(sector_colors_nodup), legendColors = sector_colors_nodup, stringsAsFactors = F)
        
        return(list(all_taxa_colors = taxColorMap, all_sector_colors = legendDF))
}



quickCorFilter = function(mtx_obj, padj_method = "BH", cor_cutoff = 0.6, padj_cutoff = 0.01){
        
        cor_obj = Hmisc::rcorr(mtx_obj, type=c("spearman"))
        
        cor_df = cordf(cor_obj$r, cor_obj$P)
        cor_df$padj = p.adjust(cor_df$p, method = padj_method)
        
        # filter by rho and p
        cor_df.filtered = cor_df[which(abs(cor_df$r) >= cor_cutoff), ]
        cor_df.filtered = cor_df.filtered[which(cor_df.filtered$padj < padj_cutoff), ]
        
        return(cor_df.filtered)
}



quickTaxaCombine = function(cor_df.filtered, tax_df, tax_level){
        
        # cor_df.filtered: from, to, r, p, padj
        rownames(tax_df) = tax_df$OTUID
        
        tax_df[[taxlevel]] = paste0(substr(tax_df$OTUID, 1, 1), ":", tax_df[[taxlevel]])
        
        new_df = data.frame(from = tax_df[cor_df.filtered$from, taxlevel],
                            to = tax_df[cor_df.filtered$to, taxlevel],
                            r = cor_df.filtered$r)
        
        return(new_df)
}

ObjSave <- function(..., folder) {
        objects <- list(...)
        object_names <- sapply(substitute(list(...))[-1], deparse)
        sapply(1:length(objects), function(i) {
                filename = paste0(folder, "/", object_names[i], ".rds")
                saveRDS(objects[i], filename)
        })
}
