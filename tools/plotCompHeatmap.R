searchUP = function(taxaTable, taxLevel){
        if(!taxLevel %in% colnames(taxaTable)){
                stop("taxa level not found in the ranking names in taxa table")
        }
        
        taxaTable[taxaTable == "unassigned"] = NA
        
        rank_pos = which(colnames(taxaTable) == taxLevel)
        
        out = apply(taxaTable[, 1:rank_pos], 1, function(x){
                max(which(complete.cases(x)))
        })
        
        return(out)
}

#' @param taxaTable tax_table object
#' @param taxLevel tax rank selected
#' @param method "moveup" or "unassigned

reassignNARank = function(taxaTable, taxLevel, method = "moveup"){
        # search up tax lineage to find higher level tax rank
        # return a vector
        
        # check column names


        if(!taxLevel %in% colnames(taxaTable)){
                stop("taxa level not found in the ranking names in taxa table")
        }
        
        taxaTable = data.frame(taxaTable, stringsAsFactors = F)
        

        out = searchUP(taxaTable, taxLevel)

        
        col_names = colnames(taxaTable)
        

        
        selected_rank = as.vector(taxaTable[, taxLevel])

        sp = vector()

        if(method=="unassign"){
                sp = sapply(selected_rank, function(x){
                        ifelse(is.na(x), "unassign", x)
                })
        } else if(method == "moveup"){
                
                
                for(i in 1:nrow(taxaTable)){
                        rowN = taxaTable[i,] # row
                        names_row = rownames(rowN) # get OTU index from data
                        rankLevel = out[names_row] # pick the OTU index
                        pickSp = ifelse(col_names[rankLevel] != taxLevel, paste0("[", substr(col_names[rankLevel][1],1,1), ":", rowN[rankLevel], "] unassigned"), rowN[rankLevel]) # if the corresponding tax rank not equals to the selected rank, we would know it is NA
                        sp = c(sp, unlist(pickSp))
                }
        }

        return(sp)
}



#' @param x a vector with whose length we can generate a list for each level
#' @return color lists

generateCols = function(x, meta, colCategory="qual",rng.seed=123){
        
        n_len = length(x)
        colorSets = brewer.pal.info
        
        whichset = rownames(colorSets[colorSets$category == "qual",])
        set.seed(rng.seed)
        sampleColor = sample(whichset, n_len)
        
        assignColor = list()
        
        for(i in 1:n_len){
                
                grp = unique(meta[[x[i]]])
                seqLong = seq_along(grp)
                names(seqLong) = grp
                
                colorSet = brewer.pal(length(grp),sampleColor[i])
                
                grp_colr = sapply(seqLong, function(x)colorSet[x])
                
                assignColor[[x[i]]] = grp_colr 
        }
        
        return(assignColor)
        
}




#' @description plot heatmap of compostional data
#' @param obj phyloseq obj
#' @param show_n int, top n to show 
#' @param norm logical, normalized the data. e.g turn counts to proportion 
#' @param use_hclust character, use vegdist distance method, e.g. "bray
#' @param hclust_method set cluster method for hclust
#' @param tranform how do you normalize the data, see microbiome
#' @param target how do you normalize the data, see micobiome
#' @param shift see microbiome
#' @param scale see microbiome
#' @param assignMethod how to demonstrate row labels(spieces) when encounter NA, moveup, search tax ranking up a rung.
#' @param fun how you choose top show_n, default mean
#' @param hp_transform if you want to transform your matrix, for plotting
#' @param percent if proportion data needs to convert to percentage
#' @param taxa_rank if choose anything other than OTU, e.g. Phylum, it will collapse taxa
#' @param userMeta if user provides metadata
#' @param metainput if userMeta == T, then this can not be NULL and has to be a dataframe
#' @param anno_list add annotation bar, input a list specify position as name list(column="period", row="...")
#' @param anno_names a list to change the names of annotation bars
#' @param dend_col show dendrogram on column names
#' @param dend_row show dendrogram on row names
#' @param hp_colr customize heatmpa color named list
#' @param anno_colr customize annotation bar color
#' @param rowsplit split rows by what factor
#' @param clusterRowSlice Logical, cluster sliced rows
#' @param colsplit split columns by what factor
#' @param clusterColSlice Logical, cluster sliced columns
#' @param change_colnames Character, "seq" - change the labels of column names to S#
#' @param ... extra arguments parsed to ComplexHeatmap object 

plotCompHeatMap = function(obj, 
                           show_n=10, 
                           norm=T, 
                           use_hclust = NULL, 
                           hclust_pos = NULL,
                           hclust_method="complete", 
                           transform="compositional", 
                           target="OTU", 
                           shift=0, 
                           scale=1,
                           assignMethod = "moveup",
                           fun=mean, 
                           hp_transform = NULL, 
                           hp_colr_alpha = 1,
                           percent = T, 
                           taxa_rank = "OTU", 
                           userMeta = F,
                           metaInput = NULL,
                           anno_list = NULL, 
                           anno_names = NULL, 
                           dend_col = T, 
                           dend_row = T, 
                           hp_colr = NULL, 
                           anno_colr = NULL, 
                           rowsplit = NULL, 
                           clusterRowSlice = T, 
                           colsplit=NULL, 
                           clusterColSlice = T, 
                           change_colnames = NULL,
                           ...){
        
        # plot heat map
        require(ComplexHeatmap)
        require(phyloseq)
        require(microbiome)
        require(tidyverse)
        require(RColorBrewer)
        require(viridis)


        
        cat("Preparing OTU matrix\n")
        # process data
        if(class(obj) != "phyloseq"){
                stop("Not a phyloseq object")
        }
        
        if(!taxa_are_rows(obj)){
                warning("Taxa are not rows")
                cat("Making taxa as rows...\n")
                otu_table(obj) = otu_table(as(t(otu_table(obj), "matrix")))
        }
        
        
        phy = obj
        rankNames = rank_names(phy)
        

        # collapse taxa
        if(taxa_rank %in% rankNames){
                        phy = phy %>%
                                tax_glom(taxa_rank, NArm = F)}
                
        
        if(norm){
                phy = microbiome::transform(phy, transform = transform, target = target, shift = shift, scale=scale)
                if(percent){
                        
                        phy = phy %>%
                                transform_sample_counts(function(x)100*x)
                }
                
        }
        
        # turn OTU into matrix
        
        # remove OTU not in the samples
        phy = phy %>%
                filter_taxa(function(x)sum(x)>0, prune = T)
        

        otu_mtx = as(otu_table(phy), "matrix")
        otuNames = taxa_names(phy)
        
        
        # apply function to each feature so we can order them
        otuFun = apply(otu_mtx, 1, fun)
        # select show_n
        uniqueTaxa_n = length(unique(tax_table(phy)[, taxa_rank]))
        if(show_n > uniqueTaxa_n){
                stop("show_n larger than the number of the taxa, pick a smaller show_n")
        }
        
        keepOTU = otuNames[order(otuFun, decreasing = T)[1:show_n]]
        
        
        taxTable = tax_table(phy)[keepOTU,]
        otu_mtx_2 = otu_mtx[keepOTU,]
        
        
        
        # transform for heatmap plotting if neccessary
        if(!is.null(hp_transform)){
                if(class(hp_transform) != "function"){
                        stop("hp_transform must be a function")
                }
                
                otu_mtx_2 = hp_transform(otu_mtx_2)
        }
        
        if(userMeta){
                if(is.null(metaInput)){
                        stop("Meta input is null. Please provide metadata")
                }
                metaTable = metaInput
        } else {
                metaTable = data.frame(sample_data(phy))
        }
        

        if(taxa_rank %in% rankNames){
                rownames(otu_mtx_2) = reassignNARank(taxaTable = taxTable, taxLevel = taxa_rank, assignMethod)
        }
        
        if(!is.null(change_colnames)){
                if(change_colnames == "seq"){
                        cat("Column names changed to sequences S#\n")
                        newcol_labels = paste0("S", seq(ncol(otu_mtx_2)))
                        colnames(otu_mtx_2) = newcol_labels
                        rownames(metaTable) = newcol_labels
                }
        }

        
        
        # making heat map
        
        cat("Making heatmap\n")
        # pickCol = rev(brewer.pal(11,"RdBu"))
        # colPalette = colorRampPalette(pickCol)
        

        
        # split col or row
        if(!is.null(colsplit)){
                splitCol = metaTable[[colsplit]]
        } else {
                splitCol = NULL
        }
        
        if(!is.null(rowsplit)){
                splitRow = metaTable[[rowsplit]]
        } else  {
                splitRow = NULL
        }
        
        col_ha = NULL
        row_ha = NULL
        if(!is.null(anno_list)){
                
                if(class(anno_list) != "list"){
                        stop("anno_list must be a list")
                }
                
                if(! toupper( names(anno_list) ) %in% toupper ( c("column", "row") )){
                        stop("anno_list format issue, not found 'column' or 'row' in names")
                }
                
                
                
                ColAnn = anno_list[["column"]]
                RowAnn = anno_list[["row"]]
                meta_cols = colnames(metaTable)
                
                
                
                TopColor = NULL
                SideColor = NULL
                userAnnColor = F
                
                if(!is.null(anno_colr)){
                        # check if anno_colr is a list
                        if(class(anno_colr) != "list"){
                                stop("anno_color must be a list")
                        }
                        
                        if(!names(anno_colr) %in% c("topCol", "sidCol")){
                                stop("anno_colr should contain names topCol or sidCol")
                        }
                        
                        TopColor = anno_colr$topCol
                        SideColor = anno_colr$sidCol
                        
                        userAnnColor = T
                        
                }
                
                
                # column
                if(!is.null(ColAnn)){
                        
                        
                        # change annotation color if user not specified
                        if(!userAnnColor){
                                TopColor = generateCols(ColAnn, metaTable)
                        } 
                                
                        
                        if(sum(ColAnn %in% meta_cols) != length(ColAnn)){
                                stop("Annotation group not found in the metadata")
                        }
                        
                        df_col = metaTable[, ColAnn]
                        
                        dots = lapply(ColAnn, as.symbol)
                        
                        dfcol_order = df_col %>%
                                tibble::rownames_to_column() %>%
                                group_by(.dots = dots) %>%
                                arrange(.by_group = T)
                                

                        if( (!is.null(anno_names)) && (names(anno_names) %in% c("cols", "rows"))){
                                if(!is.null(anno_names$cols)){
                                       colnames(df_col) = anno_names$cols
                                       names(TopColor) = anno_names$cols
                                }
                        }
                        

                        
                        col_ha = HeatmapAnnotation(df = df_col, which = "col", col = TopColor)
                        
                        
                } 
                
                # rows
                if(!is.null(RowAnn)){
                        
                        # change annotation color if user not specified
                        if(!userAnnColor){
                                SideColor = generateCols(RowAnn, metaTable)
                        } 
                        
                        if(sum(RowAnn %in% meta_cols != length(RowAnn))){
                                stop("Annotation group not found in the metadata")
                        }
                        
                        df_row = metaTable[, RowAnn]
                        

                        dfrow_order = df_row %>%
                                tibble::rownames_to_column() %>%
                                group_by(.dots = dots) %>%
                                arrange(.by_group = T)

                        
                        
                        
                        
                        if( (!is.null(anno_names)) && (names(anno_names) %in% c("cols", "rows"))){
                                if(!is.null(anno_names$rows)){
                                        colnames(df_row) = anno_names$rows
                                        names(SideColor) = anno_names$rows
                                }
                        }
                        
                        


                        
                        row_ha = HeatmapAnnotation(df = df_col, which="row", col = SideColor)
                        
                }
                
                

        }
        
        

        
        # color
        if(is.null(hp_colr)){
                col_fun = circlize::colorRamp2(seq(min(otu_mtx_2), max(otu_mtx_2),length=show_n), viridis::plasma(show_n, alpha = hp_colr_alpha))
        } else {
                rampFun = colorRampPalette(hp_colr, alpha = hp_colr_alpha)
                col_fun = circlize::colorRamp2(seq(min(otu_mtx_2), max(otu_mtx_2),length=show_n), rampFun(show_n))
        }

        # dendrogram clustering
        if(!is.null(use_hclust) && !is.null(hclust_pos)){
                whole_dist = vegan::vegdist(t(otu_mtx), method=use_hclust)
                sample.clus = as.dendrogram(hclust(whole_dist, method = hclust_method))
                otu_dist = vegan::vegdist(otu_mtx_2, method = use_hclust)
                species.clus = as.dendrogram(hclust(otu_dist, method = hclust_method))
                
                if(! toupper(hclust_pos) %in% toupper(c("Rows", "row", "Cols", "col", "both"))){
                        stop("hclust_pos must be in rows, cols, or both")
                }
                
                if(toupper(hclust_pos) %in% toupper(c("Rows", "row"))){
                        cat("hclust rows...\n")
                        
                        if(!dend_row){
                                
                                
                                cat("Setting dend_row to FALSE overide the hclust\n")
                                
                                
                                hmap = Heatmap(otu_mtx_2, 
                                               top_annotation = col_ha, 
                                               right_annotation = row_ha, 
                                               cluster_columns = dend_col, 
                                               cluster_rows = dend_row,
                                               row_order = dfrow_order$rowname,
                                               col = col_fun, 
                                               column_split = splitCol,
                                               row_split = splitRow,
                                               cluster_column_slices = clusterColSlice, 
                                               cluster_row_slices = clusterRowSlice,
                                               ...
                                )
                        } else {
                                hmap = Heatmap(otu_mtx_2, 
                                               top_annotation = col_ha, 
                                               right_annotation = row_ha, 
                                               cluster_columns = dend_col, 
                                               cluster_rows = species.clus, 
                                               col = col_fun, 
                                               column_split = splitCol,
                                               row_split = splitRow,
                                               cluster_column_slices = clusterColSlice, 
                                               cluster_row_slices = clusterRowSlice,
                                               ...
                                )    
                        }

                } else if(toupper(hclust_pos) %in% toupper(c("Cols", "col"))){
                        cat("hclust cols...\n")
                        
                        if(!dend_col){
                                cat("Setting dend_col to FALSE overide the hclust\n")
                                
                                hmap = Heatmap(otu_mtx_2, 
                                               top_annotation = col_ha, 
                                               right_annotation = row_ha, 
                                               cluster_columns = dend_col,
                                               column_order = dfcol_order$rowname,
                                               cluster_rows = dend_row, 
                                               col = col_fun, 
                                               column_split = splitCol,
                                               row_split = splitRow,
                                               cluster_column_slices = clusterColSlice, 
                                               cluster_row_slices = clusterRowSlice,
                                               ...
                                )
                        } else {
                        
                                hmap = Heatmap(otu_mtx_2, 
                                               top_annotation = col_ha, 
                                               right_annotation = row_ha, 
                                               cluster_columns = sample.clus, 
                                               cluster_rows = dend_row, 
                                               col = col_fun, 
                                               column_split = splitCol,
                                               row_split = splitRow,
                                               cluster_column_slices = clusterColSlice, 
                                               cluster_row_slices = clusterRowSlice,
                                               ...
                                )
                        }
                } else if (toupper(use_hclust) == toupper("both")){
                        
                        cat("hclust both row and col\n")
                        if((!dend_col) & (!dend_row)){
                                cat("Setting dend_col and dend_row to FALSE overide the hclust\n")
                                hmap = Heatmap(otu_mtx_2, 
                                               top_annotation = col_ha, 
                                               right_annotation = row_ha, 
                                               cluster_columns = dend_col,
                                               column_order = dfcol_order$rowname,
                                               row_order = dfrow_order$rowname,
                                               cluster_rows = dend_row, 
                                               col = col_fun, 
                                               column_split = splitCol,
                                               row_split = splitRow,
                                               cluster_column_slices = clusterColSlice, 
                                               cluster_row_slices = clusterRowSlice,
                                               ...  
                                )                              
                        } else {
                                hmap = Heatmap(otu_mtx_2, 
                                               top_annotation = col_ha, 
                                               right_annotation = row_ha, 
                                               cluster_columns = sample.clus, 
                                               cluster_rows = species.clus, 
                                               col = col_fun, 
                                               column_split = splitCol,
                                               row_split = splitRow,
                                               cluster_column_slices = clusterColSlice, 
                                               cluster_row_slices = clusterRowSlice,
                                               ...  
                                )
                        }

                }
                
        } else {
                
                cat("Default complexheatmpa cluster methods\n")
                hmap = Heatmap(otu_mtx_2, 
                               top_annotation = col_ha, 
                               right_annotation = row_ha, 
                               cluster_columns = dend_col, 
                               cluster_rows = dend_row, 
                               col = col_fun, 
                               column_split = splitCol,
                               row_split = splitRow,
                               cluster_column_slices = clusterColSlice, 
                               cluster_row_slices = clusterRowSlice,
                               ...
                )
        }
        
        
        


        cat("Plotting completed.\n")

        return(hmap)
}
