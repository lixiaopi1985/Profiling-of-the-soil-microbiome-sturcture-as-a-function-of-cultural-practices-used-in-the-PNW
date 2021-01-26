output_phylo = function(phylo, otuoutfmt = "biom", otuout_transpose = F, collapse_tax = F, getAbund=T,  methods = NULL, scale = 100, selectCols = c("OTU", "Abundance"), metaCols = NULL, groups = "OTU", outputdir = NULL, outfilenames = c("otu", "tax", "metasamples"), indexCol = "SampleID", coremember = T, ...){
        require(phyloseq)
        require(biomformat)
        require(microbiome)
        require(tidyverse)
        require(rlang)
        require(pime) 
        
        outfolder = "phylo_out"
        saveOUT = ""
        
        bycol = syms(selectCols)
        bygroup = syms(groups)
        
        if(!is.null(metaCols)){
                metacol = sym(metaCols)
        }
        
        
        
        if(is.null(outputdir)){
                outputdir = getwd()
                saveOUT = normalizePath(file.path(outputdir, outfolder))
                if(!dir.exists(saveOUT)){
                        cat("\nCreating save folder ... \n")
                        dir.create(saveOUT, recursive = T)
                } else {
                        cat("\n", paste(saveOUT, "existed!\n"))
                }
                
                
        } else {
                
                saveOUT = normalizePath(file.path(outputdir, outfolder))
                if(!dir.exists(saveOUT)){
                        cat("Creating save folder ... \n")
                        dir.create(saveOUT, recursive = T)
                } else {
                        cat(paste(saveOUT, "existed!\n"))
                }
                
        }
        
        # output OTU format :
        # 1. biom
        # 2. tsv
        
        # output taxonomy format:
        # collapsed taxonomy lineaage
        # non collapsed taxonomy lineaage
        
        # output metadata format:
        # melted and non melted
        
        # features can be saved in an specified output folder
        
        # to select best prevalence
        
        
        if(otuout_transpose){
                otu_m = as(t(otu_table(phylo)), "matrix")

        } else {
                otu_m = as(otu_table(phylo), "matrix")

        }


        otu_cols = colnames(otu_m)
        otu_df = as.data.frame(otu_m, stringsAsFactors = F)






        tax = as(tax_table(phylo), "matrix")
        tax[is.na(tax)] = "unassigned"
        tax_cols = colnames(tax)
        tax_df = as.data.frame(tax, stringsAsFactors=F)
        tax_df$OTUID = rownames(tax_df)

        # meta
        sample_df = data.frame(sample_data(phylo))

        # otu table
        if(otuoutfmt == "biom"){
                cat("write biom...\n")
                otu_biom = make_biom(data=otu_m)
                outfile = normalizePath(file.path(saveOUT, paste0(outfilenames[1], ".biom")))
                write_biom(otu_biom, outfile)
        } else if (otuoutfmt == "tsv") {
                cat("write otu matrix...\n")
                otu_df$OTUID = rownames(otu_df)
                otu_df = otu_df[c("OTUID", otu_cols)]
                outfile = normalizePath(file.path(saveOUT, paste0(outfilenames[1], ".tsv")))
                write.table(otu_df, file = outfile, sep = "\t", row.names = F, col.names = T, quote = F)
        }


        if(collapse_tax){
                cat("Writing collapsed taxonomy ...\n")
                outfile_tax = normalizePath(file.path(saveOUT, paste0(outfilenames[2], ".tsv")))
                write.table(tax_df, file=outfile_tax, sep="\t", row.names = F, col.names = T, quote = F)
        } else {

                cat("writing taxonomy lineages ... \n")
                tax_df$taxonomy = do.call(paste, c(tax_df[tax_cols], sep=";"))
                tax_df = tax_df[, c("OTUID", "taxonomy"), drop = F]
                outfile_tax = normalizePath(file.path(saveOUT, paste0(outfilenames[2], ".tsv")))
                write.table(tax_df, file=outfile_tax, sep="\t", row.names = F, col.names = T, quote = F)
        }




        cat("writing metadata ... \n")
        outfile_sample = normalizePath(file.path(saveOUT, paste0(outfilenames[3], ".tsv")))
        
        if(!is.null(metaCols)){
                sample_df %>%
                        select(!!!metacol) %>%
                        write.table(file=outfile_sample, sep="\t", row.names = F, col.names = T, quote = F)
        } else {
                write.table(sample_df, file=outfile_sample, sep="\t", row.names = F, col.names = T, quote = F)
        }
        
        

        if(getAbund){

                if(is.null(methods)){

                outfile_sample = normalizePath(file.path(saveOUT, paste0("Abundance_", outfilenames[3], ".tsv")))

                cat("writing Abundance metadata (Count)... \n")
                phylo %>%
                        psmelt() %>%
                        select(!!!bycol) %>%
                        group_by(!!!bygroup) %>%
                        summarise(mean_abd = mean(Abundance), total_counts = sum(Abundance)) %>%
                        write.table(file=outfile_sample, sep="\t", row.names = F, col.names = T, quote = F)

                }

                else {
                        outfile_sample = normalizePath(file.path(saveOUT, paste0("Abundance_", outfilenames[3], ".tsv")))

                        cat("writing Abundance metadata (", methods, "), scale by ", scale, " ... \n")
                        phylo %>%
                                microbiome::transform(transform = methods, scale = scale) %>%
                                psmelt() %>%
                                select(!!!bycol) %>%
                                group_by(!!!bygroup) %>%
                                summarise(mean_abd = mean(Abundance)) %>%
                                write.table(file=outfile_sample, sep="\t", row.names = F, col.names = T, quote = F)

                }
        }
        
        if(coremember){
                

                core_M = microbiome::core_members(phylo, ...)
                
                phylo_M = prune_taxa(taxa_names(phylo) %in% core_M, phylo)
                cat("Also calculating core members ...\n")
                # create folder
                dir.create(file.path(saveOUT, "core"), recursive = T)
                
                coreout = normalizePath(file.path(saveOUT, "core"))
                write.csv(data.frame(core_members=core_M), file = file.path(coreout, paste0("core_members_", outfilenames[1], ".csv")), quote = F, row.names = F)
                
                
                if(otuout_transpose){
                        otu_m_core = as(t(otu_table(phylo_M)), "matrix")
                        
                } else {
                        otu_m_core = as(otu_table(phylo_M), "matrix")
                        
                }
                
                
                otu_cols_C = colnames(otu_m_core)
                otu_df_C = as.data.frame(otu_m_core, stringsAsFactors = F)
                
                tax_C = as(tax_table(phylo_M), "matrix")
                tax_C[is.na(tax_C)] = "unassigned"
                tax_cols_C = colnames(tax_C)
                tax_df_C = as.data.frame(tax_C, stringsAsFactors=F)
                tax_df_C$OTUID = rownames(tax_df_C)
                
                # meta
                sample_df_C = data.frame(sample_data(phylo_M))
                
                # otu table
                if(otuoutfmt == "biom"){
                        cat("write biom with CORE members...\n")
                        otu_biom_C = make_biom(data=otu_m_core)
                        outfile_C = normalizePath(file.path(coreout, paste0("CORE_", outfilenames[1], ".biom")))
                        write_biom(otu_biom_C, outfile_C)
                } else if (otuoutfmt == "tsv") {
                        cat("write otu matrix with CORE members...\n")
                        otu_df_C$OTUID = rownames(otu_df_C)
                        otu_df_C = otu_df_C[c("OTUID", otu_cols_C)]
                        outfile_C = normalizePath(file.path(coreout, paste0("CORE_", outfilenames[1], ".tsv")))
                        write.table(otu_df_C, file = outfile_C, sep = "\t", row.names = F, col.names = T, quote = F)
                }
                
                
                if(collapse_tax){
                        cat("Writing collapsed taxonomy with CORE members ...\n")
                        outfile_tax_C = normalizePath(file.path(coreout, paste0("CORE_", outfilenames[2], ".tsv")))
                        write.table(tax_df_C, file=outfile_tax_C, sep="\t", row.names = F, col.names = T, quote = F)
                } else {
                        
                        cat("writing taxonomy lineages with CORE members ... \n")
                        tax_df_C$taxonomy = do.call(paste, c(tax_df_C[tax_cols_C], sep=";"))
                        tax_df_C = tax_df_C[, c("OTUID", "taxonomy"), drop = F]
                        outfile_tax_C = normalizePath(file.path(coreout, paste0("CORE_", outfilenames[2], ".tsv")))
                        write.table(tax_df_C, file=outfile_tax_C, sep="\t", row.names = F, col.names = T, quote = F)
                }
                
                
                
                
                cat("writing metadata with CORE members... \n")
                outfile_sample_C = normalizePath(file.path(coreout, paste0("CORE_", outfilenames[3], ".tsv")))
                if(!is.null(metaCols)){
                        sample_df_C %>%
                                select(!!!metacol) %>%
                                write.table(file=outfile_sample_C, sep="\t", row.names = F, col.names = T, quote = F)
                } else {
                        write.table(sample_df_C, file=outfile_sample_C, sep="\t", row.names = F, col.names = T, quote = F)
                }
                
                if(getAbund){
                        
                        if(is.null(methods)){
                                
                                outfile_sample_C = normalizePath(file.path(coreout, paste0("CORE_abundance_", outfilenames[3], ".tsv")))
                                
                                cat("writing Abundance metadata (Count) with CORE members... \n")
                                phylo_M %>% 
                                        psmelt() %>%
                                        select(!!!bycol) %>%
                                        group_by(!!!bygroup) %>%
                                        summarise(mean_abd = mean(Abundance), total_counts = sum(Abundance)) %>%
                                        write.table(file=outfile_sample_C, sep="\t", row.names = F, col.names = T, quote = F)
                                
                        }
                        
                        else {
                                outfile_sample_C = normalizePath(file.path(coreout, paste0("CORE_abundance_", outfilenames[3], ".tsv")))
                                
                                cat("writing Abundance metadata (", methods, ") with CORE members, scale", scale, " ... \n")
                                phylo %>%
                                        microbiome::transform(transform = methods, scale = scale) %>%
                                        prune_taxa(taxa=taxa_names(phylo) %in% core_M) %>%
                                        psmelt() %>%
                                        select(!!!bycol) %>%
                                        group_by(!!!bygroup) %>%
                                        summarise(mean_abd = mean(Abundance)) %>%
                                        write.table(file=outfile_sample_C, sep="\t", row.names = F, col.names = T, quote = F)
                                
                        }
                }
                
        }
        

        cat("Completed!\n")
        
}


process_bestformat = function(bstphylo, phylo){
        rownames(bstphylo@tax_table@.Data) = rownames(tax_table(phylo))
        rownames(bstphylo@otu_table@.Data) = rownames(otu_table(phylo))
        
        tax_table(bstphylo) = tax_table(bstphylo)[, 1:7]
        
        return(bstphylo)
}


phylo2matrix = function(phylo, transpose = T, norm = T, ...){
        
        if(norm){
                cat("Tranforming otu ... \n")
                phylo_norm = transform(phylo, ...)
        } else {
                
                cat("No tranformation ... \n")
                phylo_norm = phylo
        }
        
        # convert to OTU matrix
        cat("Converting OTU to matrix ... \n")
        
        if(transpose){
                otu_m = as(t(otu_table(phylo_norm)), "matrix")
        } else {
                otu_m = as(otu_table(phylo_norm), "matrix") 
        }
        
        
        
        return(otu_m)
}

# append time to columns
reAssignOTU = function(otupath, metapath, taxapath, rowreplace, suffix, taxcollapsed = F, indexCol = "OTUID"){
        
        # use this function to reconstruct otu table so that OTU label would contain a time point indicator
        require(tidyr)
        otuTB = read.table(otupath, sep="\t", header=T, row.names = indexCol)
        meta = read.table(metapath, sep="\t", header=T)
        
        taxdf = read.table(taxapath, sep="\t", header=T, row.names = indexCol)
        
        if(taxcollapsed){
                taxdf2 = taxdf
        } else {
                taxdf2 = taxdf %>% separate(taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
                
        }
        
        otuTB_merge = merge(t(otuTB), taxdf2, by = "row.names")
        colnames(otuTB_merge)[1] = "OTUID"
        otuTB_merge$OTUID = paste(otuTB_merge$OTUID, suffix, sep="_")
        otuTB_merge$Kingdom = paste(otuTB_merge$Kingdom, suffix, sep="_")
        otuTB_merge$Phylum = paste(otuTB_merge$Phylum, suffix, sep="_")
        otuTB_merge$Class = paste(otuTB_merge$Class, suffix, sep="_")
        otuTB_merge$Order = paste(otuTB_merge$Order, suffix, sep="_")
        otuTB_merge$Family = paste(otuTB_merge$Family, suffix, sep="_")
        otuTB_merge$Genus = paste(otuTB_merge$Genus, suffix, sep="_")
        otuTB_merge$Species = paste(otuTB_merge$Species, suffix, sep="_")
        
        indexDF = otuTB_merge[, c("OTUID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
        
        colnames(otuTB) = paste(colnames(otuTB), suffix, sep="_")
        rownames(otuTB) = meta[[rowreplace]]
        return(list(m = as.matrix(otuTB), indx = indexDF))
        
}

stackOTUMatrix = function(..., tax_level){
        
        # stack up multiple 
        
        inputlist = list(...)
        
        mtrx = lapply(inputlist, function(x)x$m)
        tax_map = lapply(inputlist, function(x)x$indx)
        
        if(length(mtrx) > 1){
                rowN = lapply(mtrx, function(x)rownames(x))
                shared = Reduce(intersect, rowN)
                newMtrix = lapply(mtrx, function(x)x[shared, ])
                newTaxa = lapply(tax_map, function(x)x[,c("OTUID", tax_level)])
                

                stacked_m = do.call(cbind, newMtrix)
                stacked_tax = do.call(rbind, newTaxa)
                
                return(list(m=stacked_m, taxa=stacked_tax))

        } else if(length(mtrx) == 1) {
                return(...)
        }
        
}

