#' @description plot venn diagram for the core OTU from phyloseq objects
#' @param phylo_in Input phyloseq object
#' @param level1_grp Grouping factor 1 default null
#' @param level2_grp Grouping factor 2 default null
#' @param normalise if normalize the count data to relative abundance
#' @param abund_cutoff cut off value for filtering relative abundance in proportion, default 0.0001
#' @param prev_cutoff define core ASVs in how many samples
#' @param category_names define category labels for each circle of the venn diagram
#' @param label_alpha set 0 to rid of background color for each label
#' @param combinedPlots if combine all the venn diagram for each level1 group together
#' @return a list of ASVs in specified groups for plot venn diagram, and ggplot objects
#' @author Xiaoping Li


plotCoreVenn = function(phylo_in, level1_grp = NULL, level2_grp=NULL, normalise=F, abund_cutoff =0, prev_cuoff = 50/100, category_names = NULL, plotCaption = T, label_alpha=0, combinedPlots = T, common.legend = NULL, legend.pos = "right", nRow=NULL, align="hv", nCol=NULL, ...){
        
        require("ggVennDiagram")
        require("ggpubr")
        require("gridExtra")
        
        # plot core OTU shared by groups
        
        # phylo_in: input phyloseq object
        # level1_grp: first factor you want to group
        # level2-grp: second factor you want to group
        # normalise: if turn count into relative abundance within each sample
        # abund_cutoff: remove OTU in each sample less than this proportion
        # prev_cutoff: remove the OTUs with presence across all samples less than this threshold
        
        
        
        # grouping
        if(class(phylo_in) != "phyloseq"){
                cat("Not phyloseq object")
                return()
        }
        
        if(normalise){
                phylo = phylo_in %>%
                        transform_sample_counts(function(x)x/sum(x)) %>%
                        filter_taxa(function(x) mean(x) > abund_cutoff, prune = T) # of all seqs
        } else {
                phylo = phylo_in %>%
                        filter_taxa(function(x) mean(x) > abund_cutoff, prune = T)
        }
        
        meta = as(sample_data(phylo), "data.frame")
        
        

        level1_list = list()
        
        # groupby subset
        
        if(!is.null(level1_grp)){
                
                
                uniq_grp = as.character(unique(meta[, level1_grp]))
                                
                # loop through each facet
                for(i in 1:length(uniq_grp)){
                        
                        phylo_copy = phylo
                        grp = uniq_grp[i]
                        newdf = meta[which(meta[, level1_grp]==grp),]
                        sample_data(phylo_copy) = sample_data(newdf)
                        level1_list[[grp]] = list() # store phyloseq
                        level1_list[[grp]][[1]] = phylo_copy
                        level1_list[[grp]][[2]] = newdf
                        
                }
                
                
                
        } else {
                
                level1_list[[1]] = list()
                level1_list[[1]][[1]] = phylo
                level1_list[[1]][[2]] = meta
                
        }
        
        
        # check second level of grouping, subset each within each facet
        level2_list = list()
        if(!is.null(level2_grp)){
                for(i in 1:length(level1_list)){
                        # each phyloseq subset

                        # get howmany groups on 2nd level
                        
                        
                        facet_phylo = level1_list[[i]][[1]]
                        facet_meta = level1_list[[i]][[2]]
                        facet_name = names(level1_list)[i]


                        

                        groups = as.character(unique(facet_meta[, level2_grp]))

                        for(j in 1:length(groups)){

                                
                                facet_phylo_copy = facet_phylo
                                each_grp = groups[j]

                                l2_meta = facet_meta[which(facet_meta[, level2_grp]==each_grp),]
                                
                                sample_data(facet_phylo_copy) = sample_data(l2_meta)
                                
                                level2_list[[facet_name]][[each_grp]] = list()
                                level2_list[[facet_name]][[each_grp]] = facet_phylo_copy

                        }


                }
        } else {
                level2_list = level1_list
        }
        
        
        # prevalence and abundance
        
        
        core_list = list()
        for(i in 1:length(level2_list)){
                
                l1_name = names(level2_list[i])
                
                core_list[[l1_name]] = list()
                
                for(k in 1:length(level2_list[i][[1]])){
                        
                        each_phylo = level2_list[i][[1]][[k]]
                        
                        l2_name = names(level2_list[i][[1]])[k]
                        
                        # core microbiome defined as present in about 100% of samples
                        ps_core = filter_taxa(each_phylo, function(x) {sum(x>0) >= (prev_cuoff * length(x))}, prune = T)
                        
                        OTUs = taxa_names(ps_core)
                        core_list[[l1_name]][[l2_name]] = OTUs
                        
                }
        }
        
        
        # plotting

        plot_list = list()
        name_core = names(core_list)
        for(i in 1:length(core_list)){
                core_n = core_list[[i]]
                name_n = name_core[i]
                
                if(!is.null(category_names)){
                                venn_core = ggVennDiagram::ggVennDiagram(core_n, category.names = category_names, label_alpha=label_alpha)
                } else {
                        
                                venn_core = ggVennDiagram::ggVennDiagram(core_n, label_alpha=label_alpha)
                                
                }
                        
                
                
                if(plotCaption){
                        vennout = venn_core + labs(caption=name_n) + theme( plot.caption = element_text(hjust = 0.5, ...))
                } else {
                        vennout = venn_core
                }
                
                
                plot_list[[name_n]] = vennout
        }
        
        
        
        if(combinedPlots){
                print(ggarrange(plotlist = plot_list, nrow = nRow, ncol = nCol, align = align, legend = legend.pos, common.legend = common.legend))
        } 
        
        return(list(coreASVs=core_list, vennplots=plot_list))
        
}



