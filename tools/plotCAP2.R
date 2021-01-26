
plotsp = function(cap,
                  env,
                  taxcol,
                  labelsp = F, 
                  labeln = 10, 
                  tax_level = "Label", 
                  axes=c("CAP1", "CAP2"), 
                  spcolor = "darkgreen", 
                  spshape = "triangle", 
                  alpha=0.3,
                  pointsize = 4,
                  xadj=0, 
                  yadj=0, 
                  nudgex=0,
                  nudgey=0,
                  sp.biplot = F,
                  scaling = 3,
                  showCentroidGroup = "all"){
        
        
        cap.sp = ggvegan:::fortify.cca(cap, display="sp")
        cap.sa = ggvegan:::fortify.cca(cap, display="sites") %>% bind_cols(env)
        cap.bp = ggvegan:::fortify.cca(cap, display="bp")
        cap.cn = ggvegan:::fortify.cca(cap, display="cn")
        cap.sp = cbind(cap.sp, tax_cols)
        cap.cn$groups =  gsub("[A-Z|0-9].+", "", cap.cn$Label)
        
        
        
        # make arrows
        if(scaling==3){
                mul = ggvegan:::arrowMul(cap.bp[, axes], rbind(cap.sa[, axes], cap.sp[, axes]))
                arrows = cap.bp[,axes]
                arrows = arrows*mul # scaling the arrow
                arrows$Label = cap.bp$Label # add arrow labels
        } else if(scaling==2) {
                mul = ggvegan:::arrowMul(cap.bp[, axes], cap.sp[, axes])
                arrows = cap.bp[,axes]
                arrows = arrows*mul # scaling the arrow
                arrows$Label = cap.bp$Label # add arrow labels
        } else if(scalling==1){
                mul = ggvegan:::arrowMul(cap.bp[, axes], cap.sa[, axes])
                arrows = cap.bp[,axes]
                arrows = arrows*mul # scaling the arrow
                arrows$Label = cap.bp$Label # add arrow labels        
        }
        
        # get summary data to label x and y axes
        summ = summary(cap)
        # each axis explained variation / each axis explained in constrained variation
        xlab.cap = paste0(axes[1], " (", round(summ$concont$importance[1,1]/summ$tot.chi, 4)*100, "% of total/", round(summ$concont$importance[2,1], 4)*100, "% of constrained)")
        ylab.cap = paste0(axes[2], "(", round(summ$concont$importance[1,2]/summ$tot.chi, 4)*100, "% of total/", round(summ$concont$importance[2,2], 4)*100, "% if constrained)")
        
        totalInertia = cap$tot.chi
        cons_inertia = sum(cap$CCA$eig)
        prop_inertia = 100*cons_inertia/totalInertia
        label_prop = paste("Total inertia: ", round(totalInertia,2), "\nConstrained:", round(cons_inertia,2),"\nPercentage:", round(prop_inertia, 2), "%")
        
        
        
        # plot species
        
        g = ggplot()
        
        if(labelsp){
                
                # label most distant species?
                dist_m = as.matrix(cap.sp[, axes])
                dist_sp = sqrt(dist_m[,1]^2 + dist_m[,2]^2)
                far_sp = names(dist_sp)[order(dist_sp, decreasing = T)]
                
                
                g.sp = g +
                        # 0,0 lines
                        geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                        geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                        
                        # species
                        geom_point(data=cap.sp, aes_string(x=axes[1], y=axes[2]), shape = spshape, color = spcolor, size=2, alpha=alpha) +
                        geom_text_repel(data=cap.sp[ rownames(cap.sp) %in% far_sp[1:labeln], ], aes_string(x=axes[1], y=axes[2], label = tax_level), color = spcolor, nudge_x = nudgex, nudge_y = nudgey) +
                        geom_text(aes(x=max(cap.sp[, axes[1]]) -  xadj, y=max(cap.sp[, axes[2]]) - yadj, label=label_prop, hjust = 0)) +
                        coord_fixed() +
                        theme_bw() +
                        labs(color="Species", x = xlab.cap, y= ylab.cap)    
        } else {
                
                g.sp = g +
                        # 0,0 lines
                        geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                        geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                        
                        # species
                        geom_point(data=cap.sp, aes_string(x=axes[1], y=axes[2]), shape = spshape, color = spcolor, size=2, alpha=alpha) +
                        geom_text(aes(x=max(cap.sp[, axes[1]]) -  xadj, y=max(cap.sp[, axes[2]]) - yadj, label=label_prop, hjust = 0)) +
                        coord_fixed() +
                        theme_bw() +
                        labs(color="Species", x = xlab.cap, y = ylab.cap)                 
        }
        
        
        if(sp.biplot){
                
                
                
                
                if(showCentroidGroup == "all"){
                        # add arrows (numerical variables)
                        g.sp = g.sp+
                                geom_segment(data=arrows[!arrows$Label %in% cap.cn$Label, ], aes_string(x="0", y="0", xend=axes[1], yend=axes[2]), arrow = arrow(length = unit(0.01, "npc"))) +
                                geom_text_repel(data=arrows[!arrows$Label %in% cap.cn$Label, ], aes_string(x=axes[1], y=axes[2], label="Label")) +
                                #add centroids (categorical variables)
                                geom_point(data=cap.cn, aes_string(x = axes[1], y=axes[2], color = "Label", shape = "groups"), size=pointsize, alpha = 1) +
                                labs(color = "Centroids", shape = "Categorical variables")  
                } else {
                        
                        if(! tolower(showCentroidGroup) %in% tolower(cap.cn$groups)){
                                stop("Centroid group not found")
                        }
                        
                        if(length(showCentroidGroup) > 1){
                                stop("Can only show 1 centroid, or you can specify show all")
                        }
                        
                        cap.cn.select = cap.cn[ which(tolower(cap.cn$groups) %in% tolower(showCentroidGroup)), ]
                        
                        # add arrows (numerical variables)
                        g.sp = g.sp+
                                geom_segment(data=arrows[!arrows$Label %in% cap.cn$Label, ], aes_string(x="0", y="0", xend=axes[1], yend=axes[2]), arrow = arrow(length = unit(0.01, "npc"))) +
                                geom_text_repel(data=arrows[!arrows$Label %in% cap.cn$Label, ], aes_string(x=axes[1], y=axes[2], label="Label")) +
                                #add centroids (categorical variables)
                                geom_point(data=cap.cn.select, aes_string(x = axes[1], y=axes[2]), shape=15, size=pointsize, alpha = 1)+
                                geom_text_repel(data=cap.cn.select, aes_string(x = axes[1], y=axes[2], label="Label"))
                        
                        
                        
                        }
                }
                
        
        return(g.sp)
}


plotsa = function(cap,
                  env,
                  tax_cols,
                  axes = c("CAP1", "CAP2"), 
                  alpha = 0.5, 
                  pointsize = 3, 
                  samplecolor = "grey45", 
                  ellipse = F, 
                  sample_group=NULL, 
                  colorSamples = F,
                  labelsa = F,
                  labelby = "Label",
                  xadj=0, 
                  yadj=0,
                  sa.biplot = F,
                  scaling = 3,
                  showCentroidGroup = "all",
                  addSp = F,
                  labelsp = F,
                  labeln = 10,
                  spshape = "triangle",
                  spcolor = "darkgreen",
                  tax_level = "Label",
                  nudgex = 0,
                  nudgey = 0,
                  spalpha = 0.3){
        
        # sites are grey
        # categorical variables will be plotted as centroids
        # numerical variables will be plotted as arrows

        
        
        cap.sp = ggvegan:::fortify.cca(cap, display="sp")
        cap.sa = ggvegan:::fortify.cca(cap, display="sites") %>% bind_cols(env)
        cap.bp = ggvegan:::fortify.cca(cap, display="bp")
        cap.cn = ggvegan:::fortify.cca(cap, display="cn")
        cap.sp = cbind(cap.sp, tax_cols)
        cap.cn$groups =  gsub("[A-Z|0-9].+", "", cap.cn$Label)
        
        
        
        # make arrows
        if(scaling==3){
                mul = ggvegan:::arrowMul(cap.bp[, axes], rbind(cap.sa[, axes], cap.sp[, axes]))
                arrows = cap.bp[,axes]
                arrows = arrows*mul # scaling the arrow
                arrows$Label = cap.bp$Label # add arrow labels
        } else if(scaling==2) {
                mul = ggvegan:::arrowMul(cap.bp[, axes], cap.sp[, axes])
                arrows = cap.bp[,axes]
                arrows = arrows*mul # scaling the arrow
                arrows$Label = cap.bp$Label # add arrow labels
        } else if(scalling==1){
                mul = ggvegan:::arrowMul(cap.bp[, axes], cap.sa[, axes])
                arrows = cap.bp[,axes]
                arrows = arrows*mul # scaling the arrow
                arrows$Label = cap.bp$Label # add arrow labels        
        }
        
        # get summary data to label x and y axes
        summ = summary(cap)
        # each axis explained variation / each axis explained in constrained variation
        xlab.cap = paste0(axes[1], " (", round(summ$concont$importance[1,1]/summ$tot.chi, 4)*100, "% of total/", round(summ$concont$importance[2,1], 4)*100, "% of constrained)")
        ylab.cap = paste0(axes[2], " (", round(summ$concont$importance[1,2]/summ$tot.chi, 4)*100, "% of total/", round(summ$concont$importance[2,2], 4)*100, "% if constrained)")
        
        totalInertia = cap$tot.chi
        cons_inertia = sum(cap$CCA$eig)
        prop_inertia = 100*cons_inertia/totalInertia
        label_prop = paste("Total inertia: ", round(totalInertia,2), "\nConstrained:", round(cons_inertia,2),"\nPercentage:", round(prop_inertia, 2), "%")
        
        
        
        g = ggplot()
        
        # grouping or not
        if(colorSamples){
                
                # plot ellipse or not
                if(ellipse){
                        # 0,0 lines
                        g.sa = g + 
                                geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                # samples
                                geom_point(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), size=pointsize, alpha=alpha) +
                                stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "norm")+
                                stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "t", linetype=2) +
                                geom_text(aes(x=max(cap.sa[, axes[1]]) -  xadj, y=max(cap.sa[, axes[2]]) - yadj, label=label_prop, hjust = 0)) +
                                coord_fixed() +
                                theme_bw() +
                                labs(color="Sites", x = xlab.cap, y=ylab.cap) 
                } else {
                        g.sa = g + 
                                geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                # samples
                                geom_point(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), size=pointsize, alpha=alpha) +
                                
                                geom_text(aes(x=max(cap.sa[, axes[1]]) - xadj, y=max(cap.sa[, axes[2]]) - yadj, label=label_prop, hjust = 0)) +
                                coord_fixed() +
                                theme_bw() +
                                labs(color="Sites", x = xlab.cap, y=ylab.cap) 
                }
                
                          
        } else {
                g.sa = g +
                        # 0,0 lines
                        geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                        geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                        # samples
                        geom_point(data=cap.sa, aes_string(x=axes[1], y=axes[2]), color=samplecolor, size=pointsize, alpha=alpha) +
                        geom_text(aes(x=max(cap.sa[, axes[1]]) - xadj, y=max(cap.sa[, axes[2]]) - xadj, label=label_prop, hjust = 0)) +
                        coord_fixed() +
                        theme_bw() +
                        labs(color="Sites", x = xlab.cap, y=ylab.cap) 
        }
        
        
        if(labelsa){
                g.sa = g.sa +
                        geom_text(data = cap.sa, aes_string(x = axes[1], y = axes[2], label = labelby), size = 2)
        }
        
        
        
        if(addSp){
                if(labelsp){
                        
                        # label most distant species?
                        dist_m = as.matrix(cap.sp[, axes])
                        dist_sp = sqrt(dist_m[,1]^2 + dist_m[,2]^2)
                        far_sp = names(dist_sp)[order(dist_sp, decreasing = T)]
                        
                        g.sa = g.sa +
                        # species
                        geom_point(data=cap.sp, aes_string(x=axes[1], y=axes[2]), shape = spshape, color = spcolor, size=2, alpha=spalpha) +
                        geom_text_repel(data=cap.sp[ rownames(cap.sp) %in% far_sp[1:labeln], ], aes_string(x=axes[1], y=axes[2], label = tax_level), color = spcolor, nudge_x = nudgex, nudge_y = nudgey)
                        } else {
                                g.sa = g.sa +
                                # species
                                geom_point(data=cap.sp, aes_string(x=axes[1], y=axes[2]), shape = spshape, color = spcolor, size=2, alpha=0.3)
                        }
                
        }
        
        
        
        if(sa.biplot){
                
                # add arrows (numerical variables)
                if(showCentroidGroup == "all"){
                        # add arrows (numerical variables)
                        g.sa = g.sa+
                                geom_segment(data=arrows[!arrows$Label %in% cap.cn$Label, ], aes_string(x="0", y="0", xend=axes[1], yend=axes[2]), arrow = arrow(length = unit(0.01, "npc"))) +
                                geom_text_repel(data=arrows[!arrows$Label %in% cap.cn$Label, ], aes_string(x=axes[1], y=axes[2], label="Label")) +
                                #add centroids (categorical variables)
                                geom_point(data=cap.cn, aes_string(x = axes[1], y=axes[2], color = "Label", shape = "groups"), size=pointsize, alpha = 1) +
                                labs(color = "Centroids", shape = "Categorical variables")  
                } else {
                        
                        if(! tolower(showCentroidGroup) %in% tolower(cap.cn$groups)){
                                stop("Centroid group not found")
                        }
                        
                        if(length(showCentroidGroup) > 1){
                                stop("Can only show 1 centroid, or you can specify show all")
                        }
                        
                        cap.cn.select = cap.cn[ which(tolower(cap.cn$groups) %in% tolower(showCentroidGroup)), ]
                        
                        # add arrows (numerical variables)
                        g.sa = g.sa+
                                geom_segment(data=arrows[!arrows$Label %in% cap.cn$Label, ], aes_string(x="0", y="0", xend=axes[1], yend=axes[2]), arrow = arrow(length = unit(0.01, "npc"))) +
                                geom_text_repel(data=arrows[!arrows$Label %in% cap.cn$Label, ], aes_string(x=axes[1], y=axes[2], label="Label")) +
                                #add centroids (categorical variables)
                                geom_point(data=cap.cn.select, aes_string(x = axes[1], y=axes[2]), shape=15, size=pointsize, alpha = 1) +
                                geom_text_repel(data=cap.cn.select, aes_string(x = axes[1], y=axes[2], label="Label"))
                        
                        
                        
                }
        }
        
        
        return(g.sa)
}



plotCAP = function(cap, 
                   env,
                   tax_cols,
                   display = 3, 
                   splitDisplay = F, 
                   commonLeg = T, 
                   leg = "bottom",                   
                   sa.axes = c("CAP1", "CAP2"), 
                   sa.alpha = 0.5, 
                   sa.pointsize = 3, 
                   sa.samplecolor = "grey45", 
                   sa.ellipse = F, 
                   sa.sample_group=NULL, 
                   sa.colorSamples = F,
                   sa.labelsa = F,
                   sa.labelby = "Label",
                   sa.xadj=0, 
                   sa.yadj=0,
                   sa.biplot = F,
                   sa.scaling = 3,
                   sa.showCentroidGroup = "all",
                   sa.addSp = F,
                   sa.labelsp = F,
                   sa.labeln = 10,
                   sa.spshape = "triangle",
                   sa.spcolor = "darkgreen", 
                   sa.tax_level = "Label",
                   sa.nudgex = 0,
                   sa.nudgey = 0,
                   sa.spalpha = 0.3, ...){
        
        require("ggrepel")
        require("ggpubr")
        # this function is used to plot capscale results
        # return a ggplot object
        
        # cap: capscale output
        # phylo: phyloseq object
        # sample_group: group to color for samples
        # tax_level: tax_level to label
        
        
        # display
        # 1: sites
        # 2: species
        # 3: biplot
        
        # praparing data

        if(display == 1){
                # sample only
                g.out = plotsa(cap, 
                               env, 
                               tax_cols,
                               axes = sa.axes, 
                               alpha = sa.alpha, 
                               pointsize = sa.pointsize, 
                               samplecolor = sa.samplecolor, 
                               ellipse = sa.ellipse, 
                               sample_group = sa.sample_group, 
                               colorSamples = sa.colorSamples,
                               labelsa = sa.labelsa,
                               labelby = sa.labelby,
                               xadj = sa.xadj, 
                               yadj = sa.yadj,
                               sa.biplot = sa.biplot,
                               scaling = sa.scaling,
                               showCentroidGroup = sa.showCentroidGroup,
                               addSp = F,
                               sa.labelsp,
                               sa.labeln,
                               sa.spshape,
                               sa.spcolor,
                               sa.tax_level,
                               sa.nudgex,
                               sa.nudgey,
                               sa.spalpha)
        } else if(display == 2){
                # species only
                g.out = plotsp(cap, 
                               env, 
                               tax_cols, 
                               ...)
                
        } else if(display == 3){
                if(splitDisplay){
                        g.sa = plotsa(cap, 
                                      env,
                                      tax_cols,
                                      axes = sa.axes, 
                                      alpha = sa.alpha, 
                                      pointsize = sa.pointsize, 
                                      samplecolor = sa.samplecolor, 
                                      ellipse = sa.ellipse, 
                                      sample_group = sa.sample_group, 
                                      colorSamples = sa.colorSamples,
                                      labelsa = sa.labelsa,
                                      labelby = sa.labelby,
                                      xadj = sa.xadj, 
                                      yadj = sa.yadj,
                                      sa.biplot = sa.biplot,
                                      scaling = sa.scaling,
                                      showCentroidGroup = sa.showCentroidGroup,
                                      addSp = F,
                                      sa.labelsp,
                                      sa.labeln,
                                      sa.spshape,
                                      sa.spcolor,
                                      sa.tax_level,
                                      sa.nudgex,
                                      sa.nudgey,
                                      sa.spalpha)
                        
                        g.sp = plotsp(cap, env, tax_cols, ...)
                        g.out = ggarrange(g.sa, g.sp, align = "hv", common.legend = commonLeg, legend = leg, labels = c("A", "B"))
                } else {
                        g.out = plotsa( cap, 
                                        env,
                                        tax_cols,
                                        axes = sa.axes, 
                                        alpha = sa.alpha, 
                                        pointsize = sa.pointsize, 
                                        samplecolor = sa.samplecolor, 
                                        ellipse = sa.ellipse, 
                                        sample_group = sa.sample_group, 
                                        colorSamples = sa.colorSamples,
                                        labelsa = sa.labelsa,
                                        labelby = sa.labelby,
                                        xadj = sa.xadj, 
                                        yadj = sa.yadj,
                                        sa.biplot = sa.biplot,
                                        scaling = sa.scaling,
                                        showCentroidGroup = sa.showCentroidGroup,
                                        addSp = T,
                                        sa.labelsp,
                                        sa.labeln,
                                        sa.spshape,
                                        sa.spcolor,
                                        sa.tax_level ,
                                        sa.nudgex,
                                        sa.nudgey,
                                        sa.spalpha)
                }
        }
        
        
        return(g.out)
}