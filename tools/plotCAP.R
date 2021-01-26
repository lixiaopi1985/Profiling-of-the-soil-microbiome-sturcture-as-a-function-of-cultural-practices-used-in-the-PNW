





plotCAP = function(cap, phylo, sample_group, tax_level = "Phylum", display = 3, split=F, legendPosition = "bottom", alignplot = "hv",  axes=c("CAP1", "CAP2"), scaling=3, labelsp=T, alpha=0.3, shape=NULL, nudgex = 0, nudgey = 0){
        
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
        env = as(sample_data(phylo), "data.frame")
        tax_cols = as.data.frame(tax_table(phylo))
        tax_cols$ASV = rownames(tax_cols)
        
        
        cap.sp = ggvegan:::fortify.cca(cap, display="sp")
        cap.sa = ggvegan:::fortify.cca(cap, display="sites") %>% bind_cols(env)
        cap.bp = ggvegan:::fortify.cca(cap, display="bp")
        cap.cn = ggvegan:::fortify.cca(cap, display="cn")
        
        cap.sp = cbind(cap.sp, tax_cols)
        
        
        
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
        xlab.cap = paste0("CAP1 (", round(summ$concont$importance[1,1]/summ$tot.chi, 4)*100, "% of total/", round(summ$concont$importance[2,1], 4)*100, "% of constrained)")
        ylab.cap = paste0("CAP2 (", round(summ$concont$importance[1,2]/summ$tot.chi, 4)*100, "% of total/", round(summ$concont$importance[2,2], 4)*100, "% if constrained)")
        
        totalInertia = cap$tot.chi
        cons_inertia = sum(cap$CCA$eig)
        prop_inertia = 100*totalInertia/cons_inertia
        label_prop = paste("Total inertia: ", round(totalInertia,2), "\nConstrained explained:", round(prop_inertia, 2), "%")
        

        
        
        g = ggplot()
        
        if(display==3 | display=="biplot"){
                
                if(isTRUE(split)){
                        if(labelsp){
                                g.sp = g +
                                        # 0,0 lines
                                        geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                        geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                        
                                        # species
                                        geom_point(data=cap.sp, aes_string(x=axes[1], y=axes[2]), shape="triangle", color="darkgreen", size=2, alpha=alpha) +
                                        geom_text_repel(data=cap.sp[ rownames(cap.sp) %in% far_sp[1:10], ], aes_string(x=axes[1], y=axes[2], label = tax_level), color = "darkgreen", nudge_x = nudgex, nudge_y = nudgey) +
                                        # biplots
                                        geom_segment(data=arrows, aes_string(x="0", y="0", xend=axes[1], yend=axes[2]), arrow = arrow(length = unit(0.01, "npc"))) +
                                        geom_text_repel(data=arrows, aes_string(x=axes[1], y=axes[2], label="Label")) +
                                        geom_label(aes(x=1, y=1, label = label_prop)) +
                                        coord_fixed() +
                                        theme_bw() +
                                        labs(color="Sites", x = xlab.cap, y=ylab.cap)
                                
                                g.sa = g +
                                        # 0,0 lines
                                        geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                        geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                        # samples
                                        geom_point(data = cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), size=2, alpha=alpha) +
                                        stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "norm")+
                                        stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "t", linetype=2) +
                                        # biplots
                                        geom_segment(data=arrows, aes_string(x="0", y="0", xend=axes[1], yend=axes[2]), arrow = arrow(length = unit(0.01, "npc"))) +
                                        geom_text_repel(data=arrows, aes_string(x=axes[1], y=axes[2], label="Label")) +
                                        geom_label(aes(x=1, y=1, label = label_prop)) +
                                        coord_fixed() +
                                        theme_bw() +
                                        labs(color="Sites", x = xlab.cap, y=ylab.cap)
                                
                                g.out = ggarrange(g.sa, g.sp, labels = c("A", "B"), common.legend = T, legend = legendPosition, align = alignplot)
                                
                        } else {
                                g.sp = g +
                                        # 0,0 lines
                                        geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                        geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                        geom_point(data=cap.sp, aes_string(x=axes[1], y=axes[2]), shape="triangle", color="darkgreen", size=2, alpha=alpha) +
                                        geom_label(aes(x=1, y=1, label = label_prop)) +
                                        coord_fixed() +
                                        theme_bw() +
                                        labs(color="Sites", x = xlab.cap, y=ylab.cap)
                                
                                
                                g.sa = g +
                                        # 0,0 lines
                                        geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                        geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                        # samples
                                        geom_point(data = cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), size=2, alpha=alpha) +
                                        stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "norm")+
                                        stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "t", linetype=2)+
                                        
                                        # biplots
                                        geom_segment(data=arrows, aes_string(x="0", y="0", xend=axes[1], yend=axes[2]), arrow = arrow(length = unit(0.01, "npc"))) +
                                        geom_text_repel(data=arrows, aes_string(x=axes[1], y=axes[2], label="Label")) +
                                        geom_label(aes(x=1, y=1, label = label_prop)) +
                                        coord_fixed() +
                                        theme_bw() +
                                        labs(color="Sites", x = xlab.cap, y=ylab.cap)
                                
                                g.out = ggarrange(g.sa, g.sp, labels=c("A", "B"), common.legend = T, legend = legendPosition, align = alignplot)
                                
                        }
                        
                } else {
                
                if(labelsp){
                        g.out = g +
                                # 0,0 lines
                                geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                # samples
                                geom_point(data = cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), size=2, alpha=alpha) +
                                stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "norm")+
                                stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "t", linetype=2)+
                                # species
                                geom_point(data=cap.sp, aes_string(x=axes[1], y=axes[2]), shape="triangle", color="darkgreen", size=2, alpha=alpha) +
                                geom_text_repel(data=cap.sp[ rownames(cap.sp) %in% far_sp[1:10], ], aes_string(x=axes[1], y=axes[2], label = tax_level), color = "darkgreen", nudge_x = nudgex, nudge_y = nudgey) +
                                # biplots
                                geom_segment(data=arrows, aes_string(x="0", y="0", xend=axes[1], yend=axes[2]), arrow = arrow(length = unit(0.01, "npc"))) +
                                geom_text_repel(data=arrows, aes_string(x=axes[1], y=axes[2], label="Label")) +
                                geom_label(aes(x=1, y=1, label = label_prop)) +
                                coord_fixed() +
                                theme_bw() +
                                labs(color="Sites", x = xlab.cap, y=ylab.cap)
                        
                } else {
                        g.out = g +
                                # 0,0 lines
                                geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                # samples
                                geom_point(data = cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), size=2, alpha=alpha) +
                                stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "norm")+
                                stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "t", linetype=2)+
                                # species
                                geom_point(data=cap.sp, aes_string(x=axes[1], y=axes[2]), shape="triangle", color="darkgreen", size=2, alpha=alpha) +
                               
                                # biplots
                                geom_segment(data=arrows, aes_string(x="0", y="0", xend=axes[1], yend=axes[2]), arrow = arrow(length = unit(0.01, "npc"))) +
                                geom_text_repel(data=arrows, aes_string(x=axes[1], y=axes[2], label="Label")) +
                                geom_label(aes(x=1, y=1, label = label_prop)) +
                                coord_fixed() +
                                theme_bw() +
                                labs(color="Sites", x = xlab.cap, y=ylab.cap)
                        
                        
                }
        }
                
                
                
        } else if(display==2|display=="species"){
                
                if(labelsp){
                        g.out = g +
                                # 0,0 lines
                                geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                
                                # species
                                geom_point(data=cap.sp, aes_string(x=axes[1], y=axes[2], color=tax_level), size=2, alpha=alpha) +
                                geom_text_repel(data=cap.sp[ rownames(cap.sp) %in% far_sp[1:10], ], aes_string(x=axes[1], y=axes[2], label = tax_level)) +
                                geom_label(aes(x=1, y=1, label = label_prop)) +
                                coord_fixed() +
                                theme_bw() +
                                labs(color=tax_level, x = xlab.cap, y=ylab.cap)
                } else {
                        g.out = g +
                                # 0,0 lines
                                geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                
                                # species
                                geom_point(data=cap.sp, aes_string(x=axes[1], y=axes[2], color=tax_level), size=2, alpha=alpha) +
                                geom_label(aes(x=1, y=1, label = label_prop)) +
                                coord_fixed() +
                                theme_bw() +
                                labs(color=tax_level, x = xlab.cap, y=ylab.cap)
                                
                }
                
                        
                
        } else if(display==1|display=="sites"|display=="samples"){
                
                if(!is.null(shape)){
                        g.out = g +
                                # 0,0 lines
                                geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                # samples
                                geom_point(data = cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group, shape = shape), size=2, alpha=alpha) +
                                stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "norm")+
                                stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "t", linetype=2) +                                
                                geom_label(aes(x=1, y=1, label = label_prop)) +
                                coord_fixed() +
                                theme_bw() +
                                labs(color="Sites", shape=shape, x = xlab.cap, y=ylab.cap)
                        
                
                } else {
                        g.out = g +
                                # 0,0 lines
                                geom_hline(aes(yintercept=0), color="#7D7D7D", linetype=2) +
                                geom_vline(aes(xintercept=0), color="#7D7D7D", linetype=2) +
                                # samples
                                geom_point(data = cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), size=2, alpha=alpha) +
                                stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "norm")+
                                stat_ellipse(data=cap.sa, aes_string(x=axes[1], y=axes[2], color=sample_group), type = "t", linetype=2)+               
                                geom_label(aes(x=1, y=1, label = label_prop)) +
                                coord_fixed() +
                                theme_bw() +
                                labs(color="Sites", x = xlab.cap, y=ylab.cap)
                        
                        
                }
                
                
        }
        return(g.out)
       
}
