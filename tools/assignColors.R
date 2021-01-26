assignColors = function(df, group, colorset = "Set3", n = 11){
        
        # assig a color colum to the orginal
        
        if(!group %in% colnames(df)){
                stop("Group is not in data frame columns")
        }
        
        N_unique = length(unique(df[, group]))
        uniques = unique(df[, group])
        
        require(RColorBrewer)
        
        colPat = colorRampPalette(brewer.pal(n, colorset))
        getCols = colPat(N_unique)
        
        # sp1 - col1
        # sp2 - col2
        
        vcolors = c()
        for(i in 1:nrow(df)){
                x = df[i, group]
                indx = match(x, uniques)
                c = getCols[indx]
                vcolors = c(vcolors, c)
        }
        
        
        df2 = df
        df2$colrs = vcolors
        
        return(df2)
        
}
