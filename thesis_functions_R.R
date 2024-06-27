celltype_composition <- function(freq, xlabel = "timepoint", condition = "cell_type", legend = NA,
                                 save = NA, cols, order, x.order, condition.order = NA, filename = "", suffix = "",
                                 title = "Cell type composition", width = 10, height = 6){
    
    freq <- set(freq, j = xlabel, value = factor(freq[[xlabel]], levels = x.order))
    if(!is.na(condition.order)){
        freq <- set(freq, j = condition, value = factor(freq[[condition]], levels = condition.order))
    }
    leg_cols <- ifelse(length(unlist(unique(freq[, ..condition]))) < 15, 1, 2)
    
    p <- ggplot(freq, aes_string(x = xlabel, y = "value", fill = condition)) + 
         geom_bar(position = "fill", stat = "identity", show.legend = FALSE) + scale_fill_manual(values = cols) +
         labs(y = "relative frequency", x = "", title = title) + 
         theme_minimal() +
         theme(axis.line = element_line(), axis.text = element_text(size = 16, colour = "black"),
               panel.grid.minor = element_line(colour = "gray95"), axis.text.x = element_text(angle = 45, hjust = 1),
               axis.title = element_text(size = 20), axis.ticks = element_line(),
               plot.title = element_text(hjust = 0.5, size = 20), legend.text = element_text(size = 16),
               legend.title = element_text(size = 16)) + guides(fill = guide_legend(ncol = leg_cols))
        

    if(!is.na(legend)){
        p <- p + geom_bar(position = "fill", stat = "identity", show.legend = TRUE)
        }
    
    if(!is.na(save)){
        print(paste0("Saving to composition_", filename, suffix, ".pdf"))
        ggsave(plot = p, width = width, height = height, filename = paste0("composition_", filename, suffix, ".pdf"))
    }
    return(p)
}

rel_freq_boxplot <- function(freq, ct = "AT2", ct_label = "cell_type", xlabel = "time_point", condition = NA,
                             h = 5, w = 6, black_med = NA, prefix = "",
                             cols, x_order, cond_order = NA, a = 0, lwd = 1.5, save = NA, filename = ""){
    
    freq <- set(freq, j = xlabel, value = factor(freq[[xlabel]], levels = x_order))
    #freq <- set(freq, j = condition, value = factor(freq[[condition]], levels = cond_order))
    
    p <- ggplot(subset(freq, get(ct_label) == ct), aes_string(x = xlabel, y = "value", color = xlabel, fill = xlabel)) + 
        geom_boxplot(width = 0.7, lwd = lwd, alpha = a) + 
        geom_point(pch = 19, position = position_jitter(), size = 2, color = "black") + #position = position_jitterdodge()
        scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
        labs(y = "relative frequency", x = "", title = ct) + 
        theme_minimal() +
                theme(axis.line = element_line(), axis.text = element_text(size = 16), legend.position = "none",
                  panel.grid.minor = element_line(colour = "gray95"), strip.text = element_text(size = 18),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_text(size = 18), axis.ticks = element_line(),
              plot.title = element_text(hjust = 0.5, size = 18))
    
    if(!is.na(black_med)){
        dat <- ggplot_build(p)$data[[1]]
        p <- p + geom_segment(data = dat, aes(x = xmin, xend = xmax, y = middle, yend = middle), 
                              color = "black", inherit.aes = F, size = 2)
        }
    
    if(!is.na(save)){
        filename = gsub("/", "_", gsub(" ", "_", filename))
        print(paste0("Saving to ", prefix, filename))
        ggsave(plot = p, width = w, height = h, filename = paste0(prefix, filename))
    }
    return(p)
}

generate_RecLig_table <- function(markers, logFC_label = "logFC", logfc_thresh = 1, save = NA, filename = ""){
    markers <- markers[which(markers[, logFC_label] > logfc_thresh),]
    
    ## Restrict to receptors and ligands which are also cell type markers
    reclig <- lr_network[union(which(lr_network$from %in% markers$gene), which(lr_network$to %in% markers$gene)), ]
    genes <- unique(c(as.character(lr_network$from), as.character(lr_network$to)))
    ok <- intersect(genes, markers$gene)
    reclig_tmp <- reclig[intersect(which(reclig$from %in% ok), which(reclig$to %in% ok)),]
    
    # Merge both receptor and ligand markers ####
    merged <- do.call(rbind, lapply(1:nrow(reclig_tmp), function(x){
      rec_tmp <- markers[which(markers$gene == reclig_tmp$to[x]),]
      lig_tmp <- markers[which(markers$gene == reclig_tmp$from[x]),]
      #rec_tmp$pair <- paste0(lig_tmp, "_", rec_tmp)
      rec_tmp$pair <- rownames(reclig_tmp)[x]
      lig_tmp$pair <- rownames(reclig_tmp)[x]
      merge(rec_tmp, lig_tmp, by.x = "pair", by.y = "pair")
    }))
    colnames(merged) <- gsub(".x", ".rec", colnames(merged), fixed = T)
    colnames(merged) <- gsub(".y", ".lig", colnames(merged), fixed = T)
    
    ## Only keep these columns
    merged = merged[, c("gene.rec", "logFC.rec", "pval_adj.rec", "cell_type.rec",
                        "gene.lig", "logFC.lig", "pval_adj.lig", "cell_type.lig")]
    
    if(!is.na(save)){
        print(paste("Saving to", filename))
        write.table(merged, file = filename, sep = "\t", row.names = F, quote = F)        
    }
    
    return(merged)
}

## Generate adjacency table of all pairs
generate_adjacency_table <- function(merged){
    tmp <- table(paste(merged$cell_type.lig, merged$cell_type.rec, sep = "|"))
    tmp <- data.frame(do.call(rbind, lapply(names(tmp),function(x) strsplit(x, '|', fixed = T)[[1]])),
                                            as.numeric(tmp), stringsAsFactors = F)
                              
    colnames(tmp) <- c('receptor', 'ligand', 'freq')
    celltypes <- unique(c(tmp[,1], tmp[,2]))
    matr <- matrix(0, length(celltypes), length(celltypes))
    colnames(matr) <- rownames(matr) <- celltypes
    n <- lapply(1:nrow(tmp), function(x){
      try(matr[tmp[x,1], tmp[x,2]] <<- tmp[x, 3]  )
    })
    return(matr)
}
                              
## [18.Juni.21] Generate and save both adjacency matrix and rec-lig pairs with logFCs
## [07.Feb.22] added path as variable
generate_intact_tables <- function(genes, master_tab, reg = "", suffix = "", save = NA, return_intact = T,
                                   path = "Data/interaction_tables/", add_dash = "@"){
    
    ## On all cell type combinations
    all_cts = names(genes)
    intact = data.frame(matrix(0, length(genes), length(genes)), stringsAsFactors = F)
    rownames(intact) = all_cts
    colnames(intact) = all_cts
    pairs = data.frame(sender_ct = character(), receiver_ct = character(), ligand = character(),
                       receptor = character(), pairs = character(), stringsAsFactors = F)
    
    for (sender_ct in all_cts){
        expr_genes_sender = genes[[sender_ct]]
        sender_ct_col = gsub(add_dash, "_", gsub("\\+", "\\\\+", sender_ct))

        for (receiver_ct in all_cts){
            expr_genes_receiver <- genes[[receiver_ct]]
            receiver_ct_col = gsub(add_dash, "_", gsub("\\+", "\\\\+", receiver_ct))
            cur_pairs = data.frame(sender_ct = character(), receiver_ct = character(), ligand = character(),
                                   receptor = character(), pairs = character(), stringsAsFactors = F)
            count = 0

            for(lig in expr_genes_sender){
                cur_intact = unique(lr_network$to[(lr_network$from == lig) & (lr_network$to %in% expr_genes_receiver)])
                count = count + length(cur_intact)

                if(length(cur_intact) > 0){
                    for(rec in cur_intact){
                        cur_pairs <- rbind(cur_pairs, data.frame(sender_ct = sender_ct, receiver_ct = receiver_ct,
                                                                 ligand = lig, receptor = rec, 
                                                                 pairs = paste0(lig, "_", rec), stringsAsFactors = F))
        }}}
            if(nrow(cur_pairs) > 0){
                cur_pairs[, "sender_logFC"] = master_tab[cur_pairs$ligand,
                                                         grep(paste0("log2fc_", sender_ct_col, reg), colnames(master_tab))]
                cur_pairs[, "receiver_logFC"] = master_tab[cur_pairs$receptor, 
                                                           grep(paste0("log2fc_", receiver_ct_col, reg), colnames(master_tab))]
                pairs <- rbind(pairs, cur_pairs)
            }        
            intact[sender_ct, receiver_ct] = count
        }
    }
    if(!is.na(save)){
        filename = paste0(reg, suffix, ".txt")
        write.table(intact, file = paste0(path, "adj_matrix", filename), sep = "\t", quote = F)
        write.table(pairs, file = paste0(path, "reclig_pairs", filename), sep = "\t", quote = F, row.names = F)
        print(paste0("Saving to ", path, "adj_matrix", filename, " and reclig_pairs", filename))
    }
    if(!is.na(return_intact)){
        return(intact)
    }
}
                              
                              
## [03.Juni.21] Replaced by slightly prettier connectome plot version below
plot_interaction_network_obsolete <- function(intact, condition = "Bleo", min_intact = 0, save = NA, filename = "",
                                     node_size = 10, scale_weight = 20, text_size = 1.2, legend = NA){
    
    intact <- data.matrix(intact)
    vertex_col <- condition_colors[condition]
    intact[which(intact < min_intact)] <- 0
    G <- as.undirected(graph.adjacency(intact, weighted = TRUE, mode = "undirected", diag = F))
    E(G)$width <- ma_scale(E(G)$weight) * scale_weight
    E(G)$color <- color.scale(E(G)$weight, extremes = c("lightgray", vertex_col))
    V(G)$label.cex <- text_size

    #node_size <- strength(G)/100
    #order = colnames(intact) for dendrogram ordering of nodes
    coords <- layout_in_circle(G, order = order(colnames(intact))) 
    
    if(!is.na(save)){
        pdf(file = filename, width = 15, height = 10, onefile = FALSE)
        print(paste("Saving to", filename))
    }

    plot(G, edge.width = E(G)$width, layout = coords, vertex.color = condition_nodes_colors[condition],
         vertex.size = node_size, vertex.label.family = "Helvetica", vertex.label.color = "black")
    title(condition, cex.main = 3) # paste(condition, "- log2FC", log), cex.main = 3)
    if(!is.na(legend)){
        color.legend(1.12, -1, 1.2, -0.25, align = "rb", cex = text_size, gradient = "y",
                     legend = c(min(E(G)$weight), round(mean(E(G)$weight), 1), max(E(G)$weight)),
                     colorRampPalette(c("lightgray", vertex_col))(100))
    }
}
   
radian.rescale <- function(x, start = 0, direction = 1) {
      c.rotate <- function(x) (x + start) %% (2 * pi) * direction
      c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
    }

                              
## Network Plotting Function - Edge = Number of Receptor Ligand Interactions across al cell types
ma_scale <- function(x){(x - min(x)) / (max(x) - min(x))}
                              
ma_scale_range <- function(x, max_scale){(x - min(x)) / (max_scale - min(x))}
          
## [06.Dez.21] This function was used in my Thesis, but replaced with version below after exploration in
## /short_scripts/Bleo_Aging/induced_connectome/induced_connectome_meta_V2.ipynb
plot_interaction_network_obsolete <- function(intact, condition = "Bleo", min_intact = 0, save = NA, filename = "", ct_map,
                                     node_size = 10, scale_weight = 20, text_size = 1.75, legend = NA, title = NA,
                                     w = 16, h = 12, main_size = 3, max_scale = NA){
    if(is.na(title)){
        title = condition
    }
    
    intact <- data.matrix(intact)
    dark_col <- condition_colors[condition]
    intact[which(intact < min_intact)] <- 0
    G <- as.undirected(graph.adjacency(intact, weighted = TRUE, mode = "undirected", diag = F))
    
    if(is.na(max_scale)){
        max_scale = max(E(G)$weight)
        plot_legend = c(min(E(G)$weight), round(mean(E(G)$weight), 1), max(E(G)$weight))
    } else{
        plot_legend = c(min(E(G)$weight), max_scale)
    }
    E(G)$width <- ma_scale_range(E(G)$weight, max_scale) * scale_weight
    #print(E(G)$width)
    #print(E(G)$weight)
    #print(str(E(G)$width))
    #print(c(min(E(G)$width), max(E(G)$width)))
    E(G)$color <- color.scale(E(G)$width / scale_weight, xrange = c(0, 1),
                              extremes = c("lightgray", condition_colors[paste0(condition, "_vertex")], dark_col))
    
    node_col <- condition_colors[paste0(condition, "_nodes")]
    coords <- layout_in_circle(G)#, order = order(colnames(intact))) 
    
    ## Shorten some labels
    #V(G)$name[V(G)$name == receiver_ct] <- ""
    V(G)$name <- ifelse(V(G)$name %in% names(ct_map), ct_map[V(G)$name], sub(" ", "\n", V(G)$name))
    V(G)$label.cex <- text_size 
    
    ## Get Node labels outside of Nodes
    node_locs <- radian.rescale(x = 1:ncol(intact), direction = -1, start = 0)
    
    if(!is.na(save)){
        pdf(file = filename, width = w, height = h, onefile = FALSE)
        print(paste("Saving to", filename))
    }

    plot(G, edge.width = E(G)$width, layout = coords,
         vertex.frame.color = node_col, vertex.color = node_col, 
         vertex.label.degree = node_locs, vertex.label.dist = 2.1,
         vertex.size = node_size, vertex.label.family = "Helvetica", vertex.label.color = "black")
    title(title, cex.main = main_size) # paste(condition, "- log2FC", log), cex.main = 3)
    if(!is.na(legend)){
        color.legend(1.1, -1.1, 1.18, -0.6, align = "rb", cex = text_size, gradient = "y",
                     legend = plot_legend,
                     colorRampPalette(c("lightgray", condition_colors[paste0(condition, "_vertex")], dark_col))(100))
    }
}

## [06.Dez.21] change mode to "plus" and change xrange to be set automatically, max edges don’t appear sometimes
plot_interaction_network <- function(intact, condition = "Bleo", min_intact = 0, save = NA, filename = "", ct_map,
                                     node_size = 10, scale_weight = 20, text_size = 1.75, legend = NA, title = NA,
                                     w = 16, h = 12, main_size = 3, max_scale = NA,
                                     xl = 1.1, yb = -1.1, xr = 1.18, yt = -0.6){
    if(is.na(title)){
        title = condition
    }
    
    intact <- data.matrix(intact)
    dark_col <- condition_colors[condition]
    intact[which(intact < min_intact)] <- 0
    G <- as.undirected(graph.adjacency(intact, weighted = TRUE, mode = "plus", diag = F))
    
    if(is.na(max_scale)){
        max_scale = max(E(G)$weight)
        plot_legend = c(min(E(G)$weight), round(mean(E(G)$weight), 1), max(E(G)$weight))
    } else{
        plot_legend = c(min(E(G)$weight), max_scale)
    }
    E(G)$width <- ma_scale_range(E(G)$weight, max_scale) * scale_weight

    max_range = max(E(G)$width / scale_weight) + 0.01 ## New
    E(G)$color <- color.scale(E(G)$width / scale_weight, xrange = c(0, max_range),
                              extremes = c("lightgray", condition_colors[paste0(condition, "_vertex")], dark_col))
    
    node_col <- condition_colors[paste0(condition, "_nodes")]
    coords <- layout_in_circle(G)#, order = order(colnames(intact))) 
    
    ## Shorten some labels
    #V(G)$name[V(G)$name == receiver_ct] <- ""
    V(G)$name <- ifelse(V(G)$name %in% names(ct_map), ct_map[V(G)$name], sub(" ", "\n", V(G)$name))
    V(G)$label.cex <- text_size 
    
    ## Get Node labels outside of Nodes
    node_locs <- radian.rescale(x = 1:ncol(intact), direction = -1, start = 0)
    
    if(!is.na(save)){
        pdf(file = filename, width = w, height = h, onefile = FALSE)
        print(paste("Saving to", filename))
    }

    plot(G, edge.width = E(G)$width, layout = coords,
         vertex.frame.color = node_col, vertex.color = node_col, 
         vertex.label.degree = node_locs, vertex.label.dist = 2.1,
         vertex.size = node_size, vertex.label.family = "Helvetica", vertex.label.color = "black")
    title(title, cex.main = main_size) # paste(condition, "- log2FC", log), cex.main = 3)
    if(!is.na(legend)){
        color.legend(xl, yb, xr, yt, align = "rb", cex = text_size, gradient = "y",
                     legend = plot_legend,
                     colorRampPalette(c("lightgray", condition_colors[paste0(condition, "_vertex")], dark_col))(100))
    }
}


## [08.Mai.21] Calculate and save PCA plots on Synthetic bulks (using voom normalization)
voom_normalize <- function(bulk, cutoff = 1){
    
    ## Create DGEList object and calculate Normalization Factors
    d0 <- DGEList(bulk)
    d0 <- calcNormFactors(d0)

    ## Filter low-expressed genes
    drop <- which(apply(cpm(d0), 1, max) < cutoff)
    d <- d0[-drop,] 
    
    ## Counts are transformed to log2 CPM, where “per million reads” is based on the normalization factors
    bulk <- voom(d)$E
    return(bulk)
}

filter_variable <- function(bulk, n_genes){
    vars <- apply(sqrt(bulk), 1, var)
    ok_genes <- names(vars[order(vars, decreasing = T)])[1:n_genes]
    expr <- as.data.frame(t(bulk[ok_genes,]))
    return(expr)
}

## plot_PCA(bulk.norm, condition = "treatment", order = order, cols = cols, comp = c(1, 2),
##          save = T, filename = "pca_bleo_wholeLung_12.pdf", width = 8, height = 7)
plot_PCA <- function(expr, order, cols, condition = "treatment", comp = c(1, 2), s = 4, 
                     ellipse = F, save = NA, filename = "", width = 8, height = 8){
    expr[, condition] <- factor(treat[, condition], levels = order)
    
    ## Perform PCA
    res.pca <- PCA(expr, graph = FALSE, scale.unit = T, quali.sup = ncol(expr))
    p <- fviz_pca_ind(res.pca, geom = "point", axes = comp, habillage = expr[, condition], mean.point = F,
                      addEllipses = ellipse, palette = cols, pointshape = 19, pointsize = s,
                      title = "Synthetic Bulk PCA")
    p <- p + theme_minimal() +
         theme(axis.line = element_line(), axis.text = element_text(size = 16, colour = "black"),
               panel.grid.minor = element_line(colour = "gray95"), #axis.text.x = element_text(angle = 45, hjust = 1),
               axis.title = element_text(size = 16), axis.ticks = element_line(),
               plot.title = element_text(hjust = 1, size = 16, face = "italic"),
               legend.text = element_text(size = 16), legend.title = element_text(size = 16)) 
    
    if(!is.na(save)){
        print(paste("Saving to", filename))
        ggsave(plot = p, width = width, height = height, filename = filename)
    }
    
    plot(p)
}

## [10.Feb.22] Added this function, to enable directed graphs, in which arrow is coloured by sending ct
plot_interaction_directed <- function(intact, condition = "Bleo", min_intact = 0, save = NA, filename = "", ct_map,
                                     node_size = 10, scale_weight = 20, text_size = 1.75, legend = NA, title = NA,
                                     w = 16, h = 12, main_size = 3, max_scale = NA, ct_cols = NA,
                                     xl = 1.1, yb = -1.1, xr = 1.18, yt = -0.6, arrow_size = 0.6){
    
    if(is.na(title)){title = condition}
    intact <- data.matrix(intact)
    dark_col <- condition_colors[condition]
    
    diag(intact) = 0 ## do not consider autocrine feedback
    intact[which(intact < min_intact)] <- 0
    
    #G <- as.undirected(graph.adjacency(intact, weighted = TRUE, mode = "plus", diag = F))
    G <- graph_from_adjacency_matrix(intact, mode = c("directed"), weighted = TRUE)
    
    if(is.na(max_scale)){
        max_scale = max(E(G)$weight)
        # plot_legend = c(min(E(G)$weight), round(mean(E(G)$weight), 1), max(E(G)$weight))
        plot_legend = c(min_intact, round(mean(E(G)$weight), 1), max(E(G)$weight))
      } else{
        plot_legend = c(min_intact, max_scale)
        }
    
    # Here the edge width are scaled so that all the networks are comparable based on specific scale_weight
    #print(E(G)$weight)
    E(G)$width <- ma_scale_range(E(G)$weight, max_scale) * scale_weight
    #print(E(G)$width)
    
    #max_range = max(E(G)$width / scale_weight) + 0.01 
    #print(max_range)
    
    #E(G)$color <- color.scale(E(G)$width / scale_weight, xrange = c(0, max_range),
    #                          extremes = c("lightgray", condition_colors[paste0(condition, "_vertex")], dark_col))
    
    ## cell type specific colouring necessary
    edges <- sapply(attr(E(G), "vnames"), function(x){strsplit(x, "\\|")[[1]][1]})
    #print(edges)
    E(G)$color <- ct_cols[edges]
    
    #node_col <- condition_colors[paste0(condition, "_nodes")]
    node_col = ct_cols[ct_order]
    coords <- layout_in_circle(G)#, order = order(colnames(intact))) 
    
    ## Shorten some labels
    V(G)$name <- ifelse(V(G)$name %in% names(ct_map), ct_map[V(G)$name], sub(" ", "\n", V(G)$name))
    V(G)$label.cex <- text_size 
    
    ## Get Node labels outside of Nodes
    node_locs <- radian.rescale(x = 1:ncol(intact), direction = -1, start = 0)
    
    if(!is.na(save)){
        pdf(file = filename, width = w, height = h, onefile = FALSE)
        print(paste("Saving to", filename))
    }

   plot(G, edge.width = E(G)$width, layout = coords,
         vertex.frame.color = node_col, vertex.color = node_col, 
         vertex.label.degree = node_locs, vertex.label.dist = 2.1,
         edge.curved = .3, edge.arrow.size = arrow_size,
         vertex.size = node_size, vertex.label.family = "Helvetica", vertex.label.color = "black")
    title(title, cex.main = main_size) # paste(condition, "- log2FC", log), cex.main = 3)
    if(!is.na(legend)){
        color.legend(xl, yb, xr, yt, align = "rb", cex = text_size, gradient = "y",
                     legend = plot_legend, colorRampPalette(c("lightgray", "black"))(100))
    }
    #print(G)
    # Return the adjacency matrix
    #return(as_adjacency_matrix(G, attr = "weight"))  
    #return(plotx)  
}

## [02.März.22] Added plotting code for substractive Plots
plot_difference_directed <- function(intact, old_mat, young_mat, save = NA, filename = "", ct_map, omit_diff = 2.5,
                                     node_size = 10, scale_weight = 20, text_size = 1.75, legend = NA, title = NA,
                                     w = 16, h = 12, main_size = 3, max_scale = NA, ct_cols = NA,
                                     xl = 1.1, yb = -1.1, xr = 1.18, yt = -0.6, arrow_size = 0.6){
    
    if(is.na(title)){title = ""}
    intact <- data.matrix(intact)
    diag(intact) = 0 ## do not consider autocrine feedback
    intact[which((intact > -omit_diff) & (intact < omit_diff))] <- 0
    G <- graph_from_adjacency_matrix(intact, mode = c("directed"), weighted = TRUE)
    
    if(is.na(max_scale)){
        max_scale = max(E(G)$weight)
         plot_legend = c(min(E(G)$weight), 0, max(E(G)$weight))
    } else{
        plot_legend = c(min(E(G)$weight), max_scale)
    }

    #E(G)$width <- ma_scale_range(E(G)$weight, max_scale) * scale_weight
    E(G)$width <- ma_scale_range(abs(E(G)$weight), max_scale) * scale_weight
    
    max_range = max(E(G)$width / scale_weight) + 0.01 
    E(G)$color <- color.scale(E(G)$weight / scale_weight, extremes = c("steelblue3", "gray90", "firebrick"))
    
    node_col = ct_cols[ct_order]
    coords <- layout_in_circle(G)
    
    ## Shorten some labels
    V(G)$name <- ifelse(V(G)$name %in% names(ct_map), ct_map[V(G)$name], sub(" ", "\n", V(G)$name))
    V(G)$label.cex <- text_size 
    
    ## Get Node labels outside of Nodes
    node_locs <- radian.rescale(x = 1:ncol(intact), direction = -1, start = 0)
    
    if(!is.na(save)){
        pdf(file = filename, width = w, height = h, onefile = FALSE)
        print(paste("Saving to", filename))
    }

    plot(G, edge.width = E(G)$width, layout = coords,
         vertex.frame.color = node_col, vertex.color = node_col, 
         vertex.label.degree = node_locs, vertex.label.dist = 2.1,
         edge.curved = .3, edge.arrow.size = arrow_size,
         vertex.size = node_size, vertex.label.family = "Helvetica", vertex.label.color = "black")
    title(title, cex.main = main_size)
    if(!is.na(legend)){
        color.legend(xl, yb, xr, yt, align = "rb", cex = text_size, gradient = "y",
                     legend = plot_legend, colorRampPalette(c("steelblue3", "gray90", "firebrick"))(100))
    }
    
    # Return the modified matrices and the graph object
    # return(list(old_mat = old_mat, young_mat = young_mat, Diff_G = G))
    #return(G)
}

## [02.März.22] Added plotting code for substractive Plots
plot_difference_directed_modified <- function(intact1, intact2, min_intact1, min_intact2, save = NA, filename = "", ct_map, omit_diff = 2.5,
                                     node_size = 10, scale_weight = 20, text_size = 1.75, legend = NA, title = NA,
                                     w = 16, h = 12, main_size = 3, max_scale = NA, ct_cols = NA,
                                     xl = 1.1, yb = -1.1, xr = 1.18, yt = -0.6, arrow_size = 0.6){
  intact1 <- data.matrix(intact1)
  intact2 <- data.matrix(intact2)
  
  diag(intact1) = 0 ## do not consider autocrine feedback
  #intact1[intact1 < min_intact1] <- 0
  intact1[which(intact1 < min_intact1)] <- 0
  
  diag(intact2) = 0 ## do not consider autocrine feedback
  #intact2[intact2 < min_intact2] <- 0
  intact2[which(intact2 < min_intact2)] <- 0
  
  intact=intact1-intact2
  
  #G <- as.undirected(graph.adjacency(intact, weighted = TRUE, mode = "plus", diag = F))
  G <- graph_from_adjacency_matrix(intact, mode = c("directed"), weighted = TRUE)
  
  if(is.na(max_scale)){
    max_scale = max(E(G)$weight)
    plot_legend = c(min(E(G)$weight), 0, max(E(G)$weight))
  } else{
    plot_legend = c(min(E(G)$weight), max_scale)
  }
  
  # Here the edge width are scaled so that all the networks are comparable based on specific scale_weight
  #print(E(G)$weight)
 # E(G)$width <- ma_scale_range(E(G)$weight, max_scale) * scale_weight
  
  E(G)$width <- ma_scale_range(abs(E(G)$weight), max_scale) * scale_weight
  
  #print(E(G)$width)
  ## cell type specific colouring necessary
  edges <- sapply(attr(E(G), "vnames"), function(x){strsplit(x, "\\|")[[1]][1]})
  # print("These are edges")
  # print(edges)
  # print(class(edges))
  # print(ct_cols[edges])
  # 
  #E(G)$color <- ct_cols[edges]
  E(G)$color <- color.scale(E(G)$weight / scale_weight, extremes = c("steelblue3", "gray90", "firebrick"))
  #print(E(G)$color)
  
  #node_col <- condition_colors[paste0(condition, "_nodes")]
  node_col = ct_cols[ct_order]
  coords <- layout_in_circle(G)#, order = order(colnames(intact))) 
  
  ## Shorten some labels
  V(G)$name <- ifelse(V(G)$name %in% names(ct_map), ct_map[V(G)$name], sub(" ", "\n", V(G)$name))
  V(G)$label.cex <- text_size 
  
  ## Get Node labels outside of Nodes
  node_locs <- radian.rescale(x = 1:ncol(intact), direction = -1, start = 0)
  
  if(!is.na(save)){
    pdf(file = filename, width = w, height = h, onefile = FALSE)
    print(paste("Saving to", filename))
  }
  
  plot(G, edge.width = E(G)$width, layout = coords,
       vertex.frame.color = node_col, vertex.color = node_col, 
       vertex.label.degree = node_locs, vertex.label.dist = 2.1,
       edge.curved = .3, edge.arrow.size = arrow_size,
       vertex.size = node_size, vertex.label.family = "Helvetica", vertex.label.color = "black")
  title(title, cex.main = main_size) # paste(condition, "- log2FC", log), cex.main = 3)
  if(!is.na(legend)){
    color.legend(xl, yb, xr, yt, align = "rb", cex = text_size, gradient = "y",
                 legend = plot_legend, colorRampPalette(c("lightgray", "black"))(100))
  }
  
  # Return the modified matrices and the graph object
  return(list(old_mat = intact1, young_mat = intact2, Diff_G = intact))
  #return(G)
}

                              