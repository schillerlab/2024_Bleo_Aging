## As I do not want to mess up my thesis functions, generate a new file and adapt those here
## Original: /Meta_Analysis/integrate_human/NicheNet/NicheNet_functions.R

### Infer Regulatory Potential with NicheNet ###

source("/home/ramadatta/Analysis/1_Schiller_Lab/for_collaborators/Janine/1_Bleomycin_ageing/4_Scripts_Meshal/ligand_receptor_inference.R")

## Get upregulated genes in disease only for receiver + subset to NicheNet ligands
get_upregulated_genes_receiver <- function(condition = "ILD", receiver_ct = "AM"){
    
    ## To deal with cell types with plus symbols in the name
    receiver_ct = gsub("\\+", "\\\\+", gsub(" ", "_", receiver_ct))
    print(receiver_ct)
    
    tab <- de_genes_tab[grep(paste0(receiver_ct, "_up"), rownames(de_genes_tab)), ]
    print(tab)
    genes <- strsplit(tab, ",")[[1]]
    geneset_oi = genes %>% .[. %in% rownames(ligand_target_matrix)]
    print(paste0("Using ", length(geneset_oi), " upregulated genes in ", condition))
    return(geneset_oi)
}


add_expressed_celltypes <- function(tab, lig_act, reg = NA, condition = "ILD", min.pct = 0.1){
    
    #reg = paste0("pct\\.|", paste0("_", condition))
    #tab <- percentages_dot[, grep(condition, colnames(percentages_dot)), drop = F]
    
    reg = ifelse(is.na(reg), paste0("pct\\.|", paste0(condition, "_")), reg)
    mat = data.frame(matrix(0, nrow(lig_act), 1,
                            dimnames = list(lig_act$test_ligand, "expressed_in")), stringsAsFactors = F)

    for(gene in lig_act$test_ligand){
        #cts <- colnames(tab[, tab[gene, ] > min.pct])
        cts <- colnames(tab[which(tab[gene, ] > min.pct)])
        mat[gene, 1] <- paste(gsub(reg, "", cts), collapse = ", ")
    }
    lig_act$expressed_in <- mat$expressed_in
    return(lig_act)
}

add_upregulated_celltypes <- function(genes, lig_act, reg = NA){
    mat = data.frame(matrix(0, nrow(lig_act), 1,
                            dimnames = list(lig_act$test_ligand, "upregulated_in")), stringsAsFactors = F)
    
    for(gene in lig_act$test_ligand){
        cts <- sapply(names(genes), function(x){if (gene %in% genes[[x]]) return(x) else {NA}})
        mat[gene, 1] = paste(names(cts[which(!is.na(cts))]), collapse = ", ")
        }
    if(is.na(reg)){    
        lig_act$upregulated_in <- mat$upregulated_in
    }
    else{
        lig_act$upregulated_in <- gsub(reg, "", mat$upregulated_in)
    }
    return(lig_act)
}

## percentages_meta and percentage_dot can be the same table, if I don’t want to split cell types further
## percentage_meta for meta cell types (receiver)
## percentage_dot for all cell types (sender, with condition specific ct separately)
ligand_activity <- function(percentages_meta, percentages_dot, geneset_oi, receiver_ct = "AM",
                            condition = "ILD", pct.thresh = 0.1, path=path){

    #receiver = paste(receiver_ct, condition, sep = "_")
    expr_genes_receiver = rownames(percentages_meta[percentages_meta[, grep(receiver_ct, colnames(percentages_meta))] > pct.thresh, ])
    #print(expr_genes_receiver)
    ## Use all other cell types as sender (only cells from current disease)
    sender = grep(condition, colnames(percentages_dot), value = T)
    list_expr_genes_sender = lapply(sender, function(x){rownames(percentages_dot[percentages_dot[, x] > pct.thresh, ])})
    expr_genes_sender = list_expr_genes_sender %>% unlist() %>% unique()

    background_expr_genes = expr_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

    ## Define a set of potential ligands
    ligands = lr_network %>% pull(from) %>% unique()
    receptors = lr_network %>% pull(to) %>% unique()
    expr_ligands = intersect(ligands, expr_genes_sender)
    expr_receptors = intersect(receptors, expr_genes_receiver)
    print(paste0("Expressed Ligands ", length(expr_ligands), " Expressed Receptors ", length(expr_receptors)))

    potential_ligands = lr_network %>% filter(from %in% expr_ligands & to %in% expr_receptors) %>%
                        pull(from) %>% unique()
    print(paste0("Potential Ligands ", length(potential_ligands)))
    
    ## Perform NicheNet ligand activity analysis
    ligand_activities = predict_ligand_activities(geneset = geneset_oi,
                                                  background_expressed_genes = background_expr_genes,
                                                  ligand_target_matrix = ligand_target_matrix,
                                                  potential_ligands = potential_ligands)
    
    
    ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
    
    ligand_receptor_inference(best_upstream_ligands=ligand_activities, expressed_receptors=expr_receptors, receiver_ct, path)
    
    return(ligand_activities)

}



get_regulatory_potential <- function(percentages_meta, percentages_dot, geneset_oi, reg = NA, receiver_ct = "AM",
                                     condition = "ILD", save = NA, path = "Data/potential_ligands/", suffix = ""){
    
    ligand_act <- ligand_activity(percentages_meta, percentages_dot, geneset_oi, 
                                  receiver_ct = gsub("\\+", "\\\\+", receiver_ct), 
                                  condition = condition, pct.thresh = 0.1, path=path)
    
    ## Add information in which cell types the genes are expressed / upregulated in disease
    ligand_act <- add_expressed_celltypes(percentages_dot, ligand_act, reg = reg, condition = condition)
    ligand_act <- add_upregulated_celltypes(all_upregulated, ligand_act, reg = reg)
    
    if(!is.na(save)){
        filename <- paste0(path, condition, "_", gsub("/| ", "_", receiver_ct), "_reg_potential", suffix, ".txt")
        print(paste0("Saving to ", filename))
        write.table(ligand_act, filename, quote = F, sep = "\t", row.names = F)
    }
    return(ligand_act)
}

### Plotting ###

scale.func <- switch(EXPR = "radius", 'size' = scale_size, 'radius' = scale_radius,
                     stop("'scale.by' must be either 'size' or 'radius'"))

## unchanged to previous version
## [24.Feb.22] Added a line to colour predicted ligands, that do not pass correlation thresh, differently
plot_ligand_targets <- function(ligand_act, genes_oi, top = 30, condition = "ILD",
                                receiver_ct = "AM", bottom = NA, corr_thresh = NA,path=path){
    
    geneset_oi = genes_oi %>% .[. %in% rownames(ligand_target_matrix)]
    print(paste0("number of genes of interest ",length(geneset_oi)))
      
    if(length(geneset_oi) == 1){
        print("Only one potential Ligand, skipping")
        return(NA)
    }
    
    best_upstream_ligands = ligand_act %>% top_n(top, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
    #print(paste0("number of best_upstream_ligands ",length(best_upstream_ligands)))
    
    active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
    #print(paste0("number of active_ligand_target_links_df ",length(active_ligand_target_links_df)))
    #print(active_ligand_target_links_df)
    
    active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
    #print(paste0("number of active_ligand_target_links ",length(active_ligand_target_links)))
    #print(active_ligand_target_links)
    
    order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
    #print(paste0("order_ligands ", order_ligands))
    
    order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
    #print(paste0("order_targets ", order_targets))
    
    # print("rownames and colnames")
    # print(rownames(active_ligand_target_links))
    # print(colnames(active_ligand_target_links))
    # 
    rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
    colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()
    
    # print(rownames(active_ligand_target_links))
    # print(colnames(active_ligand_target_links))
    # print(class(active_ligand_target_links))
    # print(active_ligand_target_links)
    #print(paste0("active_ligand_target_links", active_ligand_target_links))
    
     #vis_ligand_target = active_ligand_target_links[order_targets, order_ligands] #%>% t()
     #print("vis_ligand_target without transpose")
     #print(vis_ligand_target)
    #vis_ligand_target = active_ligand_target_links %>% t()
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
    #print("vis_ligand_target with transpose")
    #print(vis_ligand_target)
    
   #results_dir="/home/ramadatta/Analysis/1_Schiller_Lab/for_collaborators/Janine/1_Bleomycin_ageing/Bleomycin_datta_output_scripts/2_Targeted_diffxpy_dge_regpotential_GzmkT_AT2/Nichenet_output/"
   write.table(ligand_act, paste0(path, '/','nichenet_ligand_activity_', receiver_ct,'_',  condition, '.txt'), sep='\t')
   write.table(vis_ligand_target, paste0(path, '/','nichenet_ligand_target_matrix_', receiver_ct,'_', condition,  '.txt', sep='\t'))
    
    
    p_ligand_target_network = vis_ligand_target %>%
            make_heatmap_ggplot("Prioritized ligands",
                                paste0("Predicted target genes (", receiver_ct, " ", condition, ")"),
                                legend_position = "top", x_axis_position = "top",
                                legend_title = "Regulatory potential")  + 
            theme(axis.text = element_text(size = 18), title = element_text(size = 20), 
                  axis.text.x = element_text(colour = "black"),
                  plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 18)) +
            scale_fill_gradient2(low = "whitesmoke", high = condition_colors[paste0(condition, "_vertex")], 
                                 breaks = c(0, 0.0045, 0.0090))
    
    if(!is.na(bottom)){
        p_ligand_target_network <- p_ligand_target_network + theme(plot.margin = unit(c(1, 1, bottom, 1.2), "cm"))
    }
    
    if(!is.na(corr_thresh)){
        genes = unique(as.character(p_ligand_target_network$data$y))
        highlight <- ifelse(ligand_act$pearson > corr_thresh, "black", "gray60")[1:length(genes)]
        p_ligand_target_network <- p_ligand_target_network + theme(axis.text.y = element_text(colour = rev(highlight)))
    }
    
    width = round(length(unique(p_ligand_target_network$data$x)) / 3)
    if(width < 20){
        p_ligand_target_network <- p_ligand_target_network + theme(plot.margin = unit(c(1, 20 - width, bottom, 1.2), "cm"))
    }  
    return(p_ligand_target_network)
}

## Re-write to not include control anymore, won’t be possible to make the split dotplots
ligand_dotplot <- function(percentages_dot, avg_expr_dot, lig_act, genes, order, condition = "ILD", colour = "ILD",
                           scale.min = NA, reg = NA, scale.max = NA, star = NA, col.min = -2.5, col.max = 2.5, corr_thresh = NA,
                           scale = TRUE, dot.min = 0.1, dot.scale = 10, top = 30, top_pad = 5.5, title_add = ""){
    
    data.plot = t(cbind(percentages_dot[genes, ], avg_expr_dot[genes, ]))
    
    #reg = paste0("pct\\.|_", baseline, "|", paste0("_", condition))
    reg = ifelse(is.na(reg), paste0("pct\\.|", paste0(condition, "_")), reg)
    data.plot <- cbind(data.plot, data.frame(type = str_match(rownames(data.plot), "avgExpr|pct"),
                                             cell_type = gsub(reg, "", rownames(data.plot)), 
                                             #condition = str_match(rownames(data.plot), paste0(baseline, "|", condition)),
                                             row.names = rownames(data.plot)))    
    data.plot <- reshape2::melt(data.plot, id.vars = c("type", "cell_type"))
    
    colnames(data.plot) <- c("type", "cell_type", "gene", "value")
    data.plot$cell_type <- gsub("avgExpr_", "", data.plot$cell_type) ## don’t know why this happens during melting
    
    data.plot = reshape2::dcast(data.plot, cell_type + gene ~ type, value.var = "value", fun.aggregate = sum)
    data.plot$cell_type <- factor(data.plot$cell_type, levels = order)    
    #nudge <- ifelse(is.na(split), 0.15, 0)          
    
    p <- ggplot(data = data.plot, mapping = aes(x = cell_type, y = gene)) +
            geom_point(data = data.plot, mapping = aes(size = pct, color = avgExpr)) + 
            # position = position_nudge(x = nudge, y = 0)) +
            scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
            scale_colour_gradient(low = condition_colors[paste0(colour, "_mid")],
                                  high = condition_colors[colour]) + #, midpoint = 0.05) +
            theme_minimal() + labs(title = paste("Expression in", condition, title_add)) +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                  title = element_text(size = 20), plot.title = element_text(hjust = 0.5),
                  axis.text = element_text(size = 18), legend.text = element_text(size = 18),
                  panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
                  plot.margin = unit(c(top_pad, 1, 1.5, 2), "cm")) +
            guides(size = guide_legend(title = 'pct expressed'))
    
    #if(!is.na(split)){
        ## Add healthy values to the right (separate, as else it screws up my ordering)
    #    p <- p + geom_point(data = data.plot[data.plot$condition == baseline, ], shape = 21, 
    #                      mapping = aes(size = pct, color = avgExpr),
    #                      position = position_nudge(x = -0.15, y = 0)) + scale_x_discrete(limits = order)
    #    }
       
    if(!is.na(corr_thresh)){
        highlight <- ifelse(lig_act$pearson > corr_thresh, "black", "gray60")[1:length(genes)]
        p <- p + theme(axis.text.y = element_text(colour = rev(highlight)))
    }
    
    if((!is.na(star)) & (nrow(data.plot) > 0)){
        data.plot$highlight <- 0
        #highlight <- rev(lig_act$upregulated_in[1:top])
        highlight <- lig_act$upregulated_in[match(genes, lig_act$test_ligand)]
        
        names(highlight) <- genes
        if(length(genes) > 0){
            for(gene in genes){
                high <- strsplit(highlight[gene], ", ")[[1]]
                data.plot[(data.plot$gene == gene) & (data.plot$cell_type %in% high), "highlight"] <- 1
            }
        }
        p <- p + geom_point(data = data.plot[data.plot$highlight > 0, ], color = "black", pch = 8,
                            position = position_nudge(x = 0.07, y = 0.45), size = 3, stroke = 0.75)
    }
    return(p)
}


plot_regulatory_potential <- function(avg_pct, avg_expr, ligand_act, order, reg = NA, receiver_ct = "AM", condition = "ILD",
                                      colour = "ILD", star = T, do_return = NA, dot.min = 0.1, corr_thresh = NA,
                                      title_add = "", top = 30, top_pad = 1, bottom = 1, width = 30,
                                      show = NA, path = "Plots/potential_ligands/", suffix = "", save = NA){

    p1 <- plot_ligand_targets(ligand_act, genes = geneset_oi, top = top, condition = colour, bottom = bottom,
                              receiver_ct = receiver_ct, corr_thresh = corr_thresh, path=path)
    
    ## Only plot if there are more than 1 potential Ligands
    if(is.na(p1[1])){return(NA)}
    p2 <- ligand_dotplot(avg_pct, avg_expr, ligand_act, genes = unique(as.character(p1$data$y)), top_pad = top_pad,
                         dot.min = dot.min, top = top, order = order, reg = reg, corr_thresh = corr_thresh,
                         condition = condition, colour = colour, star = star, title_add = title_add)
    
    if(!is.na(show)){
        grid.arrange(p1, p2, nrow = 2)
    }

    if(!is.na(save)){
        height <- round(length(unique(p1$data$y)) / 2.5)
        height <- ifelse(height < 10, 10, height)

        filename <- paste0(path, condition, "_", gsub("/| ", "_", receiver_ct), "_reg_potential", suffix, ".pdf")
        #width = round(length(unique(p1$data$x)) / 3) + 12
        ## fixed width instead
        height = height 
        g <- arrangeGrob(p1, p2, nrow = 2)
        print(paste0("Saving to ", filename))
        ggsave(g, file = filename, width = width, height = 2 * height)
    }
    if(!is.na(do_return)){
        return(list(reg_potential = p1, dotplot = p2))
    }
}







