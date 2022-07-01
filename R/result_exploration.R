
#' Generate summary of specific pathways and top contributing genes given an annotation
#'
#'
#' @param PCA_object List, output from run_activity_analysis().
#' @param file Character, path to file if you want the analysis stored as xlsx file with one sheet by cluster.
#' @param meta Dataframe containing metadata, with cell labels as rownames.
#' @param annot_name Character, name of the column in meta that should be used as annotation for cells.
#' @param noise Double, value between 0 and 1 indicating the maximum specificity you allow a non-specific cluster to have. Can be increased if no pathways appear as specific.
#' @param top_genes Integer, number of top genes contributing to PC to print.
#' @param min_contrib Double, minimum contribution/weight of a gene to the PC. It overrides top_genes parameter and only top_genes contributing above min_contrib will be printed.
#'
#' @return List, one assay by cluster defined by annotation, and each displaying a table with the pathway, the top contributing genes and the pathway specificity.
#'
#'
#' @export
generate_summary_specific_modules<-function(PCA_object,file=NULL,meta,annot_name,top_genes=10,min_contrib=0.1){

    #### Check parameters ####
    stopifnot(is.list(PCA_object),
              is.data.frame(meta),is.character(annot_name),
              is.numeric(top_genes),
              is.numeric(min_contrib))

    if(!is.null(annot_name)&!(annot_name %in% colnames(meta))){
        stop("Uncorrect annot_name provided")
    }

    #### Build activity matrix and specificity table ####
    activity_mat<-build_activity_mat(PCA_object)
    specif<-specificity_table(activity_mat,meta = meta, annot_name =annot_name)
    colnames(specif)<-sapply(colnames(specif), function(y) strsplit(y, "S_")[[1]][2])

    #### Find specificity threshold ####
    nb_clusters<-length(unique(meta[,annot_name]))
    thr<-1/nb_clusters

    #### Find gene sets specific to each cluster ####
    active_pathways<-lapply(colnames(specif), function(x){
        rownames(specif)[which(specif[,x]>thr)]
    })
    names(active_pathways)<-colnames(specif)

    # specificity of each specific activation mode
    list_specif<-lapply(colnames(specif), function(x){
        specif[which(specif[,x]>thr),x]
    })
    names(list_specif)<-colnames(specif)

    #### Generate final output with top contributing genes ####
    list_summary<-lapply(names(active_pathways),function(y){
        # list of top N genes for each cluster
        genes_by_path<-lapply(active_pathways[[y]], function(x){
            pathway<-strsplit(x, "_mode[1-9]")[[1]][1]
            #pc<-as.numeric(tail(strsplit(x, "_PC")[[1]],n=1))
            pc<-tail(strsplit(x, "_mode")[[1]],n=1)
            contrib<-t(PCA_object[[pathway]]$gene_contrib)
            n<-min(top_genes,nrow(contrib))
            if(ncol(contrib)==1){
                contrib<-contrib[order(contrib,decreasing = T),]
                names<-names(contrib)
                m<-min(n,length(which(contrib>min_contrib)))
                if(m>0){
                    c(names[1:m],rep("",times=top_genes-m))
                }else{
                    rep("",times=top_genes)
                }
            }else{
                contrib<-contrib[order(contrib[,paste0("PC",pc)],decreasing = T),]
                names<-names(contrib[,paste0("PC",pc)])
                m<-min(n,length(which(contrib[,paste0("PC",pc)]>min_contrib)))
                if(m>0){
                    c(names[1:m],rep("",times=top_genes-m))
                }else{
                    rep("",times=top_genes)
                }
            }
        })
        names(genes_by_path)<-active_pathways[[y]]
        # transform into dataframe
        #aggregate_mat<-t(as.matrix(bind_rows(genes_by_path)))
        aggregate_mat<-cbind(t(as.matrix(dplyr::bind_rows(genes_by_path))),list_specif[[y]])
        if(length(list_specif[[y]])>1){aggregate_mat<-aggregate_mat[order(list_specif[[y]],decreasing = T),]}
        aggregate_mat
    })
    names(list_summary)<-names(active_pathways)
    list_summary<-list_summary[lapply(list_summary,length)!=0]

    #### Write to file if specified ####
    if(!is.null(file)){
        if (!requireNamespace("xlsx", quietly = TRUE)) {
            stop("Package \"xlsx\" needed for this function to work. Please install it.",
                 call. = FALSE)
        }else{
            for(i in 1:length(list_summary)){
                xlsx::write.xlsx(list_summary[[i]],file=file,sheetName = names(list_summary)[i],row.names = T,col.names = F,append = T)
            }
        }
    }else{
        list_summary
    }

}



#' Plots UMAP colored by cell annotation
#'
#' @param umap UMAP object
#' @param type Description of the labels that are plotted (cell type, unsupervised clusters,...)
#' @param labels Character vector, containing annotation for cells, in the same order as cells in the umap object
#' @param title Title of the plot.
#' @param colors List of colors, else default custom palette colors will be used
#'
#' @return UMAP plot
#'
#' @import ggplot2
#'
#' @export
plot_umap_annot<-function(umap,type="Cell type",labels,title=NULL,colors=NULL){

    #### Build plot ####
    df <- data.frame(x = umap$layout[,1],
                     y = umap$layout[,2],
                     Alias = labels)

    if(is.null(colors)){
        palette=generate_palette(length(unique(labels)))
        p<-ggplot(df, aes(x, y, colour = Alias)) +
            geom_point(shape=21,colour="white",aes(fill=Alias),size=2,stroke=0.5)+
            xlab("UMAP_1")+ylab("UMAP_2")+
            ggtitle(title)+theme_classic()+
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_fill_manual(values=palette,breaks=unique(labels)) +
            labs(fill = type)
    }
    else{
        p<-ggplot(df, aes(x, y, colour = Alias)) +
            geom_point(shape=21,colour="white",aes(fill=Alias),size=2,stroke=0.5)+
            xlab("UMAP_1")+ylab("UMAP_2")+
            ggtitle(title)+theme_classic()+
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_fill_manual(values=colors,breaks=unique(labels)) +
            labs(fill = type)
    }
    print(p)
}



#' Plots UMAP colored by expression level of a gene
#'
#' @param umap UMAP object
#' @param expr_mat Numeric matrix, containing counts, normalized or not.
#' @param gene Character, gene for which we want to plot expression
#'
#' @return UMAP plot
#'
#' @import ggplot2
#'
#' @export
plot_umap_gene<-function(umap,expr_mat,gene){

    #### Check parameters ####
    stopifnot(is.character(gene))

    if(!(gene %in% rownames(expr_mat))){
        stop(paste0("No expression available for ",gene))
    }

    #### Build plot ####
    df <- data.frame(x = umap$layout[,1],
                     y = umap$layout[,2],
                     gene_expr = expr_mat[gene,])
    p<-ggplot(df, aes(x, y, colour = gene_expr)) +
        geom_point(shape=21,colour="white",aes(fill=gene_expr),size=2)+
        xlab("UMAP_1")+ylab("UMAP_2")+
        ggtitle(gene)+
        theme_classic()+
        scale_fill_viridis_c()+
        theme(plot.title = element_text(hjust = 0.5))
    print(p)
}


#' Plots UMAP colored by activity level of a module
#'
#' @param umap UMAP object
#' @param PCA_object List, output from run_activity_analysis().
#' @param module Character, module for which you want to plot activity
#' @param scaled Boolean, by default uses the activity scores scaled between 0 and 1, set to FALSE if you want to use raw PCA scores.
#'
#' @return UMAP plots (one by activation mode)
#'
#' @import ggplot2
#'
#' @export
plot_umap_pathway_activity<-function(umap,PCA_object,module,scaled=T){

    #### Check parameters ####
    stopifnot(is.list(PCA_object),is.character(module),is.logical(scaled))

    if(!(module %in% names(PCA_object))){
        stop(paste0("No informative activity available for ",module))
    }

    #### Build plot ####
    # number of components for this module
    nb_comp<-length(PCA_object[[module]]$expl_var)
    # get activity matrix
    if(scaled){
        activity_mat<-PCA_object[[module]]$activity_scores
    }
    else{
        activity_mat<-PCA_object[[module]]$activity_scores_raw
    }

    # get gene contribution
    contrib<-t(PCA_object[[module]]$gene_contrib)
    if(nb_comp==1){
        comp<-strsplit(rownames(activity_mat), "mode")[[1]][2]
        # build a vector of top 3 genes for the unique activation mode
        contrib<-contrib[order(contrib,decreasing = T),]
        names<-names(contrib)
        top_genes<-paste(comp,names[1],names[2],names[3],sep="_")
    }else{
        # build a list of top 3 genes for each component
        comp<-sapply(colnames(contrib), function(x) strsplit(x, "PC")[[1]][2])
        list_top_genes<-lapply(comp,function(x) {
            contrib<-contrib[order(contrib[,paste0("PC",x)],decreasing = T),]
            names<-names(contrib[,paste0("PC",x)])
            paste(x,names[1],names[2],names[3],sep="_")
        })
        top_genes<-unlist(list_top_genes)
    }
    # plot
    df<-data.frame(x=rep(umap$layout[,1],nb_comp),y=rep(umap$layout[,2],nb_comp),Activity = as.vector(t(activity_mat)),PC=sort(rep(top_genes,times=ncol(activity_mat))))
    p<-ggplot(df, aes(x, y, colour =Activity)) +
        geom_point(shape=21,colour="white",aes(fill=Activity),size=2)+
        xlab("UMAP_1")+ylab("UMAP_2")+
        ggtitle(module)+theme_bw()+
        scale_fill_viridis_c() +
        facet_wrap(~PC)+ theme(plot.title = element_text(hjust = 0.5))
    print(p)
}

#' Plots heatmap of top contributing genes for each mode
#'
#' @param expr_mat Numeric matrix, containing counts, normalized or not.
#' @param PCA_object List, output from run_activity_analysis().
#' @param module Character, module for which you want to plot activity
#' @param n Number of top contributors to plot for each mode
#' @param meta Dataframe containing metadata, with cell labels as rownames.
#' @param annot_name Character, name of the column in meta that should be used to compute homogeneity for.
#' @param cluster_cols Boolean, set to TRUE to perform column clustering.
#' @param fontsize Integer, base fontsize for the plot.
#' @param colors_annot Vector of colors to use.
#'
#' @return heatmap plot of top contributor genes for each informative mode of activation, ordered by decreasing contribution
#'
#' @import ggplot2
#'
#' @export
plot_heatmap_pathway_top_contrib_genes<-function(expr_mat,PCA_object,module,n=10,meta=NULL,annot_name=NULL,cluster_cols=T,fontsize=7,colors_annot=NULL){

    #### Check parameters ####
    stopifnot(is.list(PCA_object),is.character(module),is.numeric(n),n>1)

    if(!(module %in% names(PCA_object))){
        stop(paste0("No informative activity available for ",module))
    }

    #### Build plot ####
    # number of components for this module
    nb_comp<-length(PCA_object[[module]]$expl_var)

    top_genes<-c()

    # get gene contribution
    contrib<-t(PCA_object[[module]]$gene_contrib)
    if(nb_comp==1){
        contrib<-contrib[order(contrib,decreasing = T),]
        names<-names(contrib)
        top_genes<-names[1:n]
    }else{
        # build a list of top 3 genes for each component
        comp<-sapply(colnames(contrib), function(x) strsplit(x, "PC")[[1]][2])
        list_top_genes<-lapply(comp,function(x) {
            contrib<-contrib[order(contrib[,paste0("PC",x)],decreasing = T),]
            names<-names(contrib[,paste0("PC",x)])
            names[1:n]
        })
        top_genes<-unlist(list_top_genes)
    }
    # plot

    if(!is.null(annot_name)){
        mat_col <- data.frame(as.factor(meta[,annot_name]))
        colnames(mat_col)<-annot_name
        rownames(mat_col) <- colnames(expr_mat)

        mat_colors <- list()
        if(!is.null(colors_annot)){
            if(length(unique(meta[,annot_name]))==length(colors_annot)){
                mat_colors[[annot_name]]<-colors_annot
                names(mat_colors[[annot_name]]) <- unique(meta[,annot_name])
            }
            else{stop(paste0("Number or colors provided in colors_annot different from number of levels of ",annot_name))}

        }
        else{
            palette=generate_palette(length(unique(meta[,annot_name])))
            names(palette)=unique(meta[,annot_name])
            mat_colors<-list(palette)
            names(mat_colors)<-annot_name
        }

    }

    if(nb_comp>1){
        gaps_row=seq(from=n, to=n*nb_comp-n,by=n)
    }
    else{gaps_row=n}

    mat=Matrix::t(scale(Matrix::t(expr_mat[top_genes,]),center=T,scale=T))

    if(!is.null(annot_name)){
        p<-pheatmap(
            mat               = mat,
            color             = inferno(10,direction=1),
            border_color      = NA,
            show_colnames     = FALSE,
            show_rownames     = TRUE,
            annotation_col    = mat_col,
            annotation_colors = mat_colors,
            drop_levels       = TRUE,
            fontsize          = fontsize,
            main              = paste0("Top contributor genes - ",module),
            cluster_rows = F,
            cluster_cols = cluster_cols,
            gaps_row = gaps_row,
            clustering_distance_cols = "euclidean",
            clustering_method = "ward.D2"

        )

    }
    else{
        p<-pheatmap(
            mat               = mat,
            color             = inferno(10,direction=1),
            border_color      = NA,
            show_colnames     = FALSE,
            show_rownames     = TRUE,
            drop_levels       = TRUE,
            fontsize          = fontsize,
            main              = paste0("Top contributor genes - ",module),
            cluster_rows = F,
            cluster_cols = cluster_cols,
            gaps_row = gaps_row,
            clustering_distance_cols = "euclidean",
            clustering_method = "ward.D2"

        )

    }


}



#' Plots heatmap showing individual gene expression of a pathway
#'
#' @param expr_mat Numeric matrix, containing counts, normalized or not.
#' @param activity_mat Numeric matrix, containing activity for each activation mode of each pathway.
#' @param meta Dataframe containing metadata, with cell labels as rownames.
#' @param annot_name Character, name of the column in meta that should be used to compute homogeneity for.
#' @param modules_list List, with module names as list names and character vectors containing genes in module. Output from read_gmt file for example.
#' @param module Character, pathway for which you want to plot expression.
#' @param cluster_cols Boolean, set to TRUE to perform column clustering.
#' @param cluster_rows Boolean, set to TRUE to perform row clustering.
#' @param clustering_distance Character, distance to use to perform clustering, one of "euclidean","correlation".
#' @param clustering_method Character, one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param fontsize Integer, base fontsize for the plot.
#' @param colors_annot Vector of colors to use.
#'
#' @return heatmap plot
#'
#' @import pheatmap
#'
#' @export
plot_heatmap_pathway_genes<-function(expr_mat,activity_mat=NULL,meta=NULL,annot_name=NULL,modules_list,module,cluster_cols=T,cluster_rows=T,clustering_distance="euclidean",clustering_method="ward.D2",fontsize=5,colors_annot=NULL){

    #### Check parameters ####
    stopifnot(is.list(modules_list),is.character(module),
              is.logical(cluster_cols),is.logical(cluster_rows),
              clustering_distance %in% c("euclidean","correlation"),
              clustering_method %in% c("ward.D", "ward.D2", "single", "complete", "average","mcquitty","median","centroid"),
              is.numeric(fontsize))

    if(!is.null(annot_name)){
        if(!(annot_name %in% colnames(meta))){
            stop("Uncorrect annot_name provided")
        }
    }

    #### Build annotation and colors ####
    if(!is.null(activity_mat)){
        modules<-unlist(lapply(rownames(activity_mat), function(x) strsplit(x, "_mode")[[1]][1]))
        if(!(module %in% modules)){
            stop(paste0("No informative activity available for ",module," - re-run without providing the activity matrix"))
        }
        nb_comp<-length(which(modules==module))
        sub_activity_mat<-activity_mat[which(modules==module),]
        if(nb_comp==1){
            comp_name<-rownames(activity_mat)[which(modules==module)]
            annot_to_plot<-data.frame(sub_activity_mat)
            colnames(annot_to_plot)<-paste("Mode",strsplit(comp_name, "_mode")[[1]][2],sep="_")
            annotation_colors<-list(viridis::viridis(5))
            names(annotation_colors)<-paste("Mode",strsplit(comp_name, "_mode")[[1]][2],sep="_")
        }else{
            annot_to_plot<-data.frame(t(sub_activity_mat))
            comp<-unlist(lapply(colnames(annot_to_plot), function(x) strsplit(x, "_mode")[[1]][2]))
            colnames(annot_to_plot)<-paste("Mode",comp,sep="_")
            annotation_colors<-lapply(comp,function(x) viridis::viridis(5))
            names(annotation_colors)<-paste("Mode",comp,sep="_")
        }
    }

    if(!is.null(annot_name)){
        mat_col <- data.frame(as.factor(meta[,annot_name]))
        colnames(mat_col)<-annot_name
        rownames(mat_col) <- colnames(expr_mat)

        mat_colors <- list()
        if(!is.null(colors_annot)){
            mat_colors[[annot_name]]<-colors_annot
            names(mat_colors[[annot_name]]) <- unique(meta[,annot_name])
        }
        else{
            # mat_colors<-list(heatmap_annot_color_categorical(meta[,annot_name]))
            # names(mat_colors)<-annot_name
            palette=generate_palette(length(unique(meta[,annot_name])))
            names(palette)=unique(meta[,annot_name])
            mat_colors<-list(palette)
            names(mat_colors)<-annot_name
        }

    }

    #### Build plot ####
    if(!is.null(annot_name) & is.null(activity_mat)){
        pheatmap(
            mat=expr_mat[intersect(rownames(expr_mat),modules_list[[module]]),],
            color             = inferno(10,direction=1),
            border_color      = NA,
            annotation_col = mat_col,
            cluster_cols = cluster_cols,
            cluster_rows=cluster_rows,
            fontsize = fontsize,
            clustering_distance_rows = clustering_distance,
            clustering_distance_cols = clustering_distance,
            clustering_method=clustering_method,
            main=module,
            show_colnames = F,
            annotation_colors = mat_colors,
            drop_levels       = TRUE)
    }

    if(is.null(annot_name) & !is.null(activity_mat)){
        pheatmap(
            mat=expr_mat[intersect(rownames(expr_mat),modules_list[[module]]),],
            color             = inferno(10,direction=1),
            border_color      = NA,
            annotation_col = annot_to_plot,
            cluster_cols = cluster_cols,
            cluster_rows=cluster_rows,
            fontsize = fontsize,
            clustering_distance_rows = clustering_distance,
            clustering_distance_cols = clustering_distance,
            clustering_method=clustering_method,
            main=module,
            show_colnames = F,
            annotation_colors = annotation_colors,
            drop_levels       = TRUE)
    }

    if(!is.null(annot_name) & !is.null(activity_mat)){
        pheatmap(
            mat=expr_mat[intersect(rownames(expr_mat),modules_list[[module]]),],
            color             = inferno(10,direction=1),
            border_color      = NA,
            annotation_col = cbind(mat_col,annot_to_plot),
            cluster_cols = cluster_cols,
            cluster_rows=cluster_rows,
            fontsize = fontsize,
            clustering_distance_rows = clustering_distance,
            clustering_distance_cols = clustering_distance,
            clustering_method=clustering_method,
            main=module,
            show_colnames = F,
            annotation_colors = c(mat_colors,annotation_colors),
            drop_levels       = TRUE)
    }

    if(is.null(annot_name) & is.null(activity_mat)){
        pheatmap(
            mat=expr_mat[intersect(rownames(expr_mat),modules_list[[module]]),],
            color             = inferno(10,direction=1),
            border_color      = NA,
            cluster_cols = cluster_cols,
            cluster_rows=cluster_rows,
            fontsize = fontsize,
            clustering_distance_rows = clustering_distance,
            clustering_distance_cols = clustering_distance,
            clustering_method=clustering_method,
            main=module,
            show_colnames = F,
            drop_levels       = TRUE)
    }

}


#' Plots heatmap showing activity scores for each activation mode
#'
#' @param activity_mat Numeric matrix, containing activity for each activation mode of each pathway
#' @param meta Dataframe containing metadata, with cell labels as rownames.
#' @param annot_name Character, name of the column in meta that should be used to compute homogeneity for.
#' @param cluster_cols Boolean, set to TRUE to perform column clustering.
#' @param cluster_rows Boolean, set to TRUE to perform row clustering.
#' @param clustering_distance Character, distance to use to perform clustering, one of "euclidean","correlation".
#' @param clustering_method Character, one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param fontsize Integer, base fontsize for the plot.
#' @param colors_annot Vector of colors to use.
#'
#' @return heatmap plot
#'
#' @import pheatmap
#'
#' @export
plot_heatmap_activity_mat<-function(activity_mat,meta=NULL,annot_name=NULL,cluster_cols=T,cluster_rows=T,clustering_distance="euclidean",clustering_method="ward.D2",fontsize=5,colors_annot=NULL,fontsize_row=5,fontsize_col=5){

    #### Check parameters ####
    stopifnot(is.numeric(activity_mat),
              is.logical(cluster_cols),is.logical(cluster_rows),
              clustering_distance %in% c("euclidean","correlation"),
              clustering_method %in% c("ward.D", "ward.D2", "single", "complete", "average","mcquitty","median","centroid"),
              is.numeric(fontsize))

    if(!is.null(annot_name)){
        if(!(annot_name %in% colnames(meta))){
            stop("Uncorrect annot_name provided")
        }
    }

    #### Build plot ####
    if(!is.null(annot_name)){
        mat_col <- data.frame(as.factor(meta[,annot_name]))
        colnames(mat_col)<-annot_name
        rownames(mat_col) <- colnames(activity_mat)

        mat_colors <- list()
        if(!is.null(colors_annot)){
            mat_colors[[annot_name]]<-colors_annot
            names(mat_colors[[annot_name]]) <- unique(meta[,annot_name])
        }
        else{
            # mat_colors<-list(heatmap_annot_color_categorical(meta[,annot_name]))
            # names(mat_colors)<-annot_name
            palette=generate_palette(length(unique(meta[,annot_name])))
            names(palette)=unique(meta[,annot_name])
            mat_colors<-list(palette)
            names(mat_colors)<-annot_name
        }
        pheatmap(
            activity_mat,
            annotation_col = mat_col,
            cluster_cols = cluster_cols,
            cluster_rows=cluster_rows,
            fontsize = fontsize,
            clustering_distance_rows = clustering_distance,
            clustering_distance_cols = clustering_distance,
            clustering_method=clustering_method,
            show_colnames = F,
            annotation_colors = mat_colors,
            color             = viridis(10,direction=1),
            main= "Activity matrix",
            fontsize_row = fontsize_row,
            fontsize_col = fontsize_col)
    }
    else{
        pheatmap(activity_mat,
                 cluster_cols = cluster_cols,
                 cluster_rows=cluster_rows,
                 fontsize = fontsize,
                 clustering_distance_rows = clustering_distance,
                 clustering_distance_cols = clustering_distance,
                 clustering_method=clustering_method,
                 show_colnames = F,
                 color = viridis(10,direction=1),
                 main= "Activity matrix",
                 fontsize_row = fontsize_row,
                 fontsize_col = fontsize_col)
    }
}


#' Plots barplot colored by cluster showing distribution of specificity between clusters for each activation mode
#'
#' @param PCA_object List, output from run_activity_analysis().
#' @param module Character, module for which we want to plot activity.
#' @param meta Dataframe containing metadata, with cell labels as rownames.
#' @param annot_name Character, name of the column in meta that should be used to compute homogeneity for.
#' @param colors_annot Vector of colors to use
#'
#' @return UMAP plots (one by activation mode)
#'
#' @import ggplot2
#'
#' @export
plot_pathway_specificity<-function(PCA_object,module,meta,annot_name,colors_annot=NULL){

    #### Check parameters ####
    stopifnot(is.list(PCA_object),is.character(module),
              is.data.frame(meta),is.character(annot_name))

    if(!(module %in% names(PCA_object))){
        stop(paste0("No informative activity available for ",module))
    }
    if(!(annot_name %in% colnames(meta))){
        stop("Uncorrect annot_name provided")
    }

    #### Build plot ####
    # number of components for this module
    nb_comp<-length(PCA_object[[module]]$expl_var)
    # get activity matrix (scaled, otherwise specificity doesn't work)
    activity_mat<-PCA_object[[module]]$activity_scores

    # compute specificity
    specif<-as.matrix(specificity_table(activity_mat = as.data.frame(activity_mat),meta = meta,annot_name = annot_name))
    colnames(specif)<-sapply(colnames(specif), function(y) strsplit(y, "S_")[[1]][2])
    comp<-sapply(rownames(specif), function(x) strsplit(x, "mode")[[1]][2])
    df<-data.frame(specificity=as.vector(t(as.matrix(specif))),Cluster=rep(colnames(specif),nb_comp),PC=sort(rep(as.numeric(comp),times=ncol(specif))))
    df <- within(df, PC <- factor(PC, levels=sort(unique(PC),decreasing=T)))

    if(!is.null(colors_annot)){
        p<-ggplot(df,aes(x=PC,y=specificity,fill=Cluster))+
            geom_bar(stat="identity")+
            coord_flip()+
            xlab("Activation mode")+
            ylab("Specificity")+
            ggtitle(module)+
            theme_bw()+
            theme(plot.title = element_text(hjust = 0.5))+
            scale_fill_manual(values=colors_annot,breaks=unique(meta[,annot_name]))
    }

    else{
        palette=generate_palette(length(unique(meta[,annot_name])))
        p<-ggplot(df,aes(x=PC,y=specificity,fill=Cluster))+
            geom_bar(stat="identity")+
            coord_flip()+
            xlab("Activation mode")+
            ylab("Specificity")+
            ggtitle(module)+
            theme_bw()+
            theme(plot.title = element_text(hjust = 0.5))+
            scale_fill_manual(values=palette,breaks=unique(meta[,annot_name]))
    }

    print(p)
}










