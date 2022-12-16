
#' Annotate individual cells
#'
#' @param Annotation Vector, containing cluster annotation for individual cells
#' @param clusters_annot Dataframe giving corresponding annotation for each cluster
#'
#' @return Vector containing annotation for individual cells
#'
#' @import ggplot2
#'
#' @export
annotate_cells<-function(annotation,cluster_annot){

    # merge dataframe with clusters for each cell with correspondance table between clusters and cell types, and return a result with cells in the same order as in the annotation vector.
    df<-data.frame(cell_cluster=annotation,order=1:length(annotation))
    df<-merge(df,cluster_annot,by.x="cell_cluster",by.y="cluster",all.x=T)
    df<-df[which(!is.na(df$order)),]
    df<-df[order(as.numeric(df$order)),]

    return(df$cell_type)
}

#' Generate cell type annotation for each cell and provide complementary results to interpret results
#'
#' @param PCA_obj Vector, containing cluster annotation for individual cells
#' @param thr Minimum average activity by cluster required to be assigned an annotation. Default 0, increase to be more stringent on assignation confidence
#'
#' @return List containing cell type annotation, matrix of average score by Leiden cluster and cell type used for attribution, cell annotation with Leiden clusters, activity matrix and UMAP computed on activity matrix.
#'
#' @import RANN
#' @import igraph
#' @import leidenbase
#'
#' @export
generate_cell_type_annotation<-function(PCA_obj,thr=0,compute_umap=T){
    # compute activity mat
    activity_mat<-build_activity_mat(PCA_obj,scaled = F)
    # clustering
    knn.matrix <- RANN::nn2(t(activity_mat), t(activity_mat), k = 20, searchtype = "standard")[[1]]

    jaccard.adj <- knn_jaccard(knn.matrix)
    set.seed(2016)
    graph<-igraph::graph.adjacency(
        jaccard.adj,
        mode = "undirected",
        weighted = TRUE
    )
    set.seed(2016)
    cluster_result <- leidenbase::leiden_find_partition(
        graph,
        verbose = FALSE,
        seed=2016,num_iter = 2, partition_type="ModularityVertexPartition",edge_weights=E(graph)$weight
    )
    # unknown threshold
    cluster_result$membership<-paste0("C",cluster_result$membership)
    tmp<-average_by_cluster(activity_mat,cluster_result$membership)
    rownames(tmp)<-sapply(rownames(tmp),function(x) strsplit(x,"_mode")[[1]][1])
    thr=thr
    attrib_df<-data.frame(cluster=colnames(tmp),cell_type=apply(tmp,2,function(x) ifelse(max(x)>thr,rownames(tmp)[which.max(x)],"unassigned")))
    attrib<-annotate_cells(cluster_result$membership,attrib_df)

    colnames(tmp)<-sapply(colnames(tmp),function(x) paste0(x,"_",attrib_df[which(attrib_df$cluster==x),"cell_type"]))

    if(compute_umap){
        message("Computing UMAP on activity matrix")
        umap<-run_umap(activity_mat)
    }
    else{umap<-NULL}

    return(list(cell_annotation=attrib,cluster_matrix=tmp,clusters_annotation=paste0(cluster_result$membership,"_",attrib),activity_matrix=activity_mat,umap=umap))
}

#' Run MAYA for unsupervised cell annotation with cell type
#'
#' @param expr_mat Raw count matrix. If normalized, set is_logcpm to T
#' @param modules_list List with cell types associated with their markers. If NULL, PanglaoDB will be loaded.
#' @param min_cells_pct Numeric between 0 and 1, minimum percentage of cells that should be above informativity threshold to consider activity score interesting.
#' @param organs Character to specify if list for a specific organ should be loaded. One of ("Pancreas","Brain","Lungs","Heart","Liver","Adrenal glands","GI tract","Reproductive","Kidney","Zygote","Thyroid","Embryo","Bone","Skin","Mammary gland","Eye","Olfactory system","Oral cavity","Thymus","Placenta","Urinary bladder"). "All" is another possible value, the 147 PanglaoDB cell types will be loaded.
#' @param is_logcpm Set to TRUE if data already normalized
#' @param nCores Number of cores to use. Set to 1 if working in a Windows environment, otherwise you can use the function detectCores() to find out how many cores are available.
#' @param thr Minimum average activity by cluster required to be assigned an annotation. Default 0, increase to be more stringent on assignation confidence.
#' @param max_contrib Numeric between 0 and 1, representing the maximum contribution to a mode allowed for a gene. Can influence how stringent mode selection is.
#' @param plot_heatmap Set to False to disable automatic display of activity matrix as a heatmap.
#'
#' @return List containing cell type annotation, matrix of average score by Leiden cluster and cell type used for attribution, cell annotation with Leiden clusters, activity matrix and UMAP computed on activity matrix.
#'
#' @export
#'
MAYA_predict_cell_types<-function(expr_mat,modules_list=NULL,min_cells_pct=0.05,organs=NULL,is_logcpm=F,nCores=1,thr=0,max_contrib=0.5,compute_umap=T,plot_heatmap=T,scale_before_pca=T){
    message("Running cell type identification")
    
    #### Check parameters ####
    stopifnot(is.logical(is_logcpm),
              is.numeric(min_cells_pct),is.numeric(nCores))
    if(!is.null(modules_list)){
        stopifnot(is.list(modules_list))
    }
    
    # load modules list in Panglao if necessary
    if(is.null(modules_list)){
        message("Loading PanglaoDB")
        
        path<-system.file("extdata", "PanglaoDB_markers_27_Mar_2020.tsv", package = "MAYA")
        panglao<-read.table(path,header=T,sep="\t",quote='')
        panglao<-panglao[which(panglao$species=="Mm Hs" | panglao$species=="Hs"),]
        
        # basic types
        panglao_basic<-panglao[which(panglao$organ %in% c("Connective tissue","Smooth muscle","Immune system","Vasculature","Blood","Epithelium","Skeletal muscle")),]
        basic_panglao<-list()
        for(cell_type in unique(panglao_basic$cell.type)){
            if(!(cell_type %in% c("Endothelial cells (blood brain barrier)","Endothelial cells (aorta)","Gamma delta T cells"))){
                basic_panglao[[cell_type]]<-panglao_basic[which(panglao_basic$cell.type==cell_type),"official.gene.symbol"]
            }
            
        }
        
        ### add lists for other organs if required
        if(!is.null(organs)){
            if(!is.null(intersect(unique(panglao$organ),organs))){
                panglao_specif<-panglao[which(panglao$organ %in% organs),]
                
                for(cell_type in unique(panglao_specif$cell.type)){
                    basic_panglao[[cell_type]]<-panglao_specif[which(panglao_specif$cell.type==cell_type),"official.gene.symbol"]
                }
                
            }
            if(organs=="all"){
                for(cell_type in unique(panglao$cell.type)){
                    basic_panglao[[cell_type]]<-panglao[which(panglao$cell.type==cell_type),"official.gene.symbol"]
                }
                
            }
            
        }
        
        modules_list<-basic_panglao
    }
    
    # run MAYA with cell type parameters
    suppressWarnings(PCA_obj<-run_activity_analysis(expr_mat = expr_mat,
                                                    modules_list =modules_list,
                                                    nb_comp_max = 1,
                                                    min_cells_pct = min_cells_pct,
                                                    min_module_size = 10,
                                                    max_contrib = max_contrib,
                                                    norm = !is_logcpm,
                                                    nCores = nCores,
                                                    scale_before_pca = scale_before_pca,
                                                    all_PCs_in_range=F))
    # predict cell annotations
    annot<-generate_cell_type_annotation(PCA_obj,thr=thr,compute_umap=compute_umap)
    
    if(plot_heatmap==T){
        print(plot_heatmap_activity_mat(activity_mat = annot$activity_matrix))
    }
    
    return(c(annot,list(PCA_obj=PCA_obj)))
}
