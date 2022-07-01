#' Build activity matrix and compute UMAP
#'
#' @param PCA_obj Vector, containing cluster annotation for individual cells
#'
#' @return List matrix of average score by Leiden cluster and pathway used for attribution, cell annotation with Leiden clusters, activity matrix and UMAP computed on activity matrix.
#'
#' @import RANN
#' @import igraph
#' @import leidenbase
#'
#' @export
#'
study_pathways<-function(PCA_obj,compute_umap=T){
    # compute activity mat
    activity_mat<-build_activity_mat(PCA_obj,scaled = F)

    # clustering and generate average score by matrix
    knn.matrix <- RANN::nn2(t(activity_mat), t(activity_mat), k = 20, searchtype = "standard")[[1]]

    jaccard.adj <- knn_jaccard(knn.matrix)

    graph<-igraph::graph.adjacency(
        jaccard.adj,
        mode = "undirected",
        weighted = TRUE
    )
    set.seed(2016)
    cluster_result <- leidenbase::leiden_find_partition(
        graph,
        verbose = FALSE,
        seed=2016,num_iter = 2, partition_type="ModularityVertexPartition",edge_weights=igraph::E(graph)$weight
    )
    cluster_result$membership<-paste0("C",cluster_result$membership)
    tmp<-average_by_cluster(activity_mat,cluster_result$membership)

    # compute umap
    if(compute_umap){
        umap<-run_umap(activity_mat)
    }
    else{umap<-NULL}

    return(list(cluster_matrix=tmp,clusters_annotation=cluster_result$membership,activity_matrix=activity_mat,umap=umap))

}

#' Run MAYA for unsupervised pathway acitivity analysis
#'
#' @param expr_mat Raw count matrix. If normalized, set is_logcpm to T
#' @param modules_list List with pathways associated with their genes. Can also set to "hallmark" or "kegg" to load corresponding MSigDB lists.
#' @param min_cells_pct Numeric between 0 and 1, minimum percentage of cells that should be above informativity threshold to consider activity score interesting.
#' @param is_logcpm Set to TRUE if data already normalized
#' @param nCores Number of cores to use. Set to 1 if working in a Windows environment, otherwise you can use the function detectCores() to find out how many cores are available.
#' @param thr Minimum average activity by cluster required to be assigned an annotation. Default 0, increase to be more stringent on assignation confidence.
#' @param max_contrib Numeric between 0 and 1, representing the maximum contribution to a PC allowed for a gene. Can influence how stringent mode selection is.
#'
#' @return List containing cell type annotation, matrix of average score by Leiden cluster and cell type used for attribution, cell annotation with MLeiden clusters, activity matrix and UMAp computed on activity matrix.
#'
#'
#' @export
#'
MAYA_pathway_analysis<-function(expr_mat,modules_list=NULL,min_cells_pct=0.05,is_logcpm=T,nCores=1,min_genes=10,max_contrib=0.5,compute_umap=T){
    message("Running pathway analysis")

    #### Check parameters ####
    stopifnot(is.logical(is_logcpm),
              is.numeric(min_cells_pct),is.numeric(nCores))

    # load modules list in Panglao if necessary
    if(is.null(modules_list)){
        message("Loading HALLMARK from MSigDB")
        path<-system.file("extdata", "h.all.v7.4.symbols.gmt", package = "MAYA")
        modules_list<-read_gmt(path)

    }
    if(is.character(modules_list)){
        if(modules_list=="hallmark"){
            message("Loading HALLMARK from MSigDB")
            path<-system.file("extdata", "h.all.v7.4.symbols.gmt", package = "MAYA")
            modules_list<-read_gmt(path)

        }
        else{
            if(modules_list=="kegg"){
                message("Loading KEGG from MSigDB")
                path<-system.file("extdata", "c2.cp.kegg.v7.4.symbols.gmt", package = "MAYA")
                modules_list<-read_gmt(path)
            }
        }

    }
    # at this stage, it can only be a list
    stopifnot(is.list(modules_list))

    # run MAYA with pathways parameters
    suppressWarnings(PCA_obj<-run_activity_analysis(expr_mat = expr_mat,
                                                                     modules_list =modules_list,
                                                                     nb_comp_max = 5,
                                                                     min_cells_pct = min_cells_pct,
                                                                     min_module_size = min_genes,
                                                                     max_contrib = max_contrib,
                                                                     norm = !is_logcpm,
                                                                     nCores = nCores))
    # analyze pathways
    annot<-study_pathways(PCA_obj,compute_umap=compute_umap)

    return(c(annot,list(PCA_obj=PCA_obj)))
}

