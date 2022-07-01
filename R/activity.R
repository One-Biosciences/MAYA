#' Compute pathway analysis using PCA spotting PCs driven by only one gene.
#'
#' @description For each module provided as input, runs PCA and tests informativity of successive PCs to find different activation scores.
#'
#'
#' @param expr_mat Numeric matrix, with genes as rows and cells as columns. This matrix should be normalized. If it is not, set parameter norm to TRUE to perform logCPM normalization.
#' @param modules_list List, with module names as list names and character vectors containing genes in module. Output from read_gmt file for example.
#' @param norm Boolean, set to TRUE if your data is not normalized. LogCPM normalization will be performed.
#' @param nb_comp_max Integer, maximum number of ways a module can be activated.
#' @param min_cells_pct Numeric between 0 and 1, minimum percentage of cells that should be above informativity threshold to consider activity score interesting.
#' @param nCores Number of cores to use. Set to 1 if working in a Windows environment, otherwise you can use the function detectCores() to find out how many cores are available.
#' @param min_module_size Minimum size of the gene set for which it makes sense to compute activity. Default is 10.
#'
#' @return List, for each module, stores activity scores (raw or scaled between 0 and 1) from informative PCs, gene contributions to each informative PC,
#' variance explained by each PC, threshold used for informativity and number of cells above threshold.
#'
#'
#' @importFrom stats prcomp
#' @import parallel
#' @import Matrix
#'
#' @export
#'
run_activity_analysis<-function(expr_mat,modules_list,norm=FALSE,nb_comp_max=5,min_cells_pct=0.05,nCores=1,min_module_size=10,max_contrib=0.5){

    message("Computing gene sets activity")

    #### Check parameters ####
    stopifnot(is.list(modules_list),
              is.logical(norm),is.numeric(nb_comp_max),
              is.numeric(min_cells_pct),is.numeric(nCores),
              is.numeric(min_module_size),is.numeric(max_contrib))
    stopifnot(is(expr_mat, 'sparseMatrix') | is(expr_mat, 'matrix'))

    expr_mat<-as(expr_mat,"dgCMatrix")

    #### Normalize data if specified ####
    if(norm){
        expr_mat<-logcpmNormalization(expr_mat)
    }
    # remove genes expressed in less than 10 cells
    expr_mat<-expr_mat[Matrix::rowSums(expr_mat!=0)>=10,]

    ### get minimum number of cells ###
    min_cells<-round(min_cells_pct*ncol(expr_mat))

    #### Run PCA for each gene set ####
    compute_gene_set_activity <- function(x) {
        common_genes<-intersect(modules_list[[x]],rownames(expr_mat))
        if(length(common_genes)>=min_module_size){
            # perform PCA, compute explained var and modify PCA sign to favor activation.
            pca<-stats::prcomp(scale(Matrix::t(expr_mat[common_genes,]),center=T,scale=T), scale = FALSE,retx = TRUE,rank. = 20)
            n=ncol(pca$x)
            ExpVar <- apply(pca$x[, 1:n], 2, var)/sum(apply(scale(Matrix::t(expr_mat[common_genes,]),center=T,scale=T), 2, var))
            # orientate PCs to favor activation (change orientation if genes with negative contributions have higher weights in absolute value than positive one)
            for(j in 1:min(min_module_size,nb_comp_max)){
                if(sum(pca$rotation[,j])<0){
                    pca$x[,j]<-(-pca$x[,j])
                    pca$rotation[,j]<-(-pca$rotation[,j])
                }
            }
            # scale projection for threshold detection
            projection<-scale_0_1(t(pca$x))
            # assess informativity of successive PCs and stop when conditions are unmet.
            comp<-c()
            list_thr<-c()
            nb_selected_cells<-c()
            for(i in 1:min(min_module_size,nb_comp_max)){
                if(length(comp)==0){
                    if(i>2){
                        break
                    }
                }
                if(length(comp)>0){
                    if(tail(comp,1)+1!=i){
                        break
                    }
                }
                if(max(abs(pca$rotation[,i]))<max_contrib){
                    thr<-.activity_assignmentThreshold(activity = projection[i,])
                    #stop when PC2 or more is uninformative
                    if(thr>1&i>=2){
                        break
                    }
                    #stop when PC2 or more is informative but doesn't distinguish enough active cells
                    if(i>=2&length(which(projection[i,]>=thr))<min_cells){
                        #if(length(which(projection[i,]>=thr))<min_cells){
                        break
                    }
                    #if(thr<1){
                    if(thr<1&length(which(projection[i,]>=thr))>=min_cells){
                        comp<-c(comp,i)
                        list_thr<-c(list_thr,thr)
                        nb_selected_cells<-c(nb_selected_cells,length(which(projection[i,]>=thr)))
                    }
                }
            }
            if(length(comp)!=0){
                nb_comp<-min(length(comp),nb_comp_max,min_module_size)
                # keep only PCs that pass threshold and keep scaled score
                if(length(comp)==1){
                    activity_scores<-as.matrix(t(projection[comp,])) # cells as columns
                    rownames(activity_scores)<-paste0("mode",comp)
                }else{
                    activity_scores<-projection[comp,]
                    rownames(activity_scores)<-paste0("mode",comp)
                }
                activity_scores_raw<-t(pca$x[,comp])
                rownames(activity_scores_raw)<-paste0("mode",comp)
                gene_contrib<-t(pca$rotation[,comp]) # genes as columns
                expl_var<-ExpVar[comp]
                list(activity_scores=activity_scores,activity_scores_raw=activity_scores_raw,gene_contrib=gene_contrib,expl_var=expl_var,list_thr=list_thr,nb_selected_cells=nb_selected_cells)
            }
        }
    }


    output_object<-parallel::mclapply(names(modules_list), FUN = compute_gene_set_activity,mc.cores = nCores)


    #### Generate final output ####
    names(output_object)<-names(modules_list)
    #remove modules with no informative scores from the output object.
    output_object<-output_object[lapply(output_object,length)!=0]

    message("Found at least one informative activation mode for ",length(output_object)," gene sets")

    return(output_object)
}



#' Extract activity matrix from PCA object
#'
#' @param PCA_object List, output from run_activity_analysis().
#' @param scaled Boolean. By default, using scores scaled between 0 and 1. Set to FALSE if you want to keep raw projections of cells on PCs.
#'
#' @return Numeric matrix, with the different activation ways as rows and cells as columns and containing the activity score.
#'
#' @importFrom dplyr bind_rows
#'
#' @export
#'
build_activity_mat<-function(PCA_object,scaled=T){

    #### Check parameters ####
    stopifnot(is.list(PCA_object),is.logical(scaled))

    #### keep only activity scores ####
    if(scaled){
        activity_list<-lapply(names(PCA_object), function(x) PCA_object[[x]]$activity_scores)
    }else{
        activity_list<-lapply(names(PCA_object), function(x) PCA_object[[x]]$activity_scores_raw)
    }
    names(activity_list)<-names(PCA_object)
    # rename PCs with module name
    activity_list<-lapply(names(activity_list), function(x) {
        rownames(activity_list[[x]])<-paste(x,rownames(activity_list[[x]]),sep="_")
        as.data.frame(activity_list[[x]])
    })

    #### Generate final output ####
    aggregate_mat<-as.matrix(dplyr::bind_rows(activity_list))
    return(aggregate_mat)
}





