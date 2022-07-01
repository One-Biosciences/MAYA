#' Compute specificity table given activity score matrix
#'
#' @description This function computes the specificity of activity of each gene set activation mode among groups of
#' cells that are determined by a given annotation.The specificity scores is comprised between 0 and 1, and the sum of
#' activities across clusters is 1. If a module is active only in one cluster, the specificity for this cluster will
#' be 1 and 0 for the others.
#'
#' @param activity_mat Numeric matrix, output from compute_activity function.
#' @param meta Dataframe containing metadata, with cell labels as rownames.
#' @param annot_name Character, name of the column in meta that should be used to compute specificity for.
#'
#' @return Dataframe containing specificity score for each cluster for each row of the activity matrix.
#'
#' @export
#'
specificity_table<-function(activity_mat,meta,annot_name){
    
    #### Check parameters ####
    stopifnot(is.data.frame(meta),is.character(annot_name))
    
    #### Get annotation ####
    cellTypes <- sort(as.factor(unique(meta[,annot_name])))
    
    #### Compute table of average score by annotation ####
    mean_table<-matrix(0, ncol=length(cellTypes),nrow=nrow(activity_mat))
    colnames(mean_table)<-cellTypes
    rownames(mean_table)<-rownames(activity_mat)
    for(celltype in cellTypes){
        cells<-rownames(meta[which(meta[,annot_name]==celltype),])
        mean_table[,celltype]<-rowMeans(activity_mat[,cells])
    }
    
    #### Compute specificity ####
    spe_table<-t(apply(mean_table, 1, function(x) x^2/sum(x^2)))
    colnames(spe_table)<-paste0("S_",colnames(spe_table))
    return(as.data.frame(spe_table))
}


#' Compute homogeneity table given activity score matrix
#'
#' @description This function computes the homogeneity of activity of each gene set activation mode among groups of
#' cells that are determined by a given annotation.The homogeneity scores is comprised between 0 and 1, 1 corresponding to a
#' highly homogeneous activity among cells from this given cluster. The standard deviation of the scores for this module and
#' cluster is compared with the standard deviation of a distribution with the same number of observations but with half values
#' equal to 0 and the other half equals to 1 (worse possible distribution in terms of homogeneity). Careful, if the module is inactive
#' in all cells, it will be highly homogeneous, homogeneity is not linked to the level of activation.
#'
#' @param activity_mat Numeric matrix, output from compute_activity function.
#' @param meta Dataframe containing metadata, with cell labels as rownames.
#' @param annot_name Character, name of the column in meta that should be used to compute homogeneity for.
#'
#' @return Dataframe containing homogeneity score for each cluster for each row of the activity matrix.
#'
#' @export
#'
homogeneity_table<-function(activity_mat,meta,annot_name){
    
    #### Check parameters ####
    stopifnot(is.data.frame(meta),is.character(annot_name))
    
    #### Get annotation ####
    cellTypes <- sort(as.factor(unique(meta[,annot_name])))
    
    #### Compute homogeneity for each cell type, and for each module #### 
    sd <- lapply(cellTypes, function(thisType)
        sapply(rownames(activity_mat), function(thisModule)
        {
            idx_to_keep<-which(meta[,annot_name]==thisType)
            pModule <- activity_mat[thisModule,idx_to_keep] # normalized vector corresponding to module
            #compute maximum sd
            n<-length(pModule)
            n_zeros<-round(n/2)
            n_ones<-n-n_zeros
            max_sd<-sd(c(rep(0,n_zeros),rep(1,n_ones)))
            #compute homogeneity score
            1-sd(pModule)/max_sd
        })
    )
    sd <- do.call(cbind, sd)
    colnames(sd) <- paste0("H_",cellTypes)
    return(as.data.frame(sd))
}
