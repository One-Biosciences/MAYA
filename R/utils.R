
#' Transform raw count matrix into logCPM matrix
#'
#' @description Converts continuous BMI variable to BMI groups.
#'
#' @param X Numeric matrix.
#'
#' @return Numeric matrix.
#'
#' @examples
#' raw_mat <- matrix(1:20, ncol = 5, nrow = 4)
#' norm_mat <- logcpmNormalization(raw_mat)
#' @import Matrix
#'
#' @export
logcpmNormalization <- function(X) {
    X <- Matrix::t(Matrix::t(X) / Matrix::colSums(X)) * 1e6
    X <- log1p(X)
    # Handle unfiltered matrix, 0/0 -> NaN is not compatible with sparse matrix operations
    if (!is.null(attr(X, "x"))) {
        X@x[is.nan(X@x)] <- 0
    } else {
        X[is.nan(X)] <- 0
    }
    # X@x <- log2(X@x + 1) if sparse
    # return(as.matrix(X))
    return(X) # to allow sparse format that works with prcomp for PCA computation
}

#' Scale vector between 0 and 1
#'
#' Changes the range of values of a vector keeping the shape of the distribution. Minimum value is set to 0 and maximum to 1.
#'
#' @param X Numeric vector.
#'
#' @return Numeric vector.
#'
#'
#' @examples
#' vect <- rnorm(n = 50, mean = 25, sd = 3)
#' cal_z_score(vect)
#' @noRd
cal_z_score <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}

#' Scale matrix rows between 0 and 1
#'
#' The minimum value of each row is set to 0 and the maximum to 1.
#'
#' @param X Numeric matrix.
#'
#' @return Numeric matrix.
#'
#'
#' @examples
#' raw_mat <- matrix(1:20, ncol = 5, nrow = 4)
#' norm_mat <- scale_0_1(raw_mat)
#' @export
#'
scale_0_1 <- function(X) {
    X <- t(apply(X, 1, cal_z_score))
    return(X)
}



#' Read GMT file
#'
#' Imports specified GMT file as a list
#'
#' @param gmt_file_name Character.
#'
#' @return List of modules, each associated with a character vector containing genes constituting the module.
#'
#' @examples
#' read_gmt(system.file("extdata", "h.all.v7.4.symbols.gmt", package = "MAYA"))
#' @export
read_gmt <- function(gmt_file_name) {
    #### Check parameters ####
    stopifnot(file.exists(gmt_file_name))

    #### Read GMT file ####
    gmt <- strsplit(readLines(gmt_file_name), "\t")
    names(gmt) <- sapply(gmt, `[`, 1)
    gmt <- lapply(gmt, function(x) {
        genes <- x[-c(1, 2)]
    })
    return(gmt)
}



#' Compute Jaccard distance from kNN matrix
#'
#' @description Compute Jaccard distance from kNN matrix
#'
#' @param m Numeric/sparse matrix.
#'
#' @return Numeric/sparse matrix.
#'
#' @import Matrix
jaccard <- function(m) {
    ## common values:
    A <- Matrix::tcrossprod(m)
    A <- as(A, "dgTMatrix")
    ## counts for each row
    b <- Matrix::rowSums(m)
    ## Jacard formula: #common / (#i + #j - #common)
    x <- A@x / (b[A@i + 1] + b[A@j + 1] - A@x)
    A@x <- x
    return(A)
}

#' Compute Jaccard distance from kNN results
#'
#' @description Compute Jaccard distance from kNN results
#'
#' @param knn Numeric matrix (cells as rows, index of neighbors as columns (20 columns if k=20))
#'
#' @return Numeric/sparse matrix.
#'
#' @import Matrix
knn_jaccard <- function(knn) {
    knn.df <- data.frame(i = rep(1:nrow(knn), ncol(knn)), j = as.vector(knn))
    knn.mat <- Matrix::sparseMatrix(i = knn.df[[1]], j = knn.df[[2]], x = 1)
    jaccard(knn.mat)
}

#' Average matrix by column
#'
#' @description Compute for each row the average value by provided annotation
#'
#' @param mat Numeric/sparse matrix (cells as columns)
#' @param annotation Vector containg cell annotation
#'
#' @return Numeric/sparse matrix.
#'
#' @import dplyr
average_by_cluster <- function(mat, annotation) {

    # get all possible clusters
    clusters <- unique(annotation)

    # compute row means for each subsetted matrix by cell type and store it as a list
    means <- lapply(clusters, function(x) {
        if (length(which(annotation == x)) > 1) {
            tmp <- mat[, which(annotation == x)]

            if (!is.null(dim(tmp)) && dim(tmp)[1] > 1) {
                rowMeans(tmp)
            } else {
                as.matrix(mean(tmp))
            }
        } else {
            mat[, which(annotation == x)]
        }
    })
    # store result as a matrix
    names(means) <- clusters
    out <- as.matrix(dplyr::bind_cols(means))
    rownames(out) <- rownames(mat)

    return(out)
}


#' Run UMAP with fixed seed for reproducibility
#'
#' @param mat Numeric matrix, output from compute_activity function.
#'
#' @return UMAP object
#'
#' @import umap
#'
#' @export
run_umap <- function(mat) {

    #### Check parameters ####
    stopifnot(is.numeric(mat))

    #### Compute UMAP with fixed seed ####
    set.seed(123)
    umap(t(mat))
}


#' Generate custom color palette
#'
#' @param length Number of colors required
#'
#' @return Vector of colors
#'
#' @import wesanderson
#'
#' @export
generate_palette <- function(length) {
    set.seed(6)
    palette <- unique(c(
        wesanderson::wes_palette("Rushmore1")[3],
        wesanderson::wes_palette("Rushmore1")[5],
        wesanderson::wes_palette("Zissou1")[1],
        wesanderson::wes_palette("Zissou1")[3],
        wesanderson::wes_palette("Darjeeling1")[4],
        wesanderson::wes_palette("GrandBudapest2")[1],
        sample(c(wesanderson::wes_palette("Darjeeling1"), wesanderson::wes_palette("Darjeeling2")[1:4], wesanderson::wes_palette("Moonrise3")[1:3], wesanderson::wes_palette("GrandBudapest1"), wesanderson::wes_palette("GrandBudapest2")))
    ))
    if (length <= length(palette)) {
        palette <- palette[1:length]
    } else {
        palette <- rep(palette, times = round(length / length(palette)) + 1)
        palette <- palette[1:length]
    }
    return(palette)
}
