require(vegan)
require(scales)
require(ape)
require(phyloseq)


#' Compute incidence matrix of a tree
#'
#' @title incidenceMatrix
#' @param phy Required. A \code{phylo} class object
#' @return An incidence matrix M of size nedges(tree) x ntaxa(tree) where
#'         M[i,j] is set to 1 if taxa derives from edge i and 0 otherwise.
#' @note incidenceMatrix assumes that the tree is rooted. If not, it roots it
#'       it at an arbitrary taxa (taxa 1). 
incidenceMatrix <- function(phy) {
    if (!is.rooted(phy)) {
        warning("Tree is not rooted, incidence matrix may be meaningless")        
    }
    ## Construct incidence matrix of the tree (taxa x edge matrix)
    ## All taxa descending from an edge are set to 1, all others to -1
    ntaxa <- length(phy$tip.label)
    nedges <- nrow(phy$edge)
    incidence <- matrix (0,
                         nrow = ntaxa,
                         ncol = nedges,
                         dimnames = list(phy$tip.label, phy$edge[, 2]))
    ## Incidence of internal edges
    phy.part <- prop.part(phy) ## clade composition indexed by (shifted) node number
    for (i in 2:length(phy.part)) { ## first clade corresponds to root node
        edge <- which(phy$edge[, 2] == i + ntaxa) ## node numbers are shifted by ntaxa 
        incidence[phy.part[[i]] , edge] <- 1
    }
    ## Incidence of pendant edges
    ## pendant.edges[i] is the edge leading to tip i. 
    pendant.edges <- match(seq_len(ntaxa), phy$edge[ , 2])
    for (i in seq_len(ntaxa)) {
        incidence[i, pendant.edges[i]] <- 1
    }
    attr(incidence, "pendant.edges") <- pendant.edges
    return(incidence)
}


#' Compute edge PCA components of a metagenomic sample
#' using the edgePCA idea developed in Matsen and Evans (2012)
#' with regularization of the edges for the PCA part as proposed in
#' Shen, H. and Huang, J. Z. (2008).
#'
#' @title edgePCA
#' @param physeq Required. A \code{\link{phyloseq-class}} class object
#' @param number.components Optional. Numeric. Number of components used in the PCA, defaults to 2.
#' @param method. Optional. PCA method. Any unambiguous abbreviation of 'normal',
#'                'sparse' and 'regularized'. Defaults to 'sparse'.
#' @param number.edges Optional. Integer. Number of non-zero loadings for each
#'                     component in the sparse PCA setting.
#' @param ...  Optional. Additional parameters passed on to sparcePca,
#'             regularizedPca and pca, such as number of iterations and
#'             convergence criterion.
#' @return An \code{edgepca} class object with components:
#' - \item{loadings}: The loadings of the principal component axis
#' - \item{scores}: The score of the samples
#' - \item{values}: A date frame with components eigenvalues, relative eigenvalues and
#'                  cumulative eigenvalues.
#' @note scores and extract_eigenvalues methods are associated with 'edgepca'. 
#' @references Shen, H. and Huang, J. Z. (2008). Sparse principal component
#' analysis via regularized low rank matrix approximation. _Journal
#' of Multivariate Analysis_ *99*, 1015-1034.
#' @references Matsen, F. A. and Evans, S. (2012). Edge Principal Components
#' and Squash Clustering: Using the Special Structure of Phylogenetic
#' Placement Data for Sample Comparison. _Plos One_ *8*(3):e56859.
#' @seealso \code{\link{pca}}, \code{\link{sparsePca}}, \code{\link{regularizedPca}},
edgePCA <- function(physeq,
                    number.components = 2,
                    method = c("normal", "sparse", "regularized"),
                    number.edges = 50,
                    ...) {

    ## Extract otu_table matrix
    x <- as(otu_table(physeq), "matrix")
    if (taxa_are_rows(physeq)) { x <- t(x) }
    phy <- phy_tree(physeq)
    
    ## Scale counts to frequencies
    x <- x/rowSums(x)

    ## Construct incidence matrix of the tree and transform it to contrasts
    incidence <- incidenceMatrix(phy)
    incidence <- 2 * (incidence - 0.5)
            
    ## Order community table according to edge order and create
    ## mass difference matrix
    x <- x[ , rownames(incidence), drop = FALSE]
    mdm <- x %*% incidence

    ## Correct number.edges if needed
    if (length(number.edges) == 1) {
        number.edges <- rep(number.edges, number.components)
    }
    if (length(number.edges) != number.components) {
        warning("number.edges is not long enough, set to repeats of number.edges[1]")
        number.edges <- rep(number.edges[1], number.components)
    }

    
    ## Compute number.component first axis from the mass difference matrix
    method <- match.arg(method)
    results <- switch(method,
                      sparse = sparsePca(mdm, number.components, number.edges, ...),
                      normal = pca(mdm, number.components, ...),
                      regularized = regularizedPca(mdm, number.components, ...))
    class(results) <- c("edgepca", class(results))
    return(results)
}




#' Perform a principal component analysis on a given data matrix using
#' Singular Value Decomposition.
#'
#' @title pca
#' @param x Required. Numeric matrix.
#' @param number.components Optional. Integer. Number of components returned by the PCA.
#'                          Defaults to NULL which returns all components
#' @param center Optional. Logical. Should x be shifted to be (column-wise) 0-centered?
#'               Defaults to TRUE. Alternatively, a vector of length equal to the
#'               number of columns in x. The value is passed on to scale.
#' @param scale. Optional. Logical. Should the variable be normalized to have unit
#'               variance. Defaults to TRUE. Alternatively, a vector of length
#'               equal to the number of columns in x. The value is passed on to scale.
#' @return A list with components
#' - \item{number.components}: The number of components used
#' - \item{loadings}: The matrix of variable loadings
#' - \item{scores}: The matrix of sample scores
#' - \item{values}: A list with components
#'     - \item{Eigenvalues}: the eigenvalues of the covariance/correlation
#'     - \item{Relative_eig}: the normalized eigenvalues (sum up to 1)
#'     - \item{Cumul_eig}: the cumulative eigenvalues (useful for scree plot)
#' - \item{center, scale}: centering and scaling vectors, used, if any.
#'                         Otherwise, FALSE.
#' @seealso \code{\link{edgepca}}, \code{\link{sparsePca}}, \code{\link{regularizedPca}},
pca <- function(x, number.components = NULL, center = TRUE, scale. = FALSE) {
    x <- as.matrix(x)
    ## Construct variable and sample names
    variable.names <- colnames(x)
    if (is.null(variable.names)) {
        variable.names <- paste("var", 1:ncol(x), sep = ".")
    }
    sample.names <- rownames(x)
    if (is.null(sample.names)) {
        sample.names <- 1:nrow(x)
    }
    ## Scale matrix
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if (any(sc == 0)) 
        stop("cannot rescale a constant/zero column to unit variance")
    ## Check for missing values
    if (any(is.na(x))) {
        stop("cannot handle missing values")
    }
    ## Optionally, correct number of components
    number.components <- min(number.components, ncol(x), nrow(x))
    ## apply svd
    x.svd <- svd(x, nu = 0, nv = number.components)
    ## Extract eigenvectors (loadings), eigenvalues and scores
    loadings <- x.svd$v
    scores <- x %*% loadings
    eigenvalues <- (x.svd$d)^2 / max(1, nrow(x) - 1)
    ## Add sensible names to loadings and scores
    dimnames(loadings) <- list(variable.names,
                               paste("Axis", 1:ncol(loadings), sep = "."))
    dimnames(scores) <- list(sample.names,
                             paste("Axis", 1:ncol(loadings), sep = "."))
    ## Compute relative and cumulative eigenvalues
    values <- list(Eigenvalues = eigenvalues,
                   Relative_eig = eigenvalues / sum(eigenvalues),
                   Cumul_eig = cumsum(eigenvalues) / sum(eigenvalues))
    ## Construct results and return it
    result <- list(number.components = number.components,
                   scores = scores,
                   loadings = loadings,
                   values = values,
                   center = ifelse(center, cen, FALSE),
                   scale = ifelse(scale., sc, FALSE))
    class(result) <- "pca"
    return(invisible(result))
}


#' Perform a regularized principal component analysis on a given data matrix using
#' Singular Value Decomposition and regularization. 
#'
#' @title regularizedPca
#' @param x Required. Numeric matrix.
#' @param number.components Optional. Integer. Number of components used for denoising.
#'                          Defaults to 10 and is automatically reduced to max(dim(x)) - 1
#'                          if 10 is too high. 
#' @param center Optional. Logical. Should x be shifted to be (column-wise) 0-centered?
#'               Defaults to TRUE. Alternatively, a vector of length equal to the
#'               number of columns in x. The value is passed on to scale.
#' @param scale. Optional. Logical. Should the variable be normalized to have unit
#'               variance. Defaults to TRUE. Alternatively, a vector of length
#'               equal to the number of columns in x. The value is passed on to scale.
#' @return A list with components
#' - \item{number.components}: The number of components used
#' - \item{loadings}: The matrix of variable loadings
#' - \item{scores}: The matrix of sample scores
#' - \item{values}: A list with components
#'     - \item{Eigenvalues}: the eigenvalues of the covariance/correlation
#'     - \item{Relative_eig}: the normalized eigenvalues (sum up to 1)
#'     - \item{Cumul_eig}: the cumulative eigenvalues (useful for scree plot)
#' - \item{center, scale}: centering and scaling vectors, used, if any.
#'                         Otherwise, FALSE.
#' @seealso \code{\link{pca}}, \code{\link{edgePCA}}, \code{\link{regularizedPca}}
regularizedPca <- function(x, number.components = 10, center = TRUE, scale. = FALSE) {
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ## Construct variable and sample names
    variable.names <- colnames(x)
    if (is.null(variable.names)) {
        variable.names <- paste("var", 1:p, sep = ".")
    }
    sample.names <- rownames(x)
    if (is.null(sample.names)) {
        sample.names <- 1:n
    }
    ## Scale matrix
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if (any(sc == 0)) 
        stop("cannot rescale a constant/zero column to unit variance")
    ## Check for missing values
    if (any(is.na(x))) {
        stop("cannot handle missing values")
    }
    ## Choose number of components
    if (is.null(number.components)) {
        stop("Please specify the number of components used for denoising")
    } else {
        number.components <- min(number.components, ncol(x) - 1, nrow(x) - 1)
    }
    ## apply svd
    x.svd <- svd(x, nu = number.components, nv = number.components)
    ## Shrink eigenvalues
    S <- number.components
    sigma.square <- sum(x.svd$d[-c(1:S)]^2)/(n*p - p - n*S - p*S + S^2)
    shrinked.eigenvalues = (x.svd$d[1:S]^2 - ifelse(p < n, n, p * (1 - 1/n)) * sigma.square) / x.svd$d[1:S]
    ## Extract eigenvectors (loadings), eigenvalues and scores
    loadings <- x.svd$v
    scores <- x.svd$u %*% diag(shrinked.eigenvalues)
    eigenvalues <- shrinked.eigenvalues
    ## Add sensible names to loadings and scores
    dimnames(loadings) <- list(variable.names,
                               paste("Axis", 1:ncol(loadings), sep = "."))
    dimnames(scores) <- list(sample.names,
                               paste("Axis", 1:ncol(loadings), sep = "."))
    ## Compute relative and cumulative eigenvalues
    values <- list(Eigenvalues = eigenvalues,
                   Relative_eig = eigenvalues / sum(eigenvalues),
                   Cumul_eig = cumsum(eigenvalues) / sum(eigenvalues))
    ## Construct results and return it
    result <- list(number.components = number.components,
                   scores = scores,
                   loadings = loadings,
                   values = values,
                   center = ifelse(center, cen, FALSE),
                   scale = ifelse(scale., sc, FALSE))
    class(result) <- c("rpca", "pca")
    return(invisible(result))
}


#' Perform a sparse principal component analysis on a given data matrix using
#' NIPALS algorithm. 
#'
#' @title sparsePca
#' @param x Required. Numeric matrix.
#' @param number.components Optional. Integer. Number of components used for denoising.
#'                          Defaults to 10 and is automatically reduced to max(dim(x)) - 1
#'                          if 10 is too high. 
#' @param center Optional. Logical. Should x be shifted to be (column-wise) 0-centered?
#'               Defaults to TRUE. Alternatively, a vector of length equal to the
#'               number of columns in x. The value is passed on to scale.
#' @param scale. Optional. Logical. Should the variable be normalized to have unit
#'               variance. Defaults to TRUE. Alternatively, a vector of length
#'               equal to the number of columns in x. The value is passed on to scale.
#' @param number.variables Optional. Integer vector of length number.components specifying
#'                         the number of non-zero loadings for each component. By default, all
#'                         loadings can be different from 0 (no sparsity).
#' @param iter.max Optional. Integer. The maximum number of iterations used in NIPALS for each
#'                 each component. Defaults to 500. 
#' @param tol    Optional. A (small) positive numeric used to check convergence of loadings
#'               for each component. Defaults to 1e-09. 
#' @return A list with components
#' - \item{number.components}: The number of components used
#' - \item{loadings}: The matrix of variable loadings
#' - \item{scores}: The matrix of sample scores
#' - \item{values}: A list with components
#'     - \item{Eigenvalues}: the eigenvalues of the covariance/correlation
#'     - \item{Relative_eig}: the normalized eigenvalues (sum up to 1)
#'     - \item{Cumul_eig}: the cumulative eigenvalues (useful for scree plot)
#' - \item{center, scale}: centering and scaling vectors, used, if any.
#'                         Otherwise, FALSE.
#' @seealso \code{\link{pca}}, \code{\link{edgePCA}}, \code{\link{regularizedPca}},
sparsePca <- function(x, number.components = 10, 
                      number.variables = rep(ncol(x), number.components),
                      center = TRUE, scale. = FALSE,
                      iter.max = 500, tol = 1e-09) {
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ## Construct variable and sample names
    variable.names <- colnames(x)
    if (is.null(variable.names)) {
        variable.names <- paste("var", 1:p, sep = ".")
    }
    sample.names <- rownames(x)
    if (is.null(sample.names)) {
        sample.names <- 1:n
    }
    ## Scale matrix
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if (any(sc == 0)) 
        stop("cannot rescale a constant/zero column to unit variance")
    ## Check for missing values
    if (any(is.na(x))) {
        stop("cannot handle missing values")
    }
    ## Check consistency of number.variables
    if (length(number.variables) != number.components) {
        stop("number.variables must be a vector of length ", number.components, ".")
    }
    if (any(number.variables > p)) {
        wrong.components <- which(number.variables > p)
        warning("You selected more than p variables for components",
                paste(wrong.components, collapse = ", "),
                ". Automatically changed to p.")
        number.variables[wrong.components] <- p
    }
    ## Check number of components
    number.components <- min(number.components, n, p)
    ## Set up u (scores), v (loadings) and d (eigenvalues) matrices for
    ## the NIPALS algorithm
    vect.d <- rep(NA, number.components)
    iterations <- rep(0, number.components)
    convergence <- rep(FALSE, number.components)
    mat.u <- matrix(NA, nrow = n, ncol = number.components)
    mat.v <- matrix(NA, nrow = p, ncol = number.components)
    x.new <- x
    for (h in 1:number.components) {
        ## Initialization 
        x.svd <- svd(x.new, nu = 1, nv = 1)
        u.new <- x.svd$u[ , 1]
        v.new <- x.svd$d[1] * x.svd$v[ , 1]
        v.conv <- u.conv <- FALSE
        if (h >= 2) {
            x.h = x %*% mat.v[, 1:(h - 1)]
        }
        variables.h <- p - number.variables[h] ## number of variables to discard (at least)
        iter <- 0
        ## while loop
        while ( (!v.conv || !u.conv) & iter < iter.max) {
            iter <- iter + 1
            v.temp <- t(x.new) %*% u.new
            u.old <- u.new
            v.old <- v.new
            ## Update and sparsify v by keeping only p - nx variables
            if (variables.h != 0) {
                threshold <- sort(abs(v.temp))[variables.h]
                v.new <- ifelse(abs(v.temp) > threshold,
                                v.temp - sign(v.temp) * threshold,
                                0)
            }
            ## Update and normalize u
            if (h == 1) {
                u.new <- as.vector(x.new %*% v.new)
                u.new <- u.new /sqrt(sum(u.new^2))
            } else {
                u.new <- lsfit(y = x %*% v.old,
                               x = x.h,
                               intercept = FALSE,
                               tolerance = tol)$res
                u.new = u.new/sqrt(sum(u.new^2))
            }
            ## Check convergence
            if (sum( (u.new - u.old)^2 ) < tol) {
                u.conv <- TRUE
            }
            if (sum( (v.new - v.old)^2 ) < tol) {
                v.conv <- TRUE
            }
        }
        ## Keep log of convergence
        iterations[h] <- iter
        convergence[h] <- (v.conv && u.conv)
        ## Norm v.new, update mat.u, mat.v
        mat.v[ , h] <- v.new / sqrt(sum(v.new^2))
        mat.u[ , h] <- u.new
        ## Deflate x.new with first axis
        x.new <- x.new - x.svd$d[1] * tcrossprod(x.svd$u, x.svd$v)
        ## Compute percentage of variance explained
        ## for x reconstructed with first h components
        x.reconstructed.h <- x %*% mat.v[ , 1:h, drop = F] %*%
            solve(t(mat.v[ , 1:h, drop = F]) %*% mat.v[ , 1:h, drop = F]) %*% t(mat.v[ , 1:h, drop = F])
        vect.d[h] <- sum(x.reconstructed.h^2)
    }
    ## Extract eigenvectors (loadings), eigenvalues and scores
    loadings <- mat.v
    scores <- mat.u
    eigenvalues <- -diff(c(0, vect.d)) 
    ## Add sensible names to loadings and scores
    dimnames(loadings) <- list(variable.names,
                               paste("Axis", 1:ncol(loadings), sep = "."))
    dimnames(scores) <- list(sample.names,
                               paste("Axis", 1:ncol(loadings), sep = "."))
    ## Compute relative and cumulative eigenvalues
    values <- list(Eigenvalues = eigenvalues,
                   Relative_eig = eigenvalues / sum(eigenvalues),
                   Cumul_eig = cumsum(eigenvalues) / sum(eigenvalues))
    ## Construct results and return it
    result <- list(number.components = number.components,
                   scores = scores,
                   loadings = loadings,
                   values = values,
                   center = ifelse(center, cen, FALSE),
                   scale = ifelse(scale., sc, FALSE))
    class(result) <- c("sparsePca", "pca")
    return(invisible(result))
}

#' Perform a sparse principal component analysis on a given data matrix using
#' iterative algorithm from Shen and Huang. 
#'
#' @title sparsePca
#' @param x Required. Numeric matrix.
#' @param number.components Optional. Integer. Number of components used for denoising.
#'                          Defaults to 10 and is automatically reduced to max(dim(x)) - 1
#'                          if 10 is too high. 
#' @param center Optional. Logical. Should x be shifted to be (column-wise) 0-centered?
#'               Defaults to TRUE. Alternatively, a vector of length equal to the
#'               number of columns in x. The value is passed on to scale.
#' @param scale. Optional. Logical. Should the variable be normalized to have unit
#'               variance. Defaults to TRUE. Alternatively, a vector of length
#'               equal to the number of columns in x. The value is passed on to scale.
#' @param lambda Required. Regularization parameter.
#' @param iter.max Optional. Integer. The maximum number of iterations used in iterative step
#'                 for each
#'                 each component. Defaults to 500. 
#' @param tol    Optional. A (small) positive numeric used to check convergence of loadings
#'               for each component. Defaults to 1e-09.
#' @references Shen, H. and Huang, J. Z. (2008). Sparse principal component
#' analysis via regularized low rank matrix approximation. _Journal
#' of Multivariate Analysis_ *99*, 1015-1034.
#' @return A list with components
#' - \item{number.components}: The number of components used
#' - \item{loadings}: The matrix of variable loadings
#' - \item{scores}: The matrix of sample scores
#' - \item{values}: A list with components
#'     - \item{Eigenvalues}: the eigenvalues of the covariance/correlation
#'     - \item{Relative_eig}: the normalized eigenvalues (sum up to 1)
#'     - \item{Cumul_eig}: the cumulative eigenvalues (useful for scree plot)
#' - \item{center, scale}: centering and scaling vectors, used, if any.
#'                         Otherwise, FALSE.
#' @seealso \code{\link{pca}}, \code{\link{edgePCA}}, \code{\link{regularizedPca}},
sparsePca.lambda <- function(x, number.components = 10, 
                             lambda,
                             center = TRUE, scale. = FALSE,
                             iter.max = 500, tol = 1e-09) {
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ## Construct variable and sample names
    variable.names <- colnames(x)
    if (is.null(variable.names)) {
        variable.names <- paste("var", 1:p, sep = ".")
    }
    sample.names <- rownames(x)
    if (is.null(sample.names)) {
        sample.names <- 1:n
    }
    ## Scale matrix
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if (any(sc == 0)) 
        stop("cannot rescale a constant/zero column to unit variance")
    ## Check for missing values
    if (any(is.na(x))) {
        stop("cannot handle missing values")
    }
    ## Check consistency of number.variables
    if (length(number.variables) != number.components) {
        stop("number.variables must be a vector of length ", number.components, ".")
    }
    if (any(number.variables > p)) {
        wrong.components <- which(number.variables > p)
        warning("You selected more than p variables for components",
                paste(wrong.components, collapse = ", "),
                ". Automatically changed to p.")
        number.variables[wrong.components] <- p
    }
    ## Check number of components
    number.components <- min(number.components, n, p)
    ## Set up u (scores), v (loadings) and d (eigenvalues) matrices for
    ## the iterative algorithm
    vect.d <- rep(NA, number.components)
    iterations <- rep(0, number.components)
    convergence <- rep(FALSE, number.components)
    mat.u <- matrix(NA, nrow = n, ncol = number.components)
    mat.v <- matrix(NA, nrow = p, ncol = number.components)
    x.new <- x
    for (h in 1:number.components) {
        ## Initialization 
        x.svd <- svd(x.new, nu = 1, nv = 1)
        ## u and v solutions of unpenalized approximation problem
        ## u %*% t(v) is the best rank one approximation of x
        u.new <- x.svd$u[ , 1]
        v.new <- x.svd$d[1] * x.svd$v[ , 1]
        v.conv <- u.conv <- FALSE
        iter <- 0
        ## while loop
        while ( (!v.conv || !u.conv) & iter < iter.max) {
            iter <- iter + 1
            u.old <- u.new
            v.old <- v.new
            ## update v using soft thresholding
            z <- t(x.new) %*% u.new
            v.new <- ifelse(abs(z) > lambda, sign(z) * (abs(z) - lambda), 0)
            ## Update and normalize u
            u.new <- as.vector(x.new %*% v.new)
            u.new <- u.new / sqrt(sum(u.new^2))            
            ## Check convergence
            if (sum( (u.new - u.old)^2 ) < tol) {
                u.conv <- TRUE
            }
            if (sum( (v.new - v.old)^2 ) < tol) {
                v.conv <- TRUE
            }
        }
        ## Keep log of convergence
        iterations[h] <- iter
        convergence[h] <- (v.conv && u.conv)
        ## Norm v (?) and update mat.u, mat.v
        ## v.new <- v.new / sqrt(sum(v.new^2))
        mat.v[ , h] <- v.new 
        mat.u[ , h] <- u.new
        ## Deflate x.new with rank one approximation 
        x.new <- x.new - tcrossprod(u.new, v.new)
        ## Compute percentage of variance explained
        ## for x reconstructed with first h components
        x.reconstructed.h <- mat.u[ , 1:h, drop = F] %*% mat.v[ , 1:h, drop = F]
        vect.d[h] <- sum(x.reconstructed.h^2) / sum(x^2)
    }
    ## Extract eigenvectors (loadings), eigenvalues and scores
    loadings <- mat.v
    scores <- mat.u
    eigenvalues <- diff(c(0, vect.d, 1)) 
    ## Add sensible names to loadings and scores
    dimnames(loadings) <- list(variable.names,
                               paste("Axis", 1:ncol(loadings), sep = "."))
    dimnames(scores) <- list(sample.names,
                               paste("Axis", 1:ncol(loadings), sep = "."))
    ## Compute relative and cumulative eigenvalues
    values <- list(Eigenvalues = eigenvalues,
                   Relative_eig = eigenvalues / sum(eigenvalues),
                   Cumul_eig = cumsum(eigenvalues) / sum(eigenvalues))
    ## Construct results and return it
    result <- list(number.components = number.components,
                   lambda = lambda, 
                   scores = scores,
                   loadings = loadings,
                   values = values,
                   center = ifelse(center, cen, FALSE),
                   scale = ifelse(scale., sc, FALSE))
    class(result) <- c("sparsePca", "pca")
    return(invisible(result))
}


#' S3 Method for extracting scores from pca object
#'
#' @title scores.pca
#' @param pca Required. A \code{\link{pca}} class object
#' @param choices Required. Ordination axes. If missing, returns all axes.
#' @param display Required. Unambiguous abbreviation of 'variables' or 'samples'.
#'                For compatibility with vegan and phyloseq, 'variables' can be replaced
#'                by 'species' or 'taxa' and 'samples' can be replaced by 'sites'.
#' @param ... Optional. Additional argument passed on from scores.default. Unused
#' @return A matrix of scores for axes in 'choices'
#' @seealso \code{\link{pca}}, \code{\link{edgePCA}}, \code{\link{regularizedPca}}, \code{\link{sparsePca}},
scores.pca <- function(pca, choices,
                       display = c("sites", "species", "samples", "variables", "taxa"),
                       ...) {
    display <- match.arg(display)
    if (display %in% c("sites", "samples")) {
        x  <- pca$scores
    }
    if (display %in% c("variables", "taxa", "species")) {
        x <- pca$loadings
    }
    if (!missing(choices)) {
        choices <- choices[choices <= ncol(x)]
        x <- x[ , choices, drop = FALSE]
    }
    return(x)
}

#' Extract supporting taxa from an edge pca
#'
#' @title support_taxa
#' @param physeq Required. A \code{\link{phyloseq-class}} class instance. 
#' @param epca Required. An \code{\link{edgepca}} class object
#' @param choices Optional. Ordination axis for which support taxa are extracted.
#'                Defaults to all axis. 
#' @param threshold Optional. Numeric. Threshold that a taxa loading
#'                  must pass to be extracted. Defaults to 0. 
#' @param number  Optional. Integer. Number of taxa to extract from each axis.
#'                Defaults to all taxa.  
#' @note  If number and threshold are specified, extracts at most \code{number} taxa whose loadings
#'        exceed threshold for each axis.
#' @return A character vector ot taxa names.
support_taxa <- function(physeq, epca, choices = seq(epca$number.components),
                         threshold = 0, number = ntaxa(physeq)) {
    ## Extract loading vectors
    x <- abs(epca$loadings)
    choices <- choices[choices <= ncol(x)]
    x <- x[ , choices, drop = FALSE]
    ## Extract tip loadings
    tree <- phy_tree(physeq)
    pendant.edges <- match(1:length(tree$tip.label) , tree$edge[ , 2])
    x <- x[pendant.edges, , drop = FALSE]
    ## Threshold loadings values
    y <- ifelse(x > threshold, 1, 0)
    x <- x*y
    ## Initialize taxa set
    taxa.set <- rep(FALSE, ntaxa(physeq))
    names(taxa.set) <- taxa_names(physeq)
    ## Extract taxa from each component
    for (comp in 1:ncol(x)) {
        z <- x[ , comp]
        names(z) <- tree$tip.label
        cur.taxa.set <- sort(z, decreasing = TRUE)[1:number]
        cur.taxa.set <- names(cur.taxa.set)[cur.taxa.set > 0]
        taxa.set[cur.taxa.set] <- TRUE
    }
    return(names(taxa.set)[taxa.set])
}


#' Function for plotting a tree with edge fattened and colored according to
#' the loadings of an edgePCA results.  
#'
#' @title plotLoadings
#' @param physeq Required. A \code{\link{phyloseq-class}} class instance. 
#' @param epca Required. An \code{\link{edgepca}} class object
#' @param choices Required. Ordination axes.
#' @param width.lim Optional. Numeric. Minimum and maximal edge after fattening.
#'                  Defaults to c(0.1, 4). Set to c(0, x) for true linear scaling.
#' @param fatten.edges.by Required. Aesthetics that will be used for fattening.
#'                        Subset of 'size' and 'alpha'.
#' @param fatten.tips Optional. Logical. Should tips be fattened like pendant edges.
#'                    Defaults to FALSE. 
#' @param color Optional. Color vector of length 2 used to distinguish positive
#'              loading edges (color[1]) and negative loading edges (color[2]).
#'              Defaults to c("#EF8A62", "#67A9CF")
#' @param color.tip.by Optional. Either a single character string matching a variable
#'              name in the corresponding tax_table of `physeq`, or a factor with
#'              the same length as the number of taxa in `physeq`. If NULL, nothing happens.
#' @param missing.color Optional. Color used for tips and edges with 0 loading.
#'                      Defaults to NULL. If NULL, nothing happens. Use "white", "transparent",
#'                      or par("bg") to remove them from plot.
#' @param ... Additional arguments passed on to plot.phylo
#' 
#' @return Nothing, this function is used for its side effect of plotting a tree
#' @seealso \code{\link{edgePCA}}
plotLoadings <- function(physeq, epca, choices, fatten.by = c("size"), fatten.tips = FALSE, 
                         width.lim = c(0.1, 4), color = c("#EF8A62", "#67A9CF"),
                         color.tip.by = NULL, missing.color = NULL, 
                         ...) {
    ## Check fatten.by
    if (any( !fatten.by %in% c("size", "alpha"))) {
        stop("fatten.by must be a subset of c(\"size\", \"alpha\")")
    }
    ## Extract loading vectors
    x <- epca$loadings
    choices <- choices[choices <= ncol(x)]
    x <- x[ , choices, drop = FALSE]
    ## Extract tree
    tree <- phy_tree(physeq)
    pendant.edges <- match(1:length(tree$tip.label) , tree$edge[ , 2])
    tip.loadings <- x[pendant.edges, , drop = FALSE]
    ## Get edge.color
    edge.color <- ifelse(x > 0, color[1], color[2])
    ## Get tip.color 
    tip.color <- edge.color[pendant.edges, , drop = FALSE]
    ## Get fattening factor
    fattening.factor <- rescale(abs(x), to = width.lim)
    ## Fatten edges by size
    if ("size" %in% fatten.by) {
        edge.width <- fattening.factor
    } else {
        edge.width <- rep(1, length(x))
    }
    dim(edge.width) <- dim(x)
    ## If alpha, fatten edges by alpha
    if ("alpha" %in% fatten.by) {
        scaled.alpha <- fattening.factor / width.lim[2]
        edge.color <- alpha(edge.color, scaled.alpha)
        dim(edge.color) <- dim(x)
    }
    ## Color tips
    if (!is.null(color.tip.by)) {
        tip.color <- color_edges(physeq, group = color.tip.by,
                                tip.only = TRUE, method = "majority")
        legend.palette <- tip.color$palette
        tip.color <- tip.color$tip
        tip.color <- matrix(rep(tip.color, ncol(x)), ncol = ncol(x))
    }
    ## Fatten tips
    if (fatten.tips) {
        tip.width <- edge.width[pendant.edges, , drop = FALSE]
        if ("alpha" %in% fatten.by) {
            dim(scaled.alpha) <- dim(x)
            scaled.alpha <- scaled.alpha[pendant.edges, , drop = FALSE]
            tip.color <- alpha(tip.color, scaled.alpha)
            dim(tip.color) <- dim(scaled.alpha)
        }
    }
    ## Set missing edges and possibly tips to missing.color
    if (!is.null(missing.color)) {
        edge.color[x == 0] <- missing.color
        tip.color[tip.loadings == 0] <- missing.color
    }
    ## Plot tree
    ## Prepare arguments for plot phylo
    args <- list(x = tree, edge.width = 1, edge.color = "black")
    args <- c(args, list(...))
    ## Prepare layout if there are multiple trees
    oldmar <- par("mar")
    n.axis <- ncol(x)
    if (n.axis > 1) {
        layout.nrow <- floor(sqrt(n.axis))
        layout.ncol <- ceiling(n.axis / layout.nrow)
        par(mfcol = c(layout.nrow, layout.ncol), mar = c(0, 0, 1, 0))
        for (i in 1:n.axis) {
            args[["edge.width"]] <- edge.width[ , i]
            args[["edge.color"]] <- edge.color[ , i]
            args[["tip.color"]] <- tip.color[ , i]
            if (fatten.tips) {
                args[["cex"]] <- tip.width[ , i]
            }
            do.call("plot.phylo", args)
            title(paste("Axis", choices[i]))
        }
    } else {
        par(mar = c(0, 0, 1, 0))
        args[["edge.width"]] <- edge.width[ , 1]
        args[["edge.color"]] <- edge.color[ , 1]
        args[["tip.color"]] <- tip.color[ , 1]
        if (fatten.tips) {
            args[["cex"]] <- tip.width[ , 1]
        }
        do.call("plot.phylo", args)
        title(paste("Axis", choices))
    }
    par(mar = oldmar)
}

################################################################################
# Define S3 method extract_eigenvalue function for pca class
# Function is used by `plot_scree` to get the eigenvalue vector from different
# types of ordination objects.
#' @keywords internal
# for pca (pca) objects
extract_eigenvalue.pca = function(ordination) ordination$values$Eigenvalues
################################################################################


#' Method for plotting an edgepca object with enhanced axes. Expands
#' \code{\link{plotLoadings}} and \code{\link{plot_ordination}}
#'
#' @title plot_edgepca
#' @param physeq Required. A \code{\link{phyloseq-class}} object
#' @param epca Required. A \code{\link{edgepca}} class object
#' @param axes Optional. Ordination axes. Defaults to c(1, 2). 
#' @param type Required. Unambiguous abbreviation of 'variables' or 'samples'.
#'                For compatibility with vegan and phyloseq, 'variables' can be replaced
#'                by 'species' or 'taxa' and 'samples' can be replaced by 'sites'.
#' @param args.loadings Optional. Additional passed on to plotLoadings. Defaults to NULL
#' @param widths Optional. A vector of 2 values for the widths of columns. Defaults to c(1, 1)
#' @param heights Optional. A vector of 2 values for the heights of rows. Defaults to c(1, 1)
#' @param ... Optional. Additional argument passed on to plot_ordination
#' @return Nothing. The function is called for its side effect. 
#' @seealso \code{\link{pca}}, \code{\link{edgePCA}}, \code{\link{regularizedPca}}, \code{\link{sparsePca}},
plot_edgepca <- function(physeq, epca, axes = c(1, 2),
                         widths = c(1, 1), heights = c(1, 1),
                         args.loadings = NULL, ...) {
    ## Setup graphic devices
    layout(matrix(c(2, 0, 3, 1), 2, 2), widths = widths, heights = heights)
    ## Plot first loading
    args1 <- list(physeq = physeq, epca = epca, choices = axes[1])
    args.loadings <- c(args1, args.loadings)
    do.call("plotLoadings", args = args.loadings)
    ## Plot second loading
    args.loadings[["choices"]] <- axes[2]
    do.call("plotLoadings", args = args.loadings)
    ## Plot ordination using plot_ordination and ggplot2
    plot.new()
    vps <- baseViewports()
    pushViewport(vps$figure)
    p <- plot_ordination(physeq = physeq, ordination = epca, axes = axes, ...)
    g <- ggplot_gtable(ggplot_build(p))
    grid.draw(g)
    upViewport()
}

