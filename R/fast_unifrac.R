#' Compute unifrac distances
#'
#' @inheritParams phyloseq::UniFrac
#' @param ... Arguments passed from other methods. Used for compatibility purpose.
#'
#' @return A distance matrix
#' @export
#'
#' @importFrom phyloseq access taxa_sums prune_taxa
#' @importFrom ape is.rooted
#'
#' @examples
#' data(food)
#' unifrac(food)
UniFrac <- function(physeq, weighted = FALSE, normalized = TRUE, ...) {
  tree <- phyloseq::access(physeq, "phy_tree")
  if (is.null(tree)) {
    stop("UniFrac distances require a phylogenetic tree.")
  }
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths. See tree$edge.length. Cannot compute UniFrac without branch lengths")
  }
  if (!ape::is.rooted(tree)) {
    stop("Tree is not rooted. Make sure your tree is rooted before attempting UniFrac calculation. See ?ape::root")
  }
  ## Remove taxa with null abundances
  if (any(x <- phyloseq::taxa_sums(physeq) > 0)) {
    physeq <- phyloseq::prune_taxa(x > 0, physeq)
  }
  fastUniFrac(physeq, weighted, normalized)
}


#' @importFrom phyloseq phy_tree nsamples sample_names taxa_are_rows otu_table
#' @importFrom ape reorder.phylo
#'
fastUniFrac <- function(physeq, weighted, normalized) {
  ## Extract components and order in pruning wise order
  tree    <- phyloseq::phy_tree(physeq) %>% ape::reorder.phylo("postorder")
  n_tips  <- length(tree$tip.label)
  n_nodes <- tree$Nnode

  ## Create counts matrix (samples in line, taxa in columns)
  counts <- matrix(data = 0,
                   nrow = phyloseq::nsamples(physeq), ncol = n_tips + n_nodes,
                   dimnames = list(phyloseq::sample_names(physeq)))

  ## Fill `counts` with counts descending from terminal branches
  if (phyloseq::taxa_are_rows(physeq)) {
    counts[, 1:n_tips] <- phyloseq::otu_table(physeq) %>%
      as("matrix") %>% t() %>% `[`(, tree$tip.label)
  } else {
    counts[, 1:n_tips] <- phyloseq::otu_table(physeq) %>%
      as("matrix") %>% `[`(, tree$tip.label)
  }
  ## Change counts to fractions
  counts[, ] <- counts / rowSums(counts)
  ## Compute fractions descending from inner branches using pruning algorithm
  for (i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child  <- tree$edge[i, 2]
    counts[, parent] <- counts[, parent] + counts[, child]
  }
  ## reorder counts to be in the same order as tree$edge.length
  ## using the fact that the root branch (first branch) corresponds to node n_tips + 1
  counts[, ] <- counts[, c(n_tips + 1, tree$edge[, 2])]
  if (!weighted) {
   counts[ , ] <- 0.0 + (counts > 0)
  }

  ## Weight the fractions by corresponding branch lengths,
  ## considering root branch to have length 0
  edge_length <- c(0, tree$edge.length)
  ## equivalent to but faster than `counts %*% diag(edge_length)`
  ## Multiply each row (fractions) by edge lengthes
  counts[ , ] <- counts * edge_length[col(counts)]

  ## Unnormalized unifrac distance is then simply the Manhattan distance
  ## between the modified vectors of all samples
  raw_unifrac <- dist(counts, method = "manhattan")
  if (!normalized) return(raw_unifrac)

  ## Compute normalization constants
  if (weighted) {
    stl <- rowSums(counts) ## Total length of each sample
    max_unifrac <- outer(stl, stl, FUN = "+") %>% as.dist()
    return(raw_unifrac / max_unifrac)
  } else {
    ## If a = phylogenetic length of sample 1, b of sample 2 and c = distinct length between sample 1 and sample 2
    ## then total length covered by samples 1 and 2 is
    ## d = a + b - (a + b - c)/2 = (a + b + c)/2
    stl <- rowSums(counts)
    max_unifrac <- (as.dist(outer(stl, stl, FUN = "+")) + raw_unifrac)/2
    return(raw_unifrac / max_unifrac)
  }

}
