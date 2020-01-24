#' Create a \link{biom-class} from a \link{phyloseq-class} object.
#'
#' @param physeq (Required). A \link{phyloseq-class} object.
#' @param biom_format (Optional). Either "frogs" (default) or "standard". Controls the way the taxonomy is written in the \link{biom-class} object. \code{biom_format = "frogs"} is intended for later use with \link{import_frogs} and \code{biom_format = "standard"} for later use with \link{import_biom}.
#'
#' @return A \link{biom-class}
#' @export
#'
#' @importFrom phyloseq otu_table access
#' @importFrom biomformat make_biom
#'
#' @seealso \link{make_biom}, \link{write_phyloseq}
#'
#' @details This functions differs from \link{make_biom} in the way it encodes the taxonomy
#' in the 'rows' metadata. The difference is intended for better compatibility with phyloseq
#' \code{import_*} class of functions.
#'
#' @examples
#' \dontrun{
#' data(food)
#' phyloseq_to_biom(food)
#' }
phyloseq_to_biom <- function(physeq, biom_format = c("frogs", "standard")) {
  ## Counts
  cdf <- phyloseq::otu_table(physeq) %>% as("matrix")
  if (!phyloseq::taxa_are_rows(physeq)) cdf <- t(cdf)
  ## Sample and observation metadata
  sdf <- phyloseq::access(physeq, "sam_data")
  if (!is.null(sdf)) sdf <- as(sdf, "data.frame")
  tdf <- phyloseq::access(physeq, "tax_table")
  if (!is.null(tdf)) tdf <- as(tdf, "matrix")
  ## Make simple biom
  biom <- biomformat::make_biom(data = cdf, sample_metadata = sdf, observation_metadata = tdf)
  ## Transform taxonomy to comply with frogs / standard biom format
  biom_format <- match.arg(biom_format)
  correct_taxonomy <- function(x) {
    if (biom_format == "frogs")    x$metadata <- list(blast_taxonomy = x$metadata)
    if (biom_format == "standard") x$metadata <- list(taxonomy = x$metadata)
    x
  }
  biom@.Data[[10]] <- lapply(biom@.Data[[10]], correct_taxonomy)
  biom
}



#' Export a \link{phyloseq-class} to several text files: biom, newick and fasta.
#'
#' @param physeq (Required). A \link{phyloseq-class} object
#' @param biom_file (Required). A character string indicating the file location of the biom formatted file.
#' @param tree_file (Optional). A character string indicating the file location of the newick formatted file. If \code{NULL} (default), the [phy_tree()] component of \code{physeq} is ignored.
#' @param fasta_file (Optional). A character string indicating the file location of the fasta formatted file. If \code{NULL} (default), the [refseq()] component of \code{physeq} is ignored.
#' @inheritParams phyloseq_to_biom
#'
#' @return Nothing. The function is used for its side effect of exporting a \link{phyloseq-class} object to text files.
#' @export
#'
#' @importFrom phyloseq access
#' @importFrom biomformat write_biom
#' @importFrom Biostrings writeXStringSet
#'
#' @examples
#' data(food)
#' tmp_biom <- tempfile()
#' tmp_tree <- tempfile()
#' write_phyloseq(food, biom_file = tmp_biom, tree_file = tmp_tree)
#' ## The output biom can be read again as a phyloseq object
#' import_frogs(tmp_biom, tmp_tree)
write_phyloseq <- function(physeq, biom_file, tree_file = NULL, fasta_file = NULL, biom_format = c("frogs", "standard")) {
  biom <- phyloseq_to_biom(physeq, biom_format = match.arg(biom_format))
  ## Write biom
  biomformat::write_biom(x = biom, biom_file = biom_file)
  ## Write fasta (if any)
  refseq <- phyloseq::access(physeq, "refseq")
  if (!is.null(refseq)) {
    if (is.null(fasta_file)) {
      warning("No fasta file provided, skipping export of refseq() component")
    } else{
      Biostrings::writeXStringSet(refseq, filepath = fasta_file, format = "fasta")
    }
  }
  ## Write tree (if any)
  tree <- phyloseq::access(physeq, "phy_tree")
  if (!is.null(tree)) {
    if (is.null(tree_file)) {
      warning("No tree file provided, skipping export of phy_seq() component")
    } else{
      ape::write.tree(phyloseq::phy_tree(physeq), file = tree_file)
    }
  }
}
