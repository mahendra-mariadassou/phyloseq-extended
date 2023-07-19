#' Create a \link{biom-class} from a \link{phyloseq-class} object.
#'
#' @param physeq (Required). A \link{phyloseq-class} object.
#' @param biom_format (Optional). Either "frogs" (default) or "standard". Controls the way the taxonomy is written in the \link{biom-class} object. \code{biom_format = "frogs"} is intended for later use with \link{import_frogs} and \code{biom_format = "standard"} for later use with \link{import_biom}.
#' @param rows_metadata (Optional). Either `NULL` (default) or a named list of rows metadata included in the final biom output. The taxonomy and blast_taxonomy fields are automatically ignored and replaced by corresponding fields from the physeq argument
#' @param matrix_type (Optional). Either "dense" (default) to store the otu count table as dense matrix or "sparse" to store it as a sparse matrix.
#'
#' @return A \link{biom-class} object
#' @export
#'
#' @importFrom biomformat make_biom
#' @importFrom methods as
#' @importFrom phyloseq access otu_table taxa_are_rows taxa_names
#' @importFrom purrr map_chr
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
#' phyloseq_to_biom(food, matrix_type = "sparse")
#' }
phyloseq_to_biom <- function(physeq, biom_format = c("frogs", "standard"), rows_metadata = NULL, matrix_type = c("dense", "sparse")) {
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
  ## Replace observation metadata
  if (!is.null(rows_metadata)) {
    names(rows_metadata) <- purrr::map_chr(rows_metadata, "id")
    rows_metadata <- rows_metadata[phyloseq::taxa_names(physeq)]
    biom@.Data[[10]] <- unname(rows_metadata)
  }
  ## Update taxonomy to comply with frogs / standard biom format and use taxonomy from the physeq object
  biom_format <- match.arg(biom_format)
  correct_taxonomy <- function(x) {
    ## Unhappy path
    if (is.null(tdf)) return(x)
    ## Happy path
    current_taxonomy <- unname(tdf[x$id, ])
    if (all(is.na(current_taxonomy))) current_taxonomy <- NA
    if (!is.list(x$metadata)) x$metadata <- list()
    if (biom_format == "standard") x$metadata$taxonomy <- current_taxonomy
    if (biom_format == "frogs")    x$metadata$blast_taxonomy <- current_taxonomy
    x
  }
  biom@.Data[[10]] <- lapply(biom@.Data[[10]], correct_taxonomy)
  ## Converts observations to integers for FROGS output
  if (biom_format == "frogs") biom@.Data[[12]] <- lapply(biom@.Data[[12]], as.integer)
  ## Convert count matrix to sparse format (if matrix_type = "sparse")
  matrix_type <- match.arg(matrix_type)
  if (matrix_type == "sparse") {
    biom@.Data[[7]] <- "sparse"
    non_zero_index <- which(cdf > 0L, arr.ind = TRUE, useNames = FALSE)
    sparse_cdf <- cbind(non_zero_index - 1L, ## 0-zero based coordinates in JSON
                        cdf[non_zero_index])
    biom@.Data[[12]] <- sparse_cdf
  }
  biom
}

#' Export a \link{phyloseq-class} to several text files: biom, newick and fasta.
#'
#' @param physeq (Required). A \link{phyloseq-class} object
#' @param biom_file (Required). A character string indicating the file location of the biom formatted file.
#' @param tree_file (Optional). A character string indicating the file location of the newick formatted file. If \code{NULL} (default), the [phy_tree()] component of \code{physeq} is ignored.
#' @param fasta_file (Optional). A character string indicating the file location of the fasta formatted file. If \code{NULL} (default), the [refseq()] component of \code{physeq} is ignored.
#' @param ... Additional arguments passed on [phyloseq_to_biom()]
#' @inheritParams phyloseq_to_biom
#'
#' @return Nothing. The function is used for its side effect of exporting a \link{phyloseq-class} object to text files.
#' @export
#'
#' @importFrom Biostrings writeXStringSet
#' @importFrom ape write.tree
#' @importFrom phyloseq access phy_tree
#'
#' @examples
#' data(food)
#' tmp_biom <- tempfile()
#' tmp_tree <- tempfile()
#' write_phyloseq(food, biom_file = tmp_biom, tree_file = tmp_tree, matrix_type = "sparse")
#' ## The output biom can be read again as a phyloseq object
#' import_frogs(tmp_biom, tmp_tree)
write_phyloseq <- function(physeq, biom_file, tree_file = NULL, fasta_file = NULL, biom_format = c("frogs", "standard"), ...) {
  biom <- phyloseq_to_biom(physeq, biom_format = match.arg(biom_format), ...)
  ## Write biom using a safe copy of biomformat::write_biom()
  write_biom_safe(biom, biom_file)
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

#' Export a \link{phyloseq-class} to a tsv file
#'
#' @param physeq (Required). A \link{phyloseq-class} object
#' @param tsv_file (Required). A character string indicating the file location of the tsv file.
#'
#' @return Nothing. The function is used for its side effect of exporting a \link{phyloseq-class} object to a tsv file.
#' @export
#'
#' @details The tsv files contains:
#' * taxa names
#' * taxa affiliation (when available)
#' * taxa sequences (when available)
#' * taxa abundances across samples (1 column per sample)
#' * taxa prevalence
#' * taxa mean abundance (conditional upon presence)
#' All metadata pertaining to the samples is lost.
#'
#' @importFrom dplyr as_tibble bind_cols mutate
#' @importFrom methods as
#' @importFrom phyloseq access otu_table taxa_are_rows taxa_names
#' @importFrom tibble tibble
#'
#' @examples
#' data(food)
#' phyloseq_to_tsv(food)
phyloseq_to_tsv <- function(physeq) {
  ## Base tibble
  results <- dplyr::tibble(OTU = phyloseq::taxa_names(physeq))
  ## Add taxonomic affiliations if any
  taxtab <- phyloseq::access(physeq, "tax_table")
  if (!is.null(taxtab)) {
    results <- dplyr::bind_cols(
      results,
      as(taxtab, "matrix") %>% dplyr::as_tibble()
    )
  }
  ## Add sequence
  refseq <- phyloseq::access(physeq, "refseq")
  if (!is.null(refseq)) {
    results <- mutate(results, Sequence = unname(as.character(refseq)))
  }
  ## Add proportion data
  cdf <- phyloseq::otu_table(physeq) %>% as("matrix")
  if (!phyloseq::taxa_are_rows(physeq)) cdf <- t(cdf)
  results <- dplyr::bind_cols(results, dplyr::as_tibble(cdf))
  ## Add prevalence and conditional mean
  results <- dplyr::mutate(results,
                    prevalence    = rowMeans(cdf > 0),
                    avg_abundance = rowMeans(cdf) / prevalence)
  results
}

#' Safe replacement of [biomformat::write_biom()]
#'
#' @description Works for biom data with only one column (sample), unlike the original version.
#'
#' @importFrom jsonlite toJSON
#' @importFrom methods slot
#'
#' @examples
#' data(food)
#' small_food <- food %>% prune_samples(sample_names(.)[1], . ) %>% prune_taxa(taxa_names(.)[1:5], .)
#' small_biom <- small_food %>% phyloseq_to_biom()
#' tmp_biom <- tempfile()
#' write_biom_safe(small_biom, tmp_biom)
#' readLines(tmp_biom)
## safe alternative to biomformat::write_biom
write_biom_safe <- function(x, biom_file) {
  cat(jsonlite::toJSON(list(
    id                     = slot(x, ".Data")[[1]],
    format                 = slot(x, ".Data")[[2]],
    format_url             = slot(x, ".Data")[[3]],
    type                   = slot(x, ".Data")[[4]],
    generated_by           = slot(x, ".Data")[[5]],
    date                   = slot(x, ".Data")[[6]],
    matrix_type            = slot(x, ".Data")[[7]],
    matrix_element_type    = slot(x, ".Data")[[8]],
    shape                  = slot(x, ".Data")[[9]],
    rows                   = slot(x, ".Data")[[10]],
    columns                = slot(x, ".Data")[[11]],
    data                   = lapply(slot(x, ".Data")[[12]], I)
    ), auto_unbox = TRUE),
    file = biom_file)
}
