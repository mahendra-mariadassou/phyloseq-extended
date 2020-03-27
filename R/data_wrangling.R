#' @title Merge samples based on a grouping variable given either as a sample variable or a factor.
#'
#' @description The purpose of this method is to merge/agglomerate the sample indices of a phyloseq object according to a categorical variable contained in a sample_data or a provided factor. Unlike the original \code{\link{merge_samples}}, sample metadata are not coerced to factors. Instead, the metadata from the first sample in each group is kept and the samples are renamed according to the levels of the factors.
#'
#' @inheritParams phyloseq::merge_samples
#' @param fun Either "sum" (default) or "mean". Controls how abundances are merged within each factor level.
#' @param update.names Logical. If `TRUE` (default) merged samples are named according to the levels of the grouping variable. Else, keep the name of the first sample in each group.
#'
#' @return A phyloseq object that has had its sample indices merged according to the factor indicated by the group argument. The output class matches x.
#'
#' @seealso \code{\link{merge_samples}}
#'
#' @export
#'
#' @examples
#' data(food)
#' mg <- merge_group(food, "EnvType")
#' mg
#' ## Note that sample data are preserved
#' sample_data(mg)
merge_group <- function(physeq, group, fun = c("sum", "mean"), update.names = TRUE) {
  fun <- match.arg(fun)

  # Build grouping factor
  if (!is.null(sample_data(physeq, FALSE))) {
    if (class(group) == "character" & length(group) == 1) {
      if (!group %in% sample_variables(physeq)) {
        stop("group not found among sample variable names.")
      }
      group <- get_variable(physeq, group)
    }
  }
  if (class(group) != "factor") {
    group <- factor(group)
  }
  # Group abundances
  cdf <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) cdf <- t(cdf)
  cdf_merged <- rowsum(cdf, group, reorder = FALSE)
  if (fun == "mean") {
    scaling <- table(group)[rownames(cdf_merged)]
    cdf_merged <- cdf_merged / matrix(scaling[row(cdf_merged)], ncol = ncol(cdf_merged))
  }
  ## Rename merged samples according to first sample in group
  rownames(cdf_merged) <- sample_names(physeq)[!duplicated(group)]

  ## Update otu_table
  otu_table(physeq) <- otu_table(cdf_merged, taxa_are_rows = FALSE)

  ## Update sample names
  if (update.names) {
    sample_names(physeq) <- unique(group)
  }

  physeq
}

#' Fast alternative to [tax_glom()]
#'
#' @param physeq Required. \code{\link{phyloseq-class}} object
#' @param taxrank A character string specifying the taxonomic level that you want to agglomerate over. Should be among the results of \code{rank_names(physeq)}. The default value is \code{rank_names(physeq)[1]}, which may agglomerate too broadly for a given experiment. You are strongly encouraged to try different values for this argument.
#' @param bad_empty (Optional). Character vector. Default to empty white spaces and tabs. Defines the bad/empty values that should be replaced by "Unknown".
#'
#' @return A taxonomically-agglomerated \code{\link{phyloseq-class}} object.
#' @export
#'
#' @details fast_tax_glom differs from [tax_glom()] in three important ways:
#' * It is based on dplyr function and thus much faster
#' * It does not preserve the phy_tree slot
#' * It only preserves taxonomic ranks up to \code{taxrank} and discard all other ones.
#' The "archetype" (OTU name) and sequence for each group are chosen as those from the most abundant taxa
#' within that group (for compatibility with [tax_glom()])
#'
#' @seealso [tax_glom()], [merge_taxa()]
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr as_tibble mutate group_by_at arrange slice select
#' @examples
#' data(food)
#' fast_tax_glom(food, "Species")
fast_tax_glom <- function(physeq, taxrank = rank_names(physeq)[1], bad_empty = c("", " ", "\t")) {
  rank_number <- match(taxrank, rank_names(physeq))
  if (is.na(rank_number)) {
    stop("Bad taxrank argument. Must be among the values of rank_names(physeq)")
  }
  ranks <- rank_names(physeq)[1:rank_number]
  ## change NA and empty to Unknown before merging
  tax <- as(tax_table(physeq), "matrix")
  tax[is.na(tax) | tax %in% bad_empty] <- "Unknown"
  ## create groups
  tax <- tax[ , ranks, drop = FALSE] %>%
    dplyr::as_tibble(tax, .name_repair = "minimal") %>%
    dplyr::mutate(Abundance = taxa_sums(physeq),
                  archetype = taxa_names(physeq)) %>%
    dplyr::group_by_at(vars(one_of(ranks))) %>%
    dplyr::mutate(group = group_indices())
  ## create new_taxonomy
  new_tax <- tax %>%
    dplyr::arrange(desc(Abundance)) %>%
    dplyr::slice(1) %>% dplyr::arrange(group) %>%
    dplyr::select(-group, -Abundance) %>%
    tibble::column_to_rownames(var = "archetype") %>% as.matrix()
  ## create new count table
  otutab <- otu_table(physeq)
  if (!taxa_are_rows(physeq)) otutab  <- t(otutab)
  otutab <- rowsum(otutab, group = tax$group, reorder = TRUE)
  rownames(otutab) <- rownames(new_tax)
  ## create new refseq
  seqs <- access(physeq, "refseq")
  if (!is.null(seqs)) seqs <- seqs[rownames(new_tax)]
  ## return merged phyloseq
  phyloseq(sample_data(physeq),
           tax_table(new_tax),
           otu_table(otutab, taxa_are_rows = TRUE),
           seqs
  )
}

