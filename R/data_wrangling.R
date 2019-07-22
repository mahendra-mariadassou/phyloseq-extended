#' @title Merge samples based on a grouping variable given either as a sample variable or a factor.
#'
#' @description The purpose of this method is to merge/agglomerate the sample indices of a phyloseq object according to a categorical variable contained in a sample_data or a provided factor. Unlike the original \code{\link{phyloseq::merge_samples}}, sample metada are not coerced to factors. Instead, the metadata from the first sample in each group is kept.
#'
#' @inheritParams phyloseq::merge_samples
#'
#' @return A phyloseq object that has had its sample indices merged according to the factor indicated by the group argument. The output class matches x.
#'
#' @seealso \code{\link{phyloseq::merge_samples}}
#'
#' @export
#'
#' @examples
#' data(food)
#' merge_group(food, "EnvType")
#'
merge_group <- function(physeq, group, fun = mean) {
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
  ## Rename merged samples according to first sample in group
  rownames(cdf_merged) <- sample_names(physeq)[!duplicated(group)]

  ## Update otu_table
  otu_table(physeq) <- otu_table(cdf_merged, taxa_are_rows = FALSE)
  physeq
}
