#' @title Performs a LefSe analyses on the data set, based on a grouping variable given either as a sample variable or a factor.
#'
#' @description The purpose of this method is to reproduce the lefse code in R. Briefly:
#' * the taxa counts are clr transformed (after adding a pseudo count of 1)
#' * A univariate Kruskall-Wallis test is performed on each taxa to assess association with the grouping structure
#' * A LDA model is fitted on the significant taxa
#'
#' @param physeq Required. \code{\link{phyloseq-class}} object
#' @param group Required. Grouping variable, used to compute the core microbiome in a groupwise fashion. Either a single character string matching a variable name in the corresponding `sample_data` of `physeq`, or a factor with the same length as the number of samples in `physeq`.
#' @param pseudocount Optional. Pseudocount added to all counts, to avoid problems during clr transformation. Defaults to 1 and automatically replaced by 0 if all counts are between 0 and 1.
#' @param padj.threshold Optional. Adjusted p-value threshold used to include taxa in the LDA model.Defaults to 0.05
#'
#' @return A data frame with three columns
#' * `OTU` : taxa name
#' * `padj`: adjusted p-value (BH correction) of the Kruskall-Wallis test
#' * `LD`  : Matrix or vector of LD scores
#'
#' @export
#'
#' @examples
#' data(food)
#' lefse <- lefse(food, "FoodType", padj.threshold = 0.001)
#' @importFrom MASS lda
#' @importFrom dplyr as_tibble filter left_join pull
#' @importFrom methods as
#' @importFrom phyloseq get_variable otu_table prune_taxa sample_data sample_variables taxa_are_rows taxa_names taxa_sums transform_sample_counts
#' @importFrom stats kruskal.test p.adjust
#' @importFrom tibble tibble
lefse <- function(physeq, group, pseudocount = 1, padj.threshold = 0.05) {
  # Remove empty OTUs
  physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)

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

  ## clr transform of the counts
  if (all(otu_table(physeq) < 1) & all(otu_table(physeq) >= 1)) {
    pseudocount <- 0
    warning("All counts are in (0, 1), suggesting you're working with proportions. Pseudocount automatically set to 0.")
  }
  physeq <- transform_sample_counts(physeq, function(x) {clr(x + pseudocount)})

  ## Robust Kruskall-Wallis test
  counts <- otu_table(physeq) %>% as("matrix")
  if (!taxa_are_rows(physeq)) counts <- t(counts)
  res <- dplyr::tibble(OTU     = taxa_names(physeq),
                       padj    = p.adjust(apply(counts, 1, function(x) kruskal.test(x, g = group)$p.value), method = "BH"))

  ## LDA on differential taxa
  da_otus <- dplyr::filter(res, padj <= padj.threshold) %>%
    # dplyr::arrange(padj) %>%
    dplyr::pull(OTU)
  # da_otus <- da_otus[1:min(length(da.otus), phyloseq::nsamples(physeq))]
  if (length(da_otus) == 0) stop("No significant taxa based on the Kruskal-Wallis test")
  counts <- scale(t(counts[da_otus, ]), center = TRUE, scale = TRUE)
  lda <- MASS::lda(counts, group)

  ## Join with kruskal-Wallis p-values
  dplyr::left_join(res,
                   lda$scaling %>% as_tibble(rownames = "OTU"),
                   by = "OTU")

}
