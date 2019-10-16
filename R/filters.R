filter_abundance <- function(physeq, frac = 0.001, A = 0.05 * nsamples(physeq)) {
  test_function <- function(x) { x >= frac}
  tokeep <- physeq %>% transform_sample_counts(function(x) { x / sum(x) }) %>% genefilter_sample(test_function, A = A)
  return(prune_taxa(tokeep, physeq))
}

#' Prevalence and abundance filters, applied independently on each group
#'
#' @param physeq A \link{phyloseq} class object
#' @param prev.thresh Minimum prevalence in at least one group level for a taxa to be kept
#' @param abund.thresh Minimum relative abundance in at least one sample
#' @param group Either the a single character string matching a variable name in
#'              the corresponding sample_data of `physeq`, or a factor with the same
#'              length as the number of samples in `physeq`.
#' @param rarefy Logical. Should samples be rarefied before abundances and prevalence are computed?
#'
#' @details Apply prevalence and abundance based filters to a \link{phyloseq} class object. The relative abundance threshold needs to be satisfied in at least one sample (i.e. abundances are not computed per group). The prevalence needs to be satisfied in at least one group level (or in the global dataset if no group is provided). The samples are rarefied by default to avoid depth-induced biases on discovery rates for the prevalence threshold.
#'
#' @return Vector of taxa that pass the filters, for further use with \link{prune_sample}.
#'
#' @export
#'
#' @examples
#' data(food) ## 507 taxa
#' ## Taxa with prevalence > 25% in any food type and abundance > 1% in at least one sample.
#' taxa_to_keep <- filter_phyloseq(food, 0.25, 0.01, group = "EnvType")
#' prune_taxa(taxa_to_keep, food) ## 217 taxa
filter_phyloseq <- function(physeq, prev.thresh = 0.5, abund.thresh = 0.001, group = rep(1, nsamples(physeq)), rarefy = TRUE) {
  physeq <- suppressMessages(physeq %>% rarefy_even_depth(trimOTUs = FALSE, rngseed = 20171201))
  ## Abundant otus
  test_function_abundance <- function(x) { (x / sum(x)) >= abund.thresh }
  abundant.otus <- taxa_names(physeq)[genefilter_sample(physeq, test_function_abundance, A = 1)]
  ## Prevalent otus
  prevalent.otus <- estimate_prevalence(physeq, group) %>%
    group_by(otu) %>% summarize(prevalence = max(prevalence)) %>%
    filter(prevalence >= prev.thresh) %>% `[[`("otu")
  tokeep <- intersect(abundant.otus, prevalent.otus)
  return(tokeep)
}
