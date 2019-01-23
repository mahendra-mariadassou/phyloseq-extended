filter_abundance <- function(physeq, frac = 0.001, A = 0.05 * nsamples(physeq)) {
  test_function <- function(x) { x >= frac}
  tokeep <- physeq %>% transform_sample_counts(function(x) { x / sum(x) }) %>% genefilter_sample(test_function, A = A)
  return(prune_taxa(tokeep, physeq))
}


## filter OTU with prevalence higher than 50% in each level of factor group 
## and highest relative abundance higher than 0.1% 
filter_phyloseq <- function(physeq, prev.thresh = 0.5, abund.thresh = 0.001, group = rep(1, nsamples(physeq))) {
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