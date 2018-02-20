## Estimate mean niche value of taxa
compute_niche <- function(physeq, niche, freq = TRUE) {
  ## Args
  ## - physeq: phyloseq class object, otu abundances are extracted from this object
  ## - niche:  Either the a single character string matching a
  ##           variable name in the corresponding sample_data of ‘physeq’, or a
  ##           numeric vector with the same length as the number of samples in ‘physeq’.
  ## - freq    Logical. Should counts be replaced by frequencies (TRUE) or kept as such (FALSE).
  ##           Defaults to TRUE. TRUE corresponding to weighting all samples equally while FALSE
  ##           weights them with their library size. 
  ## 
  ## Returns;
  ## A named vector with the mean niche of all taxa. If a taxa is present in samples 1, ..., n with
  ##  relative abundances a[1], ..., a[n] and the samples have niche values v[1], ..., v[n], the mean niche is
  ##  given by sum(a * v)/sum(a). 
  ##
  ## Note:
  ## - all 0 vectors lead to undefined (NA) niche value
  ## 
  ## Returns
  ## - niche.value: named vector of mean niche value
  
  ## Get sample niche vector 
  if (!is.null(sample_data(physeq, FALSE))) {
    if (class(niche) == "character" & length(niche) == 1) {
      if (! niche %in% sample_variables(physeq)) {
        stop("niche not found among sample variable names.")
      }
      niche <- get_variable(physeq, niche)
    }
  }
  if (class(niche) != "numeric") {
    warning("Niche variable is not numeric and was coerced to numeric, results may be meaningless.")
    niche <- as.numeric(niche)
  }
  
  ## Construct relative abundances by sample
  tdf <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) { tdf <- t(tdf) }
  if (freq) {
    tdf <- apply(tdf, 2, function(x) x/sum(x))
  }
  
  ## Create niche matrix M where M[, j] is niche value in sample j (all columns are proportionals to 1)
  niche.matrix <- matrix(rep(niche, each = ntaxa(physeq)), 
                         nrow = ntaxa(physeq))
  ## Set abundance to 0 in sample where the niche is missing and 
  ## turn relative abundances to (taxa-wise) proportions
  tdf <- tdf * !is.na(niche.matrix)
  total.abundance <- rowSums(tdf)
  absent.taxa <- (total.abundance == 0)
  tdf <- tdf / total.abundance
  
  ## Average niche over all samples (samples with NA do not contribute)
  niche.value <- rowSums(tdf * niche.matrix, na.rm = TRUE)
  
  ## Add names, remove NaN/NA and print a warning
  if (any(absent.taxa)) {
    warning(paste0(sum(absent.taxa), " taxa are absent from your dataset or present only in samples for which niche is missing, their mean niche value has been set to NA.\nConsider filtering them out."))
    niche.value[absent.taxa] <- NA
  }
  names(niche.value) <- taxa_names(physeq)
  
  return(niche.value)
}