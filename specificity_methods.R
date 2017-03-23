## Method dedicated to specificity index. Computes, tests and represents several specificity index (simpson, shannon, indicator species and Yanai's index)

## Compute specificity of a vector
## x represents (possibly not normalised) distribution of quantity
## of interest across conditions
specificity <- function(x, index = c("shannon", "simpson", "yanai", "indspec"), groupfrac = NULL) {
    ## Args:
    ## - x: matrix (vectors are automatically coerced to 1-row matrices), specificity
    ##      is computed along the rows of x
    ## - index: method used for computing specificity. yanai refers to (Yanai 200x) and Indspec to
    ##          Indicator Species (Dufrene & Legendre, 1997)
    ## - groupfrac: group fractions in which otus appear, should be the same size as x, required for indspec
    ##
    ## Note:
    ## - all 0 vectors lead to minimum specificity (unlike in 'diversity' of package 'vegan' where
    ##   they lead to minimum diversity which in turn leads to maximum specificity)
    ## 
    ## Returns
    ## - Div: Specificity vector
    index <- match.arg(index)
    if (index == "indspec") { stopifnot(!is.null(groupfrac))} 
    ## Change vector to matrix
    if (is.vector(x)) {
        x <- matrix(x, nrow = 1)
        if (index == "indspec") { groupfrac <- matrix(groupfrac, nrow = 1) }
    }
    ## change x so that empty sites return minimum specificity 
    total <- rowSums(x)
    x[total == 0, ] <- 1
    ## Compute diversity index
    if (index == "yanai") {
        x <- 1 - (x / rowMax(x))
        Div <- rowSums(x, na.rm =TRUE)/(ncol(x) - 1)
    } else { ## Normalize data 
        x <- x / rowSums(x)
        if (index == "shannon") {
            x <- -x * log(x)
            Div <- apply(x, 1, sum, na.rm = TRUE)
            Div <- 1/exp(Div)
        }
        if (index == "simpson") {
            x <- x * x
            Div <- apply(x, 1, sum, na.rm = TRUE)
        }
        if (index == "indspec") {
            x <- x * groupfrac
            Div <- apply(x, 1, max)
        }
    }
    return(Div)
}


## Compute local specificity of a vector
## x represents (possibly not normalised) distribution of quantity
## of interest across conditions
local_specificity <- function(x, index = c("fraction", "indspec"), groupfrac = NULL) {
    ## Args:
    ## - x: matrix (vectors are automatically coerced to 1-row matrices), specificity
    ##      is computed along the rows of x
    ## - index: method used for computing specificity. "fraction" corresponds to fraction
    ##          of overall total found in a given condition and "indspec" to "fraction" weighted
    ##          by the prevalence of quantity in samples from condition (found in groupfrac)
    ## - groupfrac: group fractions in which otus appear, should be the same size as x,
    ##              required for indspec
    ##
    ## Note:
    ## - all 0 vectors lead to minimum specificity (unlike in 'diversity' of package 'vegan' where
    ##   they lead to minimum diversity which in turn leads to maximum specificity)
    ## 
    ## Returns
    ## - Div: Specificity vector
    index <- match.arg(index)
    if (index == "indspec") { stopifnot(!is.null(groupfrac))} 
    ## Change vector to matrix
    if (is.vector(x)) {
        x <- matrix(x, nrow = 1)
        if (index == "indspec") { groupfrac <- matrix(groupfrac, nrow = 1) }
    }
    ## change x so that empty sites return minimum specificity 
    total <- rowSums(x)
    x[total == 0, ] <- 1
    ## Compute local specificity index
    Div <- x / rowSums(x, na.rm = TRUE)
    if (index == "indspec") {
        Div <- Div * groupfrac
    }
    return(Div)
}


## Estimate specificity of an otu to a given factor
estimate_specificity <- function(physeq, group,
                                 index = c("shannon", "simpson", "yanai", "indspec"), freq = TRUE,
                                 B = 999, se = TRUE, parallel = FALSE) {
  ## Args:
  ## - physeq: phyloseq class object, otu abundances are extracted from this object
  ## - index:  method used for computing specificity
  ## - group:  Either the a single character string matching a
  ##           variable name in the corresponding sample_data of ‘physeq’, or a
  ##           factor with the same length as the number of samples in ‘physeq’.
  ## - freq    Logical. Should counts be replaced by frequencies (TRUE) or kept as such (FALSE).
  ##           Defaults to TRUE. TRUE corresponding to weighting all samples equally while FALSE
  ##           weights them with their library size. 
  ## - se:     Logical, variability of the computed specificities be assessed
  ##           (by bootstrapping samples within levels of factor). If TRUE, quantiles
  ##           5, 25, 50, 75 and 95 of specificity coefficients are returned. Defaults to TRUE
  ## - B:      Only used if se is TRUE, number of bootstrap replicates used
  ##           to compute quantiles.
  ## - parallel: Logical, should computations be performed in parallel. Use mclapply and
  ##             parallel HPC framework
  ## 
  ## Returns;
  ## data frame with components
  ## - specificity: observed specificity
  ## - level: factor level in which otu is most abundant
  ## - abundance: overall relative abundance of otu (all levels are weighted equally)
  ## - mean, sd, quantiles 5, 25, 50, 75 and 95% if 'se' is TRUE
  ## Specificity method is used as attribute 'index'

  ## Get grouping factor 
  if (!is.null(sample_data(physeq, FALSE))) {
    if (class(group) == "character" & length(group) == 1) {
      if (! group %in% sample_variables(physeq)) {
        stop("group not found among sample variable names.")
      }
      group <- get_variable(physeq, group)
    }
  }
  if (class(group) != "factor") {
    group <- factor(group)
  }
  
  ## Construct relative abundances by sample
  tdf <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) { tdf <- t(tdf) }
  if (freq) {
      tdf <- apply(tdf, 2, function(x) x/sum(x))
  }
  
  ## Specificity after pooling by groups for one sample
  ## Change the default settings so that presence in no sample is
  ## considered as perfect evenness (instead of specificity, as is the case now)
  ## Per group averages are computed using rowsum and group size (instead of aggregate)
  ## for speed
  spec <- function(x, index, group) {
      meandf <- t(rowsum(t(x), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(x)),
                                          nrow = nrow(x))
      frac <- t(rowsum(t(0 + (x > 0)), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(x)),
                                            nrow = nrow(x))
      return(specificity(meandf, index, frac))
  }

  ## Random stratified sampling
  stratified_sampling<-function(df, group) {
    ## df is the data to sample from
    ## group is the factor vector used to group samples
    ## Order the data based on the groups
    groupContent <- split(1:ncol(df), f = group)
    groupSample <- unlist(lapply(groupContent, function(x) sample(x, length(x), TRUE)))
    return(list(samples = groupSample, group = group[order(group)]))
  }

  ## One sample specificity
  one_rand_sample_spec <- function(df, index, group) {
    strat_samp <- stratified_sampling(df, group)
    return(spec(x = df[ , strat_samp$samples], index = index, group = strat_samp$group))
  }

  ## Compute original sample diversity, add predominant level and overall abundance for each otu
  meandf <- t(rowsum(t(tdf), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(tdf)),
                                      nrow = nrow(tdf))
  domLevel <- colnames(meandf)[apply(meandf, 1, which.max)]
  res <- data.frame(otu = rownames(meandf),
                    specificity = spec(tdf, index, group),
                    level = domLevel,
                    abundance = rowMeans(meandf))
  
  ## Compute variability of specificity indexes
  if (se) {
      ## Replicate specificity estimates (optionally, in parallel using foreach)
      cat("Estimating se, may take a few minutes", sep = "\n")
      if (parallel) {
          resmat <- mclapply(1:B, function(i) {
              cat(paste("bootstrap sample", i), sep = "\n")
              one_rand_sample_spec(tdf, index, group)} )
          resmat <- do.call(cbind, resmat)
      } else {
          resmat <- replicate(n = B, one_rand_sample_spec(tdf, index, group))
      }
      distrib <- t(apply(resmat, 1, function(x) c(mean = mean(x), sd = sd(x),
                                                  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))))
      colnames(distrib) <- c("mean", "sd", "q5", "q25", "q50", "q75", "q95")
      res <- data.frame(res, distrib)
  }
  ## Add method
  attr(res, "index") <- index
  attr(res, "group") <- group
  return(res)
}


## Estimate local specificity of an otu to a given factor
estimate_local_specificity <- function(physeq, group,
                                 index = c("fraction", "indspec"), freq = TRUE,
                                 B = 999, se = TRUE, parallel = FALSE) {
  ## Args:
  ## - physeq: phyloseq class object, otu abundances are extracted from this object
  ## - index:  method used for computing local specificity
  ## - group:  Either the a single character string matching a
  ##           variable name in the corresponding sample_data of ‘physeq’, or a
  ##           factor with the same length as the number of samples in ‘physeq’.
  ## - freq    Logical. Should counts be replaced by frequencies (TRUE) or kept as such (FALSE).
  ##           Defaults to TRUE. TRUE corresponding to weighting all samples equally while FALSE
  ##           weights them with their library size. 
  ## - se:     Logical, variability of the computed specificities be assessed
  ##           (by bootstrapping samples within levels of factor). If TRUE, quantiles
  ##           5, 25, 50, 75 and 95 of specificity coefficients are returned. Defaults to TRUE
  ## - B:      Only used if se is TRUE, number of bootstrap replicates used
  ##           to compute quantiles.
  ## - parallel: Logical, should computations be performed in parallel. Use mclapply and
  ##             parallel HPC framework
  ## 
  ## Returns;
  ## data frame with components
  ## - specificity: observed local specificity
  ## - group: factor level 
  ## - abundance: local relative abundance of otu (all samples of a level are weighted equally)
  ## - mean, sd, quantiles 5, 25, 50, 75 and 95% if 'se' is TRUE
  ## Specificity method is used as attribute 'index'
  ## Get grouping factor 
  if (!is.null(sample_data(physeq, FALSE))) {
    if (class(group) == "character" & length(group) == 1) {
      x1 <- data.frame(sample_data(physeq))
      if (!group %in% colnames(x1)) {
        stop("group not found among sample variable names.")
      }
      group <- x1[, group]
    }
  }
  if (class(group) != "factor") {
    group <- factor(group)
  }
  ## Construct relative abundances by sample
  tdf <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) { tdf <- t(tdf) }
  if (freq) {
      tdf <- apply(tdf, 2, function(x) x/sum(x))
  }
  ## Specificity after pooling by groups for one sample
  ## Change the default settings so that presence in no sample is
  ## considered as perfect evenness (instead of specificity, as is the case now)
  ## Per group averages are computed using rowsum and group size (instead of aggregate)
  ## for speed
  local_spec <- function(x, index, group) {
      meandf <- t(rowsum(t(x), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(x)),
                                          nrow = nrow(x))
      frac <- t(rowsum(t(0 + (x > 0)), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(x)),
                                            nrow = nrow(x))
      result <- local_specificity(meandf, index, frac)
      return(as.vector(result[observed.otu]))
  }
  ## Random stratified sampling
  stratified_sampling<-function(df, group) {
    ## df is the data to sample from
    ## group is the factor vector used to group samples
    ## Order the data based on the groups
    groupContent <- split(1:ncol(df), f = group)
    groupSample <- unlist(lapply(groupContent, function(x) sample(x, length(x), TRUE)))
    return(list(samples = groupSample, group = group[order(group)]))
  }
  ## One sample specificity
  one_rand_sample_spec <- function(df, index, group) {
    strat_samp <- stratified_sampling(df, group)
    return(local_spec(x = df[ , strat_samp$samples], index = index, group = strat_samp$group))
  }
  ## Compute original sample diversity, add predominant level and overall abundance for each otu
  meandf <- t(rowsum(t(tdf), group, reorder = TRUE))
  relative.abundance <- apply(meandf, 2, function(x) x / sum(x))
  observed.otu <- as.vector(relative.abundance > 0)
  res <- data.frame(otu = rep(rownames(relative.abundance), ncol(relative.abundance)), 
                    level = rep(colnames(relative.abundance), each = nrow(relative.abundance)),
                    abundance = as.vector(relative.abundance))
  res <- res[observed.otu, ]
  res$specificity <- local_spec(tdf, index, group)
  ## Compute variability of specificity indexes
  if (se) {
      ## Replicate specificity estimates (optionally, in parallel using foreach)
      cat("Estimating se, may take a few minutes", sep = "\n")
      if (parallel) {
          resmat <- mclapply(1:B, function(i) {
              cat(paste("bootstrap sample", i), sep = "\n")
              one_rand_sample_spec(tdf, index, group)} )
          resmat <- do.call(cbind, resmat)
      } else {
          resmat <- replicate(n = B, one_rand_sample_spec(tdf, index, group))
      }
      distrib <- t(apply(resmat, 1, function(x) c(mean = mean(x), sd = sd(x),
                                                  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))))
      colnames(distrib) <- c("mean", "sd", "q5", "q25", "q50", "q75", "q95")
      res <- data.frame(res, distrib)
  }
  ## Add method
  attr(res, "index") <- index
  attr(res, "group") <- group
  return(res)
}


## Plot specifity against relative abundance
plot_specificity <- function(specOTUS, y = "specificity", color = "level",
                             se = TRUE, plot = TRUE) {
  p <- ggplot(specOTUS, aes_string(x = "abundance", y = y, color = color)) + geom_point(size = 2)
  p <- p + labs(x = "Overall abundance (log10)", y = y)
  p <- p + scale_x_log10()
  if (se) {
    p <- p + geom_linerange(aes(ymin = q25, ymax = q75), alpha = 0.2)
  }
  ## p <- p + facet_wrap(~level)
  p <- p + theme(strip.text.x = element_text(size = 12),
                 strip.text.y = element_text(size = 12),
                 axis.text.x  = element_text(angle=0, size=12))
  ## p <- p + guides(colour  = guide_legend("Treatment"))
  p <- p + geom_smooth(aes(color = NULL), method = "loess", se = FALSE)
  specmin <- ifelse(attr(specOTUS, "index") %in% c("simpson", "shannon"),
                       1/length(levels(specOTUS$level)),
                       0)
  p <- p + geom_hline(yintercept = specmin, lty = 2, color = "grey40")
  if (plot) {
    plot(p)
  }
  return(invisible(p))
}

## Plot local specifity against relative abundance
plot_local_specificity <- function(specOTUS, y = "specificity", formula = y~x, 
                                   color = "level", method = "loess",
                                   se = TRUE, plot = TRUE) {
  p <- ggplot(specOTUS, aes_string(x = "abundance", y = y, color = color)) + geom_point(size = 2)
  p <- p + labs(x = "Overall abundance (log10)", y = y)
  p <- p + scale_x_log10()
  p <- p + facet_wrap(~level)
  if (se) {
    p <- p + geom_linerange(aes(ymin = q25, ymax = q75), alpha = 0.2)
  }
  ## p <- p + facet_wrap(~level)
  p <- p + theme(strip.text.x = element_text(size = 12),
                 strip.text.y = element_text(size = 12),
                 axis.text.x  = element_text(angle=0, size=12))
  ## p <- p + guides(colour  = guide_legend("Treatment"))
  p <- p + geom_smooth(color = "grey45", method = method, formula = formula, se = FALSE)
  specmin <- ifelse(attr(specOTUS, "index") %in% c("simpson", "shannon"),
                       1/length(levels(specOTUS$level)),
                       0)
  p <- p + geom_hline(yintercept = specmin, lty = 2, color = "grey40")
  if (plot) {
    plot(p)
  }
  return(invisible(p))
}


## Plot specificity rank curve with (optionnally) error bars
plot_specificity_distribution <- function(specOTUS, x = "specificity", color = "level",
                                          se = TRUE, plot =TRUE) {
  ## Order specOTUS by specificity
  specrank <- rank(specOTUS[ , x], ties = "random")
  specOTUS$specrank <- specrank
  specOTUS <- specOTUS[ order(specOTUS$specrank), ]
  ## Plot a specificity/rank curve
  p <- ggplot(specOTUS, aes_string(x = x, y = "specrank", color = color))
  p <- p + labs(x = x, y = "Rank")
  ## p <- p + facet_wrap(~level)
  p <- p + geom_point(size = 2)
  p <- p + theme(axis.text.x  = element_text(angle=0, size=12))
  ## p <- p + guides(colour  = guide_legend("Treatment"))
  p <- p + geom_line(aes(color = NULL))
  if (se) {
      p <- p + geom_errorbarh(aes(xmin = q25, xmax = q75), alpha = 0.2, height = 0)
  }
  specmin <- ifelse(attr(specOTUS, "index") %in% c("simpson", "shannon"),
                    1/length(levels(specOTUS$level)),
                    0)  
  p <- p + geom_vline(xintercept = specmin, lty = 2, color = "grey40")
  if (plot) { 
    plot(p)
  }
  return(invisible(p))
}

## Expected Robust specificity for taxa presents in only
## one replicate of a given treatment condition
expected_robust_specificity <- function(group,
                                        index = c("shannon", "simpson", "yanai", "indspec")) {
  prop <- function(n) { return(1 - (1 - 1/n)^n) }
  index <- match.arg(index)
  specMinMax <- function(x, index) {
    res <- switch(index,
                  shannon = cbind(1/x, rep(1, length(x))),
                  simpson = cbind(1/x, rep(1, length(x))),
                  yanai   = cbind(rep(0, length(x)), rep(1, length(x))),
                  indspec = cbind(0, rep(1, length(x))))
    return(res)
  }
  propGroups <- prop(table(group))
  propGroups <- cbind(1 - propGroups, propGroups)
  specRange <- specMinMax(table(group), index)
  res = rowSums(specRange*propGroups)
  return(data.frame(level = names(res), specificity = res))
}


## Assess significance of observed specificities using permutation tests
test_specificity <- function(physeq, group,
                             index = c("shannon", "simpson", "yanai", "indspec"),
                             freq = TRUE, replace = FALSE, 
                             B = 999, parallel = FALSE) {
  ## Args:
  ## - physeq: phyloseq class object, otu abundances are extracted from this object
  ## - index:  method used for computing specificity
  ## - group:  Either the a single character string matching a
  ##           variable name in the corresponding sample_data of ‘physeq’, or a
  ##           factor with the same length as the number of samples in ‘physeq’.
  ## - freq    Logical. Should counts be replaced by frequencies (TRUE) or kept as such (FALSE).
  ##           Defaults to TRUE. TRUE corresponding to weighting all samples equally while FALSE
  ##           weights them with their library size.
  ## - replace Logical. Should samples be replaced when randomizing samples within group levels. 
  ## - B:      Number of permutation tests used to assess significance.
  ## 
  ## Returns;
  ## data frame with components
  ## - specificity: observed specificity
  ## - level: factor level in which otu is most abundant
  ## - abundance: overall relative abundance of otu (all levels are weighted equally)
  ## - rawp:  raw p-value
  ## - adjp: adjusted p-value (computed using "fdr")

  ## Get grouping factor
  if (!is.null(sample_data(physeq, FALSE))) {
    if (class(group) == "character" & length(group) == 1) {
      x1 <- data.frame(sample_data(physeq))
      if (!group %in% colnames(x1)) {
        stop("group not found among sample variable names.")
      }
      group <- x1[, group]
    }
  }
  if (class(group) != "factor") {
    group <- factor(group)
  }
  
  ## Construct relative abundances by sample
  tdf <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) { tdf <- t(tdf) }
  if (freq) {
      tdf <- apply(tdf, 2, function(x) x/sum(x))
  }
  
  ## Specificity after pooling by groups for one sample
  ## Per group averages are computed using rowsum and group size (instead of aggregate)
  ## for speed
  spec <- function(x, index, group) {
      meandf <- t(rowsum(t(x), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(x)),
                                          nrow = nrow(x))
      frac <- t(rowsum(t(0 + (x > 0)), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(x)),
                                            nrow = nrow(x))
      return(specificity(meandf, index, frac))
  }

  ## One sample specificity (not stratified)
  one_rand_sample_spec <- function(df, index, group) {
    return(spec(df[ , sample.int(length(group), replace = replace)], index, group))
  }

  ## Compute original sample diversity, add predominant level and overall abundance for each otu
  meandf <- t(rowsum(t(tdf), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(tdf)),
                                      nrow = nrow(tdf))
  domLevel <- colnames(meandf)[apply(meandf, 1, which.max)]
  res <- data.frame(specificity = spec(tdf, index, group),
                    otu = rownames(meandf),
                    level = domLevel,
                    abundance = rowMeans(meandf),
                    row.names = row.names(tdf))
  ## Replicate specificity estimates
  cat("Estimating p-values, may take a few minutes", sep = "\n")
  if (parallel) {
      resmat <- mclapply(1:B, function(i) {
          cat(paste("bootstrap sample", i), sep = "\n")
          one_rand_sample_spec(tdf, index, group)} )
      resmat <- do.call(cbind, resmat)
      resmat <- cbind(resmat, res$specificity)
  } else {
      resmat <- cbind(replicate(n = B, one_rand_sample_spec(tdf, index, group)), res$specificity)
  }
  distrib <- t(apply(resmat, 1, function(x) c(mean = mean(x), sd = sd(x), min = min(x),
                                              quantile(x, probs = c(0.5, 0.75, 0.9, 0.95, 0.99)))))
  colnames(distrib) <- c("bmean", "bsd", "bmin", "bq50", "bq75", "bq90", "bq95", "bq99")
  rawp <- rowMeans(sweep(resmat, 1, res$specificity, ">="))
  res <- data.frame(res, rawp = rawp, adjp = p.adjust(rawp, method = "fdr"), distrib)
  ## Add method
  attr(res, "index") <- index
  attr(res, "group") <- group
  return(res)
}


## Assess significance of observed local specificities using permutation tests
test_local_specificity <- function(physeq, group,
                             index = c("fraction", "indspec"),
                             freq = TRUE, replace = FALSE, type = c("local", "global"),
                             B = 999, parallel = FALSE) {
  ## Args:
  ## - physeq: phyloseq class object, otu abundances are extracted from this object
  ## - index:  method used for computing local specificity
  ## - group:  Either the a single character string matching a
  ##           variable name in the corresponding sample_data of ‘physeq’, or a
  ##           factor with the same length as the number of samples in ‘physeq’.
  ## - freq    Logical. Should counts be replaced by frequencies (TRUE) or kept as such (FALSE).
  ##           Defaults to TRUE. TRUE corresponding to weighting all samples equally while FALSE
  ##           weights them with their library size.
  ## - type    Unambiguous abbreviation of "global" or "local". If "local", specificity of a taxa
  ##           to an environment is compared to random specificities in that environment only. If
  ##           "global", it is compared to the maximum specificity over all environments. "global"
  ##            is more conservative than "local". Defaults to "local". 
  ## - replace Logical. Should samples be replaced when randomizing samples within group levels. 
  ## - B:      Number of permutation tests used to assess significance.
  ## 
  ## Returns;
  ## data frame with components
  ## - specificity: observed specificity
  ## - level: grouping factor level
  ## - abundance: local relative abundance of otu (all levels are weighted equally)
  ## - rawp:  raw p-value
  ## - adjp: adjusted p-value (computed using "fdr")

  ## Get grouping factor
  if (!is.null(sample_data(physeq, FALSE))) {
    if (class(group) == "character" & length(group) == 1) {
      x1 <- data.frame(sample_data(physeq))
      if (!group %in% colnames(x1)) {
        stop("group not found among sample variable names.")
      }
      group <- x1[, group]
    }
  }
  if (class(group) != "factor") {
    group <- factor(group)
  }
  
  ## Construct relative abundances by sample
  tdf <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) { tdf <- t(tdf) }
  if (freq) {
      tdf <- apply(tdf, 2, function(x) x/sum(x))
  }

  ## Get type
  type <- match.arg(type)
  
  ## Local Specificity after pooling by groups for one sample
  ## Per group averages are computed using rowsum and group size (instead of aggregate)
  ## for speed
  local_spec <- function(x, index, group, type) {
      meandf <- t(rowsum(t(x), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(x)),
                                          nrow = nrow(x))
      frac <- t(rowsum(t(0 + (x > 0)), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(x)),
                                                  nrow = nrow(x))
      result <- local_specificity(meandf, index, frac)
      if (type == "global") {
          result <- matrix(apply(result, 1, max), nrow = nrow(result), ncol = ncol(result))
      }
      return(as.vector(result[observed.otu]))
  }
  ## One sample specificity (not stratified)
  one_rand_sample_spec <- function(df, index, group, type) {
    return(local_spec(df[ , sample.int(length(group), replace = replace)], index, group, type))
  }

  ## Compute original sample diversity, add local relative abundance for each otu
  meandf <- t(rowsum(t(tdf), group, reorder = TRUE))
  relative.abundance <- apply(meandf, 2, function(x) x / sum(x))
  observed.otu <- as.vector(relative.abundance > 0)
  res <- data.frame(otu = rep(rownames(relative.abundance), ncol(relative.abundance)), 
                    level = rep(colnames(relative.abundance), each = nrow(relative.abundance)),
                    abundance = as.vector(relative.abundance))
  res <- res[observed.otu, ]
  res$specificity <- local_spec(tdf, index, group, type = "local")
  ## Replicate specificity estimates
  cat("Estimating p-values, may take a few minutes", sep = "\n")
  if (parallel) {
      resmat <- mclapply(1:B, function(i) {
          cat(paste("bootstrap sample", i), sep = "\n")
          one_rand_sample_spec(tdf, index, group, type)} )
      resmat <- do.call(cbind, resmat)
      resmat <- cbind(resmat, res$specificity)
  } else {
      resmat <- cbind(replicate(n = B, one_rand_sample_spec(tdf, index, group, type)), res$specificity)
  }
  distrib <- t(apply(resmat, 1, function(x) c(mean = mean(x), sd = sd(x), min = min(x),
                                              quantile(x, probs = c(0.5, 0.75, 0.9, 0.95, 0.99)))))
  colnames(distrib) <- c("bmean", "bsd", "bmin", "bq50", "bq75", "bq90", "bq95", "bq99")
  rawp <- rowMeans(sweep(resmat, 1, res$specificity, ">="))
  res <- data.frame(res, rawp = rawp, adjp = p.adjust(rawp, method = "fdr"), distrib)
  ## Add method
  attr(res, "index") <- index
  attr(res, "group") <- group
  return(res)
}


## Compare saturation speed of rarefaction curves for specific species against all species. 
specific_rarefaction <- function(physeq, step = 10, group, index = c("fraction", "indspec"),
                                 specificity = NULL, threshold = 0.9, pvalue = 0.05, 
                                 label = NULL, B = 999, 
                                 color = NULL, parallel = FALSE, se = TRUE, plot = TRUE) {
    ## Args:
    ## - physeq: phyloseq class object, from which abundance data are extracted
    ## - step: Step size for sample size in rarefaction curves
    ## - group:  Either the a single character string matching a
    ##           variable name in the corresponding sample_data of ‘physeq’, or a
    ##           factor with the same length as the number of samples in ‘physeq’.
    ## - specifity: Default `NULL`. Data frame with components "specificity" (observed local specificity) and 
    ##              "group" (factor level). 
    ## - index: Default "indspec". Used if specifity is `NULL`, argument passed on to
    ##          'link{estimate_local_specificity}' to compute specificity
    ## - B: (Optional) Default 999. Integer value. Number of replicates used to compute specificity adjusted p-values. 
    ## - threshold: Default 0.9. Numeric value used to label an otu as specific when
    ##              constructing rarefaction curve. ‘Specificity’ must be higher
    ##              than 'threshold'. 
    ## - pvalue: Default 0.05. Numeric value used to label an otu as specific when
    ##           constructing rarefaction curve. Used only if threshold is NULL.
    ##           ‘adjp’ must be lower than 'pvalue'.
    ## - label: (Optional) Default `NULL`. Character string. The name of the variable
    ##          to map to text labels on the plot. Similar to color option
    ##          but for plotting text.
    ## - color: (Optional). Default ‘NULL’. Character string. The name of the
    ##          variable to map to colors in the plot. This can be a sample
    ##          variable (among the set returned by
    ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
    ##          returned by ‘rank_names(physeq)’).
    ##          Finally, The color scheme is chosen automatically by
    ##          ‘link{ggplot}’, but it can be modified afterward with an
    ##          additional layer using ‘scale_color_manual’.
    ## - plot:  Logical, should the graphic be plotted.
    ## - parallel: should rarefaction be parallelized (using parallel framework)
    ## - parallel: should richness standard errors be computed
    ## require vegan
    x <- as(otu_table(physeq), "matrix")
    if (taxa_are_rows(physeq)) { x <- t(x) }

    ## Get sample data 
    if (!is.null(sample_data(physeq, FALSE))) {
        sdf <- as(sample_data(physeq), "data.frame")
        sdf$Sample <- rownames(sdf)
    } else {
        stop("Sample data required in phyloseq class argument")
    }
    
    ## Compute specificity, and if necessary adjusted p-values for
    ## local specificity
    if (is.null(specificity)) {
        specificity <- estimate_local_specificity(physeq = physeq, group = group,
                                                  index = index, freq = TRUE,
                                                  B = B, se = FALSE, parallel = parallel)
    }
    if (is.null(threshold)) {
        test_specificity <- test_local_specificity(physeq = physeq,
                                                   group = group,
                                                   index = index, freq = TRUE,
                                                   B = B, parallel = parallel)
        specificity <- merge(specificity, test_specificity)
        attr(specificity, "index") <- attr(test_specificity, "index")
        attr(specificity, "group") <- attr(test_specificity, "group")
    }
    
    ## Call specific otus on either p-value or specificity value
    if (is.null(threshold)) {
        specificity$specific <- with(specificity, adjp < pvalue)
    } else {
        specificity$specific <- with(specificity, specificity > threshold)
    }
    
    ## Local modification of vegan 'rarefy' function
    ## Used to compute estimated relative richness when limiting oneself to a
    ## species subset only. Richness is relative
    .rarefy <- function(x, sample, se = se, subset = rep(TRUE, length(x))) {
        .rarefun <- function(n) {
            J <- sum(x)
            x <- x[subset]
            x <- x[x > 0]
            ldiv <- lchoose(J, n)
            p1 <- ifelse(J - x < n, 0, exp(lchoose(J-x, n) - ldiv))
            out <- sum(1 - p1)
            if (se) {
                V <- sum(p1 * (1 - p1))
                Jxx <- J - outer(x, x, "+")
                ind <- lower.tri(Jxx)
                Jxx <- Jxx[ind]
                V <- V + 2 * sum(ifelse(Jxx < n, 0,
                                        exp(lchoose(Jxx, n) - ldiv)) - outer(p1, p1)[ind])
                out <- cbind(out, sqrt(V))
            }
            out
        }
        ## Normalization by maximum richness for comparable scales
        S.rare <- sapply(sample, function(n) .rarefun(n)) / length(x[subset & (x > 0)])
        if (se) {
            rownames(S.rare) <- c("S", "se")
        } else {
            dim(S.rare) <- c(1, length(S.rare))
            rownames(S.rare) <- c("S")
        }
        return(S.rare)
    }
    
    ## This script is adapted from vegan `rarecurve` function
    ## Computes rarefaction curves (otu richness) and se.     
    tot <- rowSums(x)
    S <- rowSums(x > 0)
    nr <- nrow(x)

    ## For complete richness
    if (parallel) {
        rarefun <- function(i) {
            cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
            n <- seq(1, tot[i], by = step)
            if (n[length(n)] != tot[i]) 
                n <- c(n, tot[i])
            y <- .rarefy(x[i, ], n, se = se)
            return(data.frame(t(y), Size = n, Sample = rownames(x)[i], Richness = "Complete"))
        }
        out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
    } else {
        out <- lapply(seq_len(nr), function(i) {
            cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
        n <- seq(1, tot[i], by = step)
            if (n[length(n)] != tot[i]) 
                n <- c(n, tot[i])
            y <- .rarefy(x[i, ], n, se = se)
            return(data.frame(t(y), Size = n, Sample = rownames(x)[i], Richness = "Complete"))
        })
    }
    df.comp <- do.call(rbind, out)

    ## For specific richness
    if (parallel) {
        rarefun <- function(i) {
            cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
            n <- seq(1, tot[i], by = step)
            if (n[length(n)] != tot[i]) 
                n <- c(n, tot[i])
            i.level <- attr(specificity, "group")[i]
            specific <- subset(specificity,
                               specific & level == i.level,
                               select = "otu")
            specific <- (colnames(x) %in% as.character(specific))
            y <- .rarefy(x[i, ], n, se = se, subset = specific)
            return(data.frame(t(y), Size = n, Sample = rownames(x)[i], Richness = "Specific"))
        }
        out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
    } else {
        out <- lapply(seq_len(nr), function(i) {
            cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
        n <- seq(1, tot[i], by = step)
            if (n[length(n)] != tot[i]) 
                n <- c(n, tot[i])
            i.level <- attr(specificity, "group")[i]
            specific <- subset(specificity,
                               specific & level == i.level,
                               select = "otu")
            specific <- (colnames(x) %in% as.character(specific$otu))
            y <- .rarefy(x[i, ], n, se = se, subset = specific)
            return(data.frame(t(y), Size = n, Sample = rownames(x)[i], Richness = "Specific"))
        })
    }
    df.spec <- do.call(rbind, out)

    ## Merge all richnesses
    df <- rbind(df.comp, df.spec)
    df$Richness <- factor(df$Richness, levels = c("Specific", "Complete"))
    
    ## Merge sample data with rarefaction curve
    sdf$Group <- attr(specificity, "group")
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = 1, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
    
    ## Add, any custom-supplied plot-mapped variables
    if( length(color) > 1 ){
        data$color <- color
        names(data)[names(data)=="color"] <- deparse(substitute(color))
        color <- deparse(substitute(color))
    }
    if( length(label) > 1 ){
        labels$label <- label
        names(labels)[names(labels)=="label"] <- deparse(substitute(label))
        label <- deparse(substitute(label))
    }
    
    p <- ggplot(data = data, aes_string(x = "Size", y = "S",
                    group = "interaction(Richness, Sample)", color = color))
    p <- p + labs(x = "Sample Size", y = "(Relative) Species Richness")
    if (!is.null(label)) {
        p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                           size = 4, hjust = 0)
    }
    p <- p + facet_wrap(~Group)
    p <- p + geom_line(aes(linetype = Richness))
    if (se) {
        p <- p + geom_ribbon(aes_string(ymin = "S - se", ymax = "S + se", color = NULL, fill = color), alpha = 0.1)
    }

    if (plot) {
        plot(p)
    }

    invisible(p)
}
