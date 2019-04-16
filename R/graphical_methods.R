## Set of graphical methods for phyloseq objects (mainly related to ordination,
## and library size after normalisation)
require(ggplot2)
require(scales)
require(reshape2)

rarecurve2 <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", label = TRUE, col = "black", ...)
  ## See documentation for vegan rarecurve, col is now used to define
  ## custom colors for lines and panels
{
  tot <- rowSums(x)
  S <- specnumber(x)
  nr <- nrow(x)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i])
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),
                                           y = z, xout = sample, rule = 1)$y)
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_len(length(out))) {
    color <- col[((ln-1) %% length(col)) + 1]
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = color, ...)
  }
  if (label) {
    ordilabel(cbind(tot, S), labels = rownames(x), col = col, ...)
  }
  invisible(out)
}

## Rarefaction curve, ggplot style
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed.
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }

  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)

  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)

  ## Get sample data
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }

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

  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

## Plot composition at taxaRank2 level within taxa taxaSet1 at taxaRank1 level
## Restricts plot to numberOfTaxa taxa
plot_composition <- function(physeq,
                             taxaRank1 = "Phylum",
                             taxaSet1 = "Proteobacteria",
                             taxaRank2 = "Family",
                             numberOfTaxa = 9, fill = NULL,
                             x = "Sample",
                             y = "Abundance", facet_grid = NULL) {
  ## Args:
  ## - physeq: phyloseq class object
  ## - taxaRank1: taxonomic level in which to do the first subsetting
  ## - taxaSet1: subset of level taxaRank1 to use
  ## - taxaRank2: taxonomic level used to agglomerate
  ## - numberOfTaxa: number of (top) taxa to keep at level taxaRank2
  ##
  ## Returns:
  ## - ggplot2 graphics
  ggdata <- ggformat(physeq, taxaRank1, taxaSet1, taxaRank2,
                     fill, numberOfTaxa)
  p <- ggplot(ggdata, aes_string(x = x, y = y, fill = fill, color = fill, group = "Sample"))
  ## Manually change color scale to assign grey to "Unknown" (if any)
  if (!is.null(fill) && any(c("Unknown", "Other") %in% unique(ggdata[, fill]))) {
      ranks <- as.character(unique(ggdata[, fill]))
      ranks <- ranks[ ! ranks %in% c("Multi-affiliation", "Unknown", "Other")]
      colvals <- c(gg_color_hue(length(ranks)),
                   "grey75", "grey45", "black")
      names(colvals) <- c(ranks, "Multi-affiliation", "Unknown", "Other")
      ## Now add the manually re-scaled layer with Unassigned as grey
      p <- p + scale_fill_manual(values=colvals) + scale_color_manual(values = colvals)

  }
  p <- p + geom_bar(stat = "identity", position = "stack")
  if ( !is.null(facet_grid)) {
    p <- p + facet_grid(facets = facet_grid, scales = "free_x")
  }
  p <- p + theme(axis.text.x=element_text(angle=90), axis.title.x=element_blank())
  p <- p + ggtitle(paste("Composition within", taxaSet1, "(", numberOfTaxa, "top", taxaRank2, ")"))
  return(p)
}

correct_levels <- function(physeq, DF, map.var) {
  oldLevels <- character(0)
  if (any(DF[ , map.var] == "samples", na.rm = TRUE)) {
    oldLevels <- unique(tax_table(physeq)[ , map.var])
  }
  if (any(DF[ , map.var] == "taxa", na.rm = TRUE)) {
    oldLevels <- levels(get_variable(physeq, map.var))
  }
  allLevels <- unique(c(as.character(DF[, map.var]), oldLevels))
  allLevels <- allLevels[!is.na(allLevels)]
  DF[ , map.var] <- factor(DF[ , map.var], levels = allLevels)
  return(DF)
}



## Find numberOfTaxa most abundant taxa at taxaRank level and return their relative
## abundance in all sample
top_taxa_abundance<- function(physeq, numberOfTaxa = 9, raw = FALSE) {
    ## Args:
    ## - physeq: phyloseq class object
    ## - raw: (Required). Defaults to 'FALSE', logical. Should abundances be transformed to
    ##        frequencies before sorting otus from most to less abundant.
    ## - numberOfTaxa: number of (top) taxa to keep
    ##
    ## Returns:
    ## - x: Data frame with numberOfTaxa rows (one per final otu) and column 'Abundance' standing
    ##      for  total counts (if raw = TRUE) or average frequency
    ##      in samples of physeq. Rownames of x are either otus' names (if taxaRank = NULL) or names
    ##      corresponding to taxonomic rank 'TaxaRank'. All NA ranks are assigned to 'Unassigned'.
  stopifnot(!is.null(tax_table(physeq, FALSE)))
  otutab <- otu_table(physeq)
  if ( !taxa_are_rows(otutab) ) {otutab = t(otutab)}
  otutab <- as(otutab, "matrix")
  if (raw) {
      Abundance <- rowSums(otutab)
  } else {
      otutab <- apply(otutab, 2, function(x) x / sum(x))
      Abundance <- rowMeans(otutab)
  }
  ## Get top taxa
  mdf <- data.frame(OTU = names(Abundance), Abundance = Abundance)
  ## Add taxonomic information
  tax <- as(tax_table(physeq), "matrix")
  tax <- data.frame(OTU = rownames(tax), tax)
  mdf <- merge(mdf, tax, by.x = "OTU")
  ## Keep only numberOfTaxa top taxa
  topTaxa <- names(sort(Abundance, decreasing = TRUE))[1:numberOfTaxa]
  mdf <- mdf[ match(topTaxa, mdf$OTU), ]
  mdf$OTU <- factor(mdf$OTU, levels = unique(mdf$OTU))
  return(mdf)
}


## Find numberOfTaxa most abundant taxa at taxaRank level
top_taxa <- function(physeq, taxaRank, numberOfTaxa = 9) {
  stopifnot(!is.null(tax_table(physeq, FALSE)))
  otutab <- otu_table(physeq)
  if ( !taxa_are_rows(otutab) ) {otutab = t(otutab)}
  otutab <- as(otutab, "matrix")
  otutab <- apply(otutab, 2, function(x) x / sum(x))
  ## Subset to OTUs belonging to taxaSet1 to fasten process
  stopifnot(taxaRank %in% colnames(tax_table(physeq)))
  mdf <- melt(data = otutab, varnames = c("OTU", "Sample"))
  colnames(mdf)[3] <- "Abundance"
  mdf <- mdf[mdf$Abundance > 0, ]
  ## Add taxonomic information
  tax <- as(tax_table(physeq), "matrix")
  tax <- data.frame(OTU = rownames(tax), tax)
  mdf <- merge(mdf, tax, by.x = "OTU")
  ## Aggregate by taxaRank and recover most abundant taxa
  abundanceByTaxa <- aggregate(as.formula(paste("Abundance ~", taxaRank)), data = mdf, FUN = sum)
  ii <- order(abundanceByTaxa$Abundance, decreasing = TRUE)
  ## Keep only numberOfTaxa top taxa
  topTaxa <- (abundanceByTaxa[ii, taxaRank])[1:min(numberOfTaxa, nrow(abundanceByTaxa))]
  return(as(topTaxa, "character"))
}

## Find numberOfConditions most abundant conditions in variable
top_conditions <- function(physeq, variable, numberOfConditions = 9) {
  stopifnot(!is.null(sample_data(physeq, FALSE)))
  var <- get_variable(physeq, variable)
  sortedCond <- names(sort(table(var), decreasing = TRUE))
  topCond <- sortedCond[1:min(numberOfConditions, length(sortedCond))]
  return(topCond)
}

## Aggregate coordinates by replicate to have "mean" coordinates
## for a given replicate condition.
replicate_means <- function(DF, replicate) {
  Axis1 <- colnames(DF)[1]
  Axis2 <- colnames(DF)[2]
  res1 <- aggregate(as.formula(paste(Axis1, "~", replicate)), data = DF, FUN = mean)
  res2 <- aggregate(as.formula(paste(Axis2, "~", replicate)), data = DF, FUN = mean)
  res <- merge(res1, res2, by.x = replicate, sort = FALSE)
  return(res)
}


## Custom modification of plot_ordination to label outliers species and add mean location
## of samples aggregated by variable. If shape maps to more than 9 symbols,
## use only symbols corresponding to 9 most abundant symbols
plot_biplot <- function(physeq, ordination, axes=c(1, 2), color = NULL, replicate = color,
                        shape = NULL, label = NULL, title = NULL, outlier = NULL) {
  DF <- plot_ordination(physeq, ordination, "biplot", axes, color, shape, label, title, TRUE)
  ## Retrieve correct levels
  if(!is.null(shape)) { DF <- correct_levels(physeq, DF, shape) }
  if(!is.null(color)) { DF <- correct_levels(physeq, DF, color) }
  ## Keep only top 9 shapes
  if (!is.null(shape) & length(levels(DF[ ,shape])) > (9+1)) {
    if ( !"taxa" %in% levels(DF[, shape])) {
      DF[ , shape] <- factor(DF[ ,shape], levels = c(top_taxa(physeq, shape), "samples"))
    } else {
      DF[ , shape] <- factor(DF[ ,shape], levels = c(top_conditions(physeq, shape), "taxa"))
    }
  }
  ## Name dimensions
  x <- colnames(DF)[1]
  y <- colnames(DF)[2]
  ## label outliers
  if ( !is.null(outlier) ) {
    outliers <- taxa_names(physeq)[taxa_sums(physeq) > sum(otu_table(physeq)) * outlier]
  }
  ## Get mean of each batch of replicate
  if (!is.null(replicate)) { sampleCoordinates <- replicate_means(DF, replicate) }
  ## Mapping section
  p <- ggplot(DF, aes_string(x = x, y = y))
  ## Plot building (do not include aes in main figure, allows displays of further layers)
  if ( is.null(color) ) {
    ## Rename color title in legend
    p <- p + geom_point(aes_string(color = "id.type", shape = shape), na.rm = TRUE)
    p <- update_labels(p, list(colour = "type"))
  } else {
    p <- p + geom_point(aes_string(size = "id.type", color = color, shape = shape), na.rm = TRUE)
    ## Check if variable is discrete
    if( is.discrete(DF[, color]) ){
      colvals <- gg_color_hue(length(levels(as(DF[, color], "factor"))))
      names(colvals) <- levels(as(DF[, color], "factor"))
      ## Now make the taxa or samples dark grey
      colvals[names(colvals) %in% c("samples", "taxa")] <- "grey45"
      ## Now add the manually re-scaled layer with taxa/samples as grey
      p <- p + scale_colour_manual(values=colvals)
    }
    ## Adjust size so that samples are bigger than taxa by default.
    p <- p + scale_size_manual("type", values=c(samples=5, taxa=2))
  }

  ## Adjust shape scale
  shape.names <- levels(DF[, shape])
  shape.scale <- 17:(17 + length(shape.names) - 1)
  names(shape.scale) <- shape.names
  shape.scale["samples"] <- 16
  p <- p + scale_shape_manual(values = shape.scale)

  ## Add the text labels
  if( !is.null(label) ){
    label_map <- aes_string(x=x, y=y, color = color, label=label, na.rm=TRUE)
    p <- p + geom_text(label_map, data=DF[ !is.na(DF[ , label]), ],
                       size=3, vjust=1.5, na.rm=TRUE)
  }

  ## Add replicates
  if( !is.null(replicate) ){
    rep_map <- aes_string(x=x, y=y, label=replicate, na.rm=TRUE)
    if (replicate == color) {
      repcols <- colvals[sampleCoordinates[ , replicate]]
    } else {
      repcols <- "black"
    }
    p <- p + geom_text(rep_map, data = sampleCoordinates,
                       colour = repcols,
                       size=4, vjust=1.5, na.rm=TRUE)
  }

  ## Optionally add a title to the plot
  if( !is.null(title) ){
    p <- p + ggtitle(title)
  }

  ## Add fraction variability to axis labels, if available
  if( length(extract_eigenvalue(ordination)[axes]) > 0 ){
      ## Only attempt to add fraction variability
      ## if extract_eigenvalue returns something
      eigvec = extract_eigenvalue(ordination)
      ## Fraction variability, fracvar
      fracvar = eigvec[axes] / sum(eigvec)
      ## Percent variability, percvar
      percvar = round(100*fracvar, 1)
      ## The string to add to each axis label, strivar
      ## Start with the curent axis labels in the plot
      strivar = as(c(p$label$x, p$label$y), "character")
      ## paste the percent variability string at the end
      strivar = paste0(strivar, "   [", percvar, "%]")
      ## Update the x-label and y-label
      p = p + xlab(strivar[1]) + ylab(strivar[2])
  }

  ## Add outliers attribute
  if ( !is.null(outlier) ) {attr(p, "outliers") <- outliers}
  if ( !is.null(replicate) ) {attr(p, "repcols") <- repcols}

  ## return the ggplot object
  return(p)
}


## Custom modification of plot_ordination to  add mean location
## of samples aggregated by variable.
plot_samples <- function(physeq, ordination, axes=c(1, 2), color = NULL,
                         replicate = color,
                         shape = NULL, label = NULL, title = NULL) {
  DF <- plot_ordination(physeq, ordination, "samples", axes, color, shape, label, title, TRUE)
  ## Retrieve correct levels
  if(!is.null(shape)) { DF <- correct_levels(physeq, DF, shape) }
  if(!is.null(color) & !is.null(replicate)) {
    if (replicate == color) {
      if (is.numeric(DF[, color]))
        message("Using quantitative variable as grouping variable, try setting 'replicate = NULL' for better results.")
      DF <- correct_levels(physeq, DF, color)
    }
  }
  ## Name dimensions
  x <- colnames(DF)[1]
  y <- colnames(DF)[2]
  ## Get mean of each batch of replicate
  if (!is.null(replicate)) { sampleCoordinates <- replicate_means(DF, replicate) }
  ## Mapping section
  p <- ggplot(DF, aes_string(x = x, y = y, color = color, shape = shape))
  ## Plot building
  p <- p + geom_point(na.rm = TRUE)

  ## Add the text labels
  if( !is.null(label) ){
    label_map <- aes_string(x=x, y=y, label=label, color = color)
    p <- p + geom_text(label_map, data=DF[!is.na(DF[ , label]) , ],
                       size=3, vjust=1.5, na.rm=TRUE)
  }

  ## Add replicates
  if( !is.null(replicate) ){
    if (color == replicate) { ## map color to replicate if color and replicate grouping agree
      rep_map <- aes_string(x=x, y=y, label=replicate, color = replicate, shape = NULL)
    } else { ## don't set color aes
      rep_map <- aes_string(x=x, y=y, label=replicate, color = NULL, shape = NULL)
    }
    p <- p + geom_text(rep_map, data = sampleCoordinates,
                       size=4, vjust=1.5)
  }

  ## Optionally add a title to the plot
  if( !is.null(title) ){
    p <- p + ggtitle(title)
  }

  ## Add fraction variability to axis labels, if available
  if( length(extract_eigenvalue(ordination)[axes]) > 0 ){
      ## Only attempt to add fraction variability
      ## if extract_eigenvalue returns something
      eigvec = extract_eigenvalue(ordination)
      ## Fraction variability, fracvar
      fracvar = eigvec[axes] / sum(eigvec)
      ## Percent variability, percvar
      percvar = round(100*fracvar, 1)
      ## The string to add to each axis label, strivar
      ## Start with the curent axis labels in the plot
      strivar = as(c(p$label$x, p$label$y), "character")
      ## paste the percent variability string at the end
      strivar = paste0(strivar, "   [", percvar, "%]")
      ## Update the x-label and y-label
      p = p + xlab(strivar[1]) + ylab(strivar[2])
  }

  ## return the ggplot object
  return(p)
}

## ggplot hue color scale
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}


## Plotting fonction once library sizes have been estimated
ggnorm <- function(physeq, cds, x = "X.SampleID", color = NULL, title = NULL) {
  ## Args:
  ## - cds: Count data set (class eSet of BioConductor)
  ## - x:  ggplot x aesthetics, used for plotting
  ## - color: ggplot color aes, used for plotting
  ## - title: optional, plot title
  ## - physeq: phyloseq class object from which metadata are extracted
  ##
  ## Returns:
  ## - ggplot2 figure with distribution of normalized log2(counts+1) on left, colored by
  ##   color and mean/variance function on right panel.
  countdf <- counts(cds, normalize = TRUE)
  countdf <- log2(countdf+1)
  countdf[countdf == 0] <- NA
  countdf <- melt(countdf, varnames = c("OTU", "X.SampleID"))
  countdf <- countdf[!is.na(countdf$value), ]
  countdf <- merge(countdf, as(sample_data(physeq), "data.frame"))
  p <- ggplot(countdf, aes_string(x = x, y = "value", color = color)) + geom_boxplot()
  p <- p + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete("Sample")
  p <- p + scale_y_discrete("log2(Normalized counts + 1)")
  if (!is.null(title)) {  p <- p + ggtitle(title) }
  ## Open graphical device
  par(no.readonly=TRUE)
  plot.new()
  ## setup layout
  gl <- grid.layout(nrow=1, ncol=2, widths=unit(c(1,1), 'null'))
  ## grid.show.layout(gl)
  ## setup viewports
  vp.1 <- viewport(layout.pos.col=1) # boxplot
  vp.2 <- viewport(layout.pos.col=2) # mean - sd plot
  ## init layout
  pushViewport(viewport(layout=gl))
  ## Access the first viewport
  pushViewport(vp.1)
  ## print our ggplot graphics here
  print(p, newpage=FALSE)
  ## done with the viewport
  popViewport()
  ## Move to the second viewport
  pushViewport(vp.2)
  ##  start new base graphics in second viewport
  par(new=TRUE, fig=gridFIG())
  meanSdPlot(log2(counts(cds, normalized = TRUE)[, ] + 1), ranks = FALSE, main = title)
  ##  done with the viewport
  popViewport()
}

## return data frame for relative abundance plots in ggplot2
## Return relative abundance of top NumberOfTaxa OTUs at the taxaRank2 level
## within taxaSet1 at taxaRank1 level
ggformat <- function(physeq, taxaRank1 = "Phylum", taxaSet1 = "Proteobacteria",
                     taxaRank2 = "Family", fill = NULL, numberOfTaxa = 9) {
    ## Args:
    ## - physeq: phyloseq class object
    ## - taxaRank1: taxonomic level in which to do the first subsetting
    ## - taxaSet1: subset of level taxaRank1 to use
    ## - taxaRank2: taxonomic level used to agglomerate
    ## - numberOfTaxa: number of (top) taxa to keep at level taxaRank2
    ##
    ## Returns:
    ## - data frame with taxonomic, sample and otu counts informations melted together
    ## Enforce orientation and transform count to relative abundances
    stopifnot(!is.null(sample_data(physeq, FALSE)),
              !is.null(tax_table(physeq, FALSE)))
    otutab <- otu_table(physeq)
    if ( !taxa_are_rows(otutab) ) {
    otutab = t(otutab)
    }
    otutab <- as(otutab, "matrix")
    otutab <- apply(otutab, 2, function(x) x / sum(x))
    ## Subset to OTUs belonging to taxaSet1 to fasten process
    if (is.null(fill)) { fill <- taxaRank2 }
    stopifnot(all(c(taxaRank1, taxaRank2, fill) %in% c(rank_names(physeq), "OTU")))
    otutab <- otutab[tax_table(physeq)[ , taxaRank1] %in% taxaSet1, , drop = FALSE]
    if (nrow(otutab) == 0) {
        stop(paste("No otu belongs to", paste(taxaSet1, collapse = ","), "\n",
                   "at taxonomic level", taxaRank1))
    }
    mdf <- melt(data = otutab, varnames = c("OTU", "Sample"))
    colnames(mdf)[3] <- "Abundance"
    ## mdf <- mdf[mdf$Abundance > 0, ] ## Remove absent taxa
    ## Add taxonomic information and replace NA and unclassified Unknown
    tax <- as(tax_table(physeq), "matrix")
    tax[is.na(tax)] <- "Unknown"
    tax[grepl("unknown", tax)] <- "Unknown"
    tax[tax %in% c("", "unclassified", "Unclassified", "NA")] <- "Unknown"
    tax <- data.frame(OTU = rownames(tax), tax)
    mdf <- merge(mdf, tax, by.x = "OTU")
    ## Aggregate by taxaRank2
    mdf <- aggregate(as.formula(paste("Abundance ~ Sample +",
                                      fill, "+", taxaRank2)),
                     data = mdf, FUN = sum)
    topTaxa <- aggregate(as.formula(paste("Abundance ~ ", taxaRank2)), data = mdf, FUN = sum)
    ## Keep only numberOfTaxa top taxa and aggregate the rest as "Other"
    topTax <- as.character(topTaxa[ order(topTaxa[ , "Abundance"], decreasing = TRUE), taxaRank2])
    topTax <- topTax[topTax != "Unknown"]
    topTax <- topTax[1:min(length(topTax), numberOfTaxa)]
    ## Change to character and correct taxonomic levels
    correct_taxonomy <- function(x) {
      c(sort(x[!x %in% c("Multi-affiliation", "Unknown", "Other")]),
        c("Multi-affiliation", "Unknown", "Other"))
    }
    mdf[ , taxaRank2] <- as.character(mdf[ , taxaRank2])
    mdf[ , fill] <- as.character(mdf[ , fill])
    ii <- (mdf[ , taxaRank2] %in% c(topTax, "Unknown"))
    mdf[!ii , taxaRank2] <- "Other"
    mdf[!ii , fill] <- "Other"
    mdf <- aggregate(as.formula(paste("Abundance ~ Sample +", fill,
                                      "+", taxaRank2)), data = mdf, FUN = sum)
    mdf[, taxaRank2] <- factor(mdf[, taxaRank2],
                               levels = correct_taxonomy(topTax))
    mdf[ , fill] <- as.character(mdf[ , fill])
    mdf[, fill] <- factor(mdf[, fill],
                          levels = correct_taxonomy(unique(mdf[, fill])))
    ## Add sample data.frame
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- sample_names(physeq)
    mdf <- merge(mdf, sdf, by.x = "Sample")
    ## Sort the entries by abundance to produce nice stacked bars in ggplot
    mdf <- mdf[ order(mdf[ , fill], mdf$Abundance, decreasing = TRUE), ]
    return(mdf)
}

## Plot a distance matrix as a heatmap with samples sorted according to
## order vector
plot_dist_as_heatmap <- function(dist, order = NULL, title = NULL,
                                 low = "#B1F756", high = "#132B13",
                                 show.names = FALSE) {
  ## Args:
  ## - dist: distance matrix (dist class)
  ## - order: (optional) ordering of the samples of dist for representation
  ## - title: (optional) graph title
  ## - low, high: (optional) Colours for low and high ends of the gradient
  ## - show.names: (optional) Logical. Should sample names be displayed in the heatmap.
  ##               Defaults to FALSE
  ##
  ## Returns:
  ## - a ggplot2 object
  data <- melt(as(dist, "matrix"))
  colnames(data) <- c("x", "y", "distance")
  if (!is.null(order)) {
    data$x <- factor(data$x, levels = order)
    data$y <- factor(data$y, levels = order)
  }
  p <- ggplot(data, aes(x = x, y = y, fill = distance)) + geom_tile()
  p <- p + theme(axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.x = element_text(angle = 90))
  if (!show.names) {
    p <- p + theme(axis.title.x = element_blank(),
                   axis.title.y = element_blank())
  }
  p <- p + scale_fill_gradient(low = low, high = high)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

## Wrapper around hclust to represent clustering tree
## with leaves colored according to some variables
plot_clust <- function(physeq, dist, method = "ward.D2", color = NULL,
                       label = NULL,
                       title = paste(method, "linkage clustering tree"),
                       legend = TRUE,
                       palette = NULL) {
  ## Args:
  ## - physeq: phyloseq class object
  ## - dist: distance matrix (dist class) or character to be used in phyloseq::distance function
  ## - method: (character) linkage method used in hclust, defaults to "ward.D2"
  ## - color: (character) variable name used to color tree leaves. Defaults to NULL
  ## - label: (character) one the sample_variable from physeq
  ## - title: (character) optional. Plot title, defaults to "method" clustering tree.
  ## - palette: (named color vector) optional. Manual color palette
  ## - legend: (logical) add legend or not
  ## Returns:
  ## - a plot object
  if (is.character(color)) {
    legend.title <- NULL
    color <- phyloseq::get_variable(physeq, color)
  } else {
    legend.title <- NULL
    color <- rep("black", nsamples(physeq))
  }
  color <- as.factor(color)
  ## compute distance
  if (is.character(dist)) {
   dist <- dist[1]
   dist <- phyloseq::distance(physeq, method = dist)
  }
  ## automatic color palette: one color per different sample type
  if (is.null(palette)) {
    palette <- hue_pal()(length(levels(color)))
  } else {
    palette <- palette[levels(color)]
  }
  tipColor = col_factor(palette, levels = levels(color))(color)
  ## Change hclust object to phylo object and plot
  clust <- as.phylo(hclust(dist, method = method))
  ## change tip label if needed
  if (!is.null(label)) {
    tip.dict <- setNames(as.character(phyloseq::get_variable(physeq, label)),
                         sample_names(physeq))
    clust$tip.label <- tip.dict[clust$tip.label]
  }
  ## plot clustering tree
  plot(clust,
       tip.color = tipColor,
       direction = "downwards",
       main = title)
  ## add legend (at figure bottom, over 4 columns)
  if (legend) {
    legend("bottom", legend = levels(color) , xpd=NA,
           fill = palette, border = palette,cex=0.8, bty="n",
           ncol=4,  inset = c(0,-0.05))
  }
}

## Extract legend from a ggplot object
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


################################################################################
# Define S3 generic extract_eigenvalue function
# Function is used by `plot_scree` to get the eigenvalue vector from different
# types of ordination objects.
# Used S3 generic in this case because many ordination objects, the input, are
# not formally-defined S4 classes, but vaguely-/un-defined S3.
#' @keywords internal
extract_eigenvalue = function(ordination) UseMethod("extract_eigenvalue", ordination)
# Default is to return NULL (e.g. for NMDS, or non-supported ordinations/classes).
extract_eigenvalue.default = function(ordination) NULL
# for pcoa objects
extract_eigenvalue.pcoa = function(ordination) ordination$values$Relative_eig
# for CCA objects
extract_eigenvalue.cca = function(ordination) c(ordination$CCA$eig, ordination$CA$eig)
# for RDA objects
extract_eigenvalue.rda = function(ordination) c(ordination$CCA$eig, ordination$CA$eig)
# for dpcoa objects
extract_eigenvalue.dpcoa = function(ordination) ordination$eig
# for decorana (dca) objects
extract_eigenvalue.decorana = function(ordination) ordination$evals
# for pca (edgePCA) objects
extract_eigenvalue.pca = function(ordination) ordination$values$Eigenvalues
