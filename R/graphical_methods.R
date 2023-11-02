# Set of graphical methods for phyloseq objects (mainly related to ordination, and library size after normalisation)


#' Rarefaction curves, ggplot styles
#'
#' @param physeq phyloseq class object, from which abundance data are extracted
#' @param step Step size for sample size in rarefaction curves
#' @param label Default `NULL`. Character string. The name of the variable to map to text labels on the plot. Similar to color option but for plotting text.
#' @param color (Optional). Default `NULL`. Character string. The name of the variable to map to colors in the plot. This can be a sample variable (among the set returned by \code{\link[=sample_variables]{sample_variables(physeq)}} or taxonomic rank (among the set returned by \code{[=rank_names]{rank_names(physeq)}}). Finally, The color scheme is chosen automatically by \link{ggplot}, but it can be modified afterward with an additional layer using \link{scale_color_manual}.
#' @param plot Logical, should the graphic be plotted.
#' @param parallel should rarefaction be parallelized (using parallel framework)
#' @param se Default TRUE. Logical. Should standard errors be computed.
#' @param verbose (Optional). Verbose (TRUE, default) or silent (FALSE) output
#'
#' @import ggplot2 scales dplyr
#' @export
#'
#' @examples
#' data(food)
#' ggrare(food, step = 100, color = "EnvType", se = FALSE)
#' @importFrom dplyr bind_rows inner_join mutate
#' @importFrom ggplot2 aes geom_line geom_ribbon geom_text ggplot labs
#' @importFrom methods as
#' @importFrom parallel mclapply
#' @importFrom phyloseq otu_table sample_data sample_names taxa_are_rows
#' @importFrom rlang sym
#' @importFrom tibble tibble
#' @importFrom vegan rarefy
ggrare <- function(physeq, step = 10, label = NULL, color = NULL,
                   plot = TRUE, parallel = FALSE, se = TRUE, verbose = TRUE) {
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) {
    x <- t(x)
  }

  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)

  rarefun <- function(i) {
    if (verbose) cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, , drop = FALSE], n, se = se)
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
  df <- dplyr::bind_rows(out)

  ## Get sample data
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame") %>%
      dplyr::mutate(Sample = sample_names(physeq))
    data <- dplyr::inner_join(df, sdf, by = "Sample")
    labels <- dplyr::tibble(x = tot, y = S, Sample = sample_names(physeq)) %>%
      inner_join(sdf, by = "Sample")
  }

  p <- ggplot(data = data, aes(x = Size, y = .S, group = Sample))
  if (!is.null(color)) p <- p + aes(color = !!sym(color), fill = !!sym(color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(
      data = labels, aes(x = x, y = y, label = !!sym(label)),
      size = 4, hjust = 0
    )
  }
  p <- p + geom_line()

  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes(ymin = .S - .se, ymax = .S + .se, color = NULL), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

#' Versatile function for plotting compositions graphs
#'
#' @param physeq phyloseq class object
#' @param taxaRank1 Taxonomic rank used for subsetting. If NULL, \code{taxaSet1} is ignored and no subsetting of taxa is performed.
#' @param taxaSet1 Subset of taxa (at rank \code{taxaRank1}) to keep for the plots. If NULL, no subsetting is applied.
#' @param taxaRank2 taxonomic level used to agglomerate, use "OTU", "ASV" or NULL to enforce no aggregation
#' @param numberOfTaxa number of (most abundant) taxa to keep at level taxaRank2
#' @param startFrom \code{startFrom - 1} is the number of (most abundant) taxa to discard before
#'                  selecting the `numberOfTaxa` most abundant taxa.
#' @param fill Taxonomic rank used for filling (should be between taxaRank1 and taxaRank2 in the hierarchy of ranks)
#' @param x Variable mapped to x-axis
#' @param y Variable mapped to y-axis
#' @param facet_grid variable used for faceting.
#' @param sampleOrder Sample order (ignored if \code{x} is different from "Samples")
#' @param taxaOrder Taxa order, either "abundance" (default) or "name" (alphabetical order)
#' @param spread If `TRUE` spread taxonomy from top to bottom using [tax_spread()] to avoid Unknown, Multi-affiliation and NA from showing up in the plot.
#' @param ... Additional arguments passed on to geom_bar
#'
#' @details Allows the user to restrict the plot to taxa in the set `taxaSet1` at level `taxaRank1` and aggregate taxa at level `taxaRank2` in the plot. The plot is limited to `numberOfTaxa` taxa (defaults to 9) and other taxa are automatically lumped in the "Other" category.
#'
#' @return a ggplot2 graphic
#' @export
#'
#' @examples
#' data(food)
#' plot_composition(food, "Kingdom", "Bacteria", "Phylum", facet_grid = "~EnvType")
#' plot_composition(food, taxaRank1 = "Family", taxaSet1 = "Flavobacteriaceae", taxaRank2 = "Species", facet_grid = "~EnvType", spread = TRUE)
#' ## Contrast with
#' plot_composition(food, taxaRank1 = "Family", taxaSet1 = "Flavobacteriaceae", taxaRank2 = "Species", facet_grid = "~EnvType")
#' plot_composition(food, taxaRank1 = "Family", taxaSet1 = "Flavobacteriaceae", taxaRank2 = "Species", fill = "Phylum", facet_grid = "~EnvType")
#' ## Change taxa ordering
#' plot_composition(food, "Kingdom", "Bacteria", "Phylum", taxaOrder = "name", facet_grid = "~EnvType")
#' @importFrom ggplot2 aes aes_string element_blank element_text expansion facet_grid geom_bar ggplot ggtitle labs scale_color_manual scale_fill_manual scale_y_continuous theme theme_bw
#' @importFrom phyloseq rank_names
#' @importFrom scales brewer_pal hue_pal
#' @importFrom stats setNames
plot_composition <- function(physeq,
                             taxaRank1 = NULL,
                             taxaSet1 = NULL,
                             taxaRank2 = "Phylum",
                             numberOfTaxa = 9,
                             fill = NULL,
                             startFrom = 1,
                             x = "Sample",
                             y = "Abundance",
                             sampleOrder = NULL,
                             taxaOrder = "abundance",
                             facet_grid = NULL,
                             spread = FALSE,
                             ...) {
  if (is.null(taxaRank2) || taxaRank2 %in% c("OTU", "ASV")) {
    taxaRank2 <- "OTU_rank"
  }
  if (is.null(taxaRank1) || is.null(taxaSet1)) {
    taxaRank1 <- rank_names(physeq)[1]
    taxaSet1 <- NULL
  }
  if (is.null(fill)) {
    fill <- taxaRank2
  }
  if (fill %in% c("OTU", "ASV")) fill <- "OTU_rank"

  ggdata <- ggformat(
    physeq, taxaRank1, taxaSet1, taxaRank2,
    fill, numberOfTaxa, startFrom, taxaOrder, spread
  ) %>% droplevels()
  if (!is.null(sampleOrder)) ggdata$Sample <- factor(ggdata$Sample, levels = sampleOrder)

  p <- ggplot(
    ggdata,
    aes_string(x = x, y = y, fill = fill, group = "Sample")
  )
  ## Manually change color scale to assign grey to "Unknown" (if any)
  if (!is.null(fill) && any(c("Unknown", "Other") %in% levels(ggdata[[fill]]))) {
    levels <- levels(ggdata[[fill]])
    proper_levels <- setdiff(levels, c("Multi-affiliation", "Unknown", "Other"))

    if (length(proper_levels) > 12) {
      warning("Too many taxa: reverting from Brewer 'Paired' scale to a hue scale.\nThe human eye has trouble picking more than 12 colors, consider showing a smaller number of taxa.")
      fill_cols <- scales::hue_pal()(length(proper_levels))
    } else {
      fill_cols <- scales::brewer_pal(palette = "Paired")(length(proper_levels))
    }

    colvals <- c(
      setNames(fill_cols, proper_levels),
      "Multi-affiliation" = "grey75",
      "Unknown"           = "grey45",
      "Other"             = "black"
    )
    colvals <- colvals[levels]
    ## Now add the manually re-scaled layer with Unassigned as grey
    p <- p +
      scale_fill_manual(values = colvals)

    if (fill != taxaRank2) {
      colvals[] <- "grey20"
      p <- p +
        aes(color = .data[[fill]]) +
        scale_color_manual(values = colvals) +
        labs(subtitle = paste("Grey boxes correspond to", taxaRank2, "within the same", fill))
    }
  }

  p <- p + geom_bar(stat = "identity", position = "stack", ...)
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid, scales = "free_x")
  }
  p <- p +
    scale_y_continuous(expand = expansion(0, 0)) +
    ggtitle(paste0(
      "Composition",
      ## Filter or not
      ifelse(!is.null(taxaSet1),
        paste0(" within ", paste0(taxaSet1, collapse = ", ")),
        ""
      ),
      ## number of taxa represented on the plot
      " (", taxaRank2, " ", startFrom, " to ",
      startFrom + numberOfTaxa - 1, ")"
    )) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90),
      axis.title.x = element_blank()
    ) +
    NULL
  p
}

#' Internal function to correct levels in plot_samples
#' @keywords internal
#' @importFrom phyloseq get_variable tax_table
correct_levels <- function(physeq, DF, map.var) {
  oldLevels <- character(0)
  if (any(DF[, map.var] == "samples", na.rm = TRUE)) {
    oldLevels <- unique(tax_table(physeq)[, map.var])
  }
  if (any(DF[, map.var] == "taxa", na.rm = TRUE)) {
    oldLevels <- levels(get_variable(physeq, map.var))
  }
  allLevels <- unique(c(as.character(DF[, map.var]), oldLevels))
  allLevels <- allLevels[!is.na(allLevels)]
  DF[, map.var] <- factor(DF[, map.var], levels = allLevels)
  return(DF)
}


#' @importFrom methods as
#' @importFrom phyloseq otu_table tax_table taxa_are_rows
top_taxa_abundance <- function(physeq, numberOfTaxa = 9, raw = FALSE) {
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
  if (!taxa_are_rows(otutab)) {
    otutab <- t(otutab)
  }
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
  mdf <- mdf[match(topTaxa, mdf$OTU), ]
  mdf$OTU <- factor(mdf$OTU, levels = unique(mdf$OTU))
  return(mdf)
}


#' Find the most abundant taxa at a given taxonomic rank
#'
#' @param physeq phyloseq class object
#' @param taxaRank taxonomic level used for agglomeration. Default NULL, equivalent to \code{taxaRank = "OTU"} and no agglomeration.
#' @param numberOfTaxa Number of top taxa to return
#'
#' @return A matrix corresponding to the `numberOfTaxa` most abundant taxa at level `taxaRank`-level
#' @export
#'
#' @examples
#' data(food)
#' top_taxa(food, "Phylum", 10)
#' @importFrom methods as
#' @importFrom phyloseq ntaxa rank_names tax_table taxa_names taxa_sums transform_sample_counts
#' @importFrom utils head
top_taxa <- function(physeq, taxaRank = NULL, numberOfTaxa = 9) {
  if (is.null(taxaRank)) {
    taxaRank <- "OTU"
  }
  tax_table(physeq) <- cbind(as(tax_table(physeq), "matrix"),
    OTU = taxa_names(physeq)
  )
  ## Normalize counts and aggregate at level taxaRank
  physeq <- physeq %>%
    transform_sample_counts(function(x) {
      x / sum(x)
    })
  ## Manual tax glom
  if (taxaRank != "OTU") {
    physeq <- fast_tax_glom(physeq, taxrank = taxaRank)
  }
  ## Most abundant taxa
  top_taxa <- taxa_sums(physeq) %>%
    sort(decreasing = TRUE) %>%
    head(n = min(numberOfTaxa, ntaxa(physeq))) %>%
    names()
  ## Corresponding names
  tax_table(physeq)[top_taxa, 1:match(taxaRank, rank_names(physeq))] %>% as("matrix")
}


## Find numberOfConditions most abundant conditions in variable
#' @importFrom phyloseq get_variable sample_data
top_conditions <- function(physeq, variable, numberOfConditions = 9) {
  stopifnot(!is.null(sample_data(physeq, FALSE)))
  var <- get_variable(physeq, variable)
  sortedCond <- names(sort(table(var), decreasing = TRUE))
  topCond <- sortedCond[1:min(numberOfConditions, length(sortedCond))]
  return(topCond)
}

## Aggregate coordinates by replicate to have "mean" coordinates
## for a given replicate condition.
#' @importFrom stats aggregate as.formula
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
#' @importFrom ggplot2 aes_string geom_point geom_text ggplot ggtitle scale_colour_manual scale_shape_manual scale_size_manual update_labels xlab ylab
#' @importFrom methods as
#' @importFrom phyloseq otu_table plot_ordination taxa_names taxa_sums
#' @importFrom plyr is.discrete
plot_biplot <- function(physeq, ordination, axes = c(1, 2), color = NULL, replicate = color,
                        shape = NULL, label = NULL, title = NULL, outlier = NULL) {
  DF <- plot_ordination(physeq, ordination, "biplot", axes, color, shape, label, title, TRUE)
  ## Retrieve correct levels
  if (!is.null(shape)) {
    DF <- correct_levels(physeq, DF, shape)
  }
  if (!is.null(color)) {
    DF <- correct_levels(physeq, DF, color)
  }
  ## Keep only top 9 shapes
  if (!is.null(shape) & length(levels(DF[, shape])) > (9 + 1)) {
    if (!"taxa" %in% levels(DF[, shape])) {
      DF[, shape] <- factor(DF[, shape], levels = c(top_taxa(physeq, shape), "samples"))
    } else {
      DF[, shape] <- factor(DF[, shape], levels = c(top_conditions(physeq, shape), "taxa"))
    }
  }
  ## Name dimensions
  x <- colnames(DF)[1]
  y <- colnames(DF)[2]
  ## label outliers
  if (!is.null(outlier)) {
    outliers <- taxa_names(physeq)[taxa_sums(physeq) > sum(otu_table(physeq)) * outlier]
  }
  ## Get mean of each batch of replicate
  if (!is.null(replicate)) {
    sampleCoordinates <- replicate_means(DF, replicate)
  }
  ## Mapping section
  p <- ggplot(DF, aes_string(x = x, y = y))
  ## Plot building (do not include aes in main figure, allows displays of further layers)
  if (is.null(color)) {
    ## Rename color title in legend
    p <- p + geom_point(aes_string(color = "id.type", shape = shape), na.rm = TRUE)
    p <- update_labels(p, list(colour = "type"))
  } else {
    p <- p + geom_point(aes_string(size = "id.type", color = color, shape = shape), na.rm = TRUE)
    ## Check if variable is discrete
    if (is.discrete(DF[, color])) {
      colvals <- gg_color_hue(length(levels(as(DF[, color], "factor"))))
      names(colvals) <- levels(as(DF[, color], "factor"))
      ## Now make the taxa or samples dark grey
      colvals[names(colvals) %in% c("samples", "taxa")] <- "grey45"
      ## Now add the manually re-scaled layer with taxa/samples as grey
      p <- p + scale_colour_manual(values = colvals)
    }
    ## Adjust size so that samples are bigger than taxa by default.
    p <- p + scale_size_manual("type", values = c(samples = 5, taxa = 2))
  }

  ## Adjust shape scale
  shape.names <- levels(DF[, shape])
  shape.scale <- 17:(17 + length(shape.names) - 1)
  names(shape.scale) <- shape.names
  shape.scale["samples"] <- 16
  p <- p + scale_shape_manual(values = shape.scale)

  ## Add the text labels
  if (!is.null(label)) {
    label_map <- aes_string(x = x, y = y, color = color, label = label, na.rm = TRUE)
    p <- p + geom_text(label_map,
      data = DF[!is.na(DF[, label]), ],
      size = 3, vjust = 1.5, na.rm = TRUE
    )
  }

  ## Add replicates
  if (!is.null(replicate)) {
    rep_map <- aes_string(x = x, y = y, label = replicate, na.rm = TRUE)
    if (replicate == color) {
      repcols <- colvals[sampleCoordinates[, replicate]]
    } else {
      repcols <- "black"
    }
    p <- p + geom_text(rep_map,
      data = sampleCoordinates,
      colour = repcols,
      size = 4, vjust = 1.5, na.rm = TRUE
    )
  }

  ## Optionally add a title to the plot
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  ## Add fraction variability to axis labels, if available
  if (length(extract_eigenvalue(ordination)[axes]) > 0) {
    ## Only attempt to add fraction variability
    ## if extract_eigenvalue returns something
    eigvec <- extract_eigenvalue(ordination)
    ## Fraction variability, fracvar
    fracvar <- eigvec[axes] / sum(eigvec)
    ## Percent variability, percvar
    percvar <- round(100 * fracvar, 1)
    ## The string to add to each axis label, strivar
    ## Start with the curent axis labels in the plot
    strivar <- as(c(p$label$x, p$label$y), "character")
    ## paste the percent variability string at the end
    strivar <- paste0(strivar, "   [", percvar, "%]")
    ## Update the x-label and y-label
    p <- p + xlab(strivar[1]) + ylab(strivar[2])
  }

  ## Add outliers attribute
  if (!is.null(outlier)) {
    attr(p, "outliers") <- outliers
  }
  if (!is.null(replicate)) {
    attr(p, "repcols") <- repcols
  }

  ## return the ggplot object
  return(p)
}


## Custom modification of plot_ordination to  add mean location
## of samples aggregated by variable.
#' @importFrom ggplot2 aes_string geom_point geom_text ggplot ggtitle xlab ylab
#' @importFrom methods as
#' @importFrom phyloseq plot_ordination
plot_samples <- function(physeq, ordination, axes = c(1, 2), color = NULL,
                         replicate = color,
                         shape = NULL, label = NULL, title = NULL) {
  DF <- plot_ordination(physeq, ordination, "samples", axes, color, shape, label, title, TRUE)
  ## Retrieve correct levels
  if (!is.null(shape)) {
    DF <- correct_levels(physeq, DF, shape)
  }
  if (!is.null(color) & !is.null(replicate)) {
    if (replicate == color) {
      if (is.numeric(DF[, color])) {
        message("Using quantitative variable as grouping variable, try setting 'replicate = NULL' for better results.")
      }
      DF <- correct_levels(physeq, DF, color)
    }
  }
  ## Name dimensions
  x <- colnames(DF)[1]
  y <- colnames(DF)[2]
  ## Get mean of each batch of replicate
  if (!is.null(replicate)) {
    sampleCoordinates <- replicate_means(DF, replicate)
  }
  ## Mapping section
  p <- ggplot(DF, aes_string(x = x, y = y, color = color, shape = shape))
  ## Plot building
  p <- p + geom_point(na.rm = TRUE)

  ## Add the text labels
  if (!is.null(label)) {
    label_map <- aes_string(x = x, y = y, label = label, color = color)
    p <- p + geom_text(label_map,
      data = DF[!is.na(DF[, label]), ],
      size = 3, vjust = 1.5, na.rm = TRUE
    )
  }

  ## Add replicates
  if (!is.null(replicate)) {
    if (color == replicate) { ## map color to replicate if color and replicate grouping agree
      rep_map <- aes_string(x = x, y = y, label = replicate, color = replicate, shape = NULL)
    } else { ## don't set color aes
      rep_map <- aes_string(x = x, y = y, label = replicate, color = NULL, shape = NULL)
    }
    p <- p + geom_text(rep_map,
      data = sampleCoordinates,
      size = 4, vjust = 1.5
    )
  }

  ## Optionally add a title to the plot
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  ## Add fraction variability to axis labels, if available
  if (length(extract_eigenvalue(ordination)[axes]) > 0) {
    ## Only attempt to add fraction variability
    ## if extract_eigenvalue returns something
    eigvec <- extract_eigenvalue(ordination)
    ## Fraction variability, fracvar
    fracvar <- eigvec[axes] / sum(eigvec)
    ## Percent variability, percvar
    percvar <- round(100 * fracvar, 1)
    ## The string to add to each axis label, strivar
    ## Start with the curent axis labels in the plot
    strivar <- as(c(p$label$x, p$label$y), "character")
    ## paste the percent variability string at the end
    strivar <- paste0(strivar, "   [", percvar, "%]")
    ## Update the x-label and y-label
    p <- p + xlab(strivar[1]) + ylab(strivar[2])
  }

  ## return the ggplot object
  return(p)
}

## ggplot hue color scale
# gg_color_hue <- function(n) {
#   if (n == 0) return(NULL)
#   scales::hue_pal()(n)
# }


#' Format phyloseq data for easy composition plots
#'
#' @inheritParams plot_composition
#'
#' @return A ggplot friendly data frame with the relative abundance of the top `NumberOfTaxa` OTUs (starting from the `startFrom`th most abundant) at the taxaRank2 level within taxaSet1 at taxaRank1 level
#' @export
#'
#' @examples
#' data(food)
#' ggformat(food)
#' ggformat(food, taxaSet1 = NULL, taxaRank2 = "Phylum")
#' @importFrom dplyr across all_of arrange as_tibble bind_cols case_when desc filter group_by if_else mutate pull row_number select ungroup
#' @importFrom methods as
#' @importFrom phyloseq ntaxa prune_taxa psmelt rank_names sample_data sample_sums tax_table taxa_names taxa_sums transform_sample_counts
#' @importFrom tibble column_to_rownames remove_rownames
ggformat <- function(physeq, taxaRank1 = "Phylum", taxaSet1 = "Proteobacteria",
                     taxaRank2 = "Family", fill = NULL,
                     numberOfTaxa = 9, startFrom = 1,
                     taxaOrder = c("abundance", "name"),
                     spread = FALSE) {
  ## Enforce orientation and transform count to relative abundances
  stopifnot(
    !is.null(sample_data(physeq, FALSE)),
    !is.null(tax_table(physeq, FALSE))
  )

  count_to_prop <- function(x) {
    x / sum(x)
  }
  ## Transform to proportions (if not already)
  if (any(sample_sums(physeq) > 1)) {
    physeq <- transform_sample_counts(physeq, count_to_prop)
  }

  ## Check that taxaranks and fill are propers ranks
  if (is.null(fill)) {
    fill <- taxaRank2
  }
  if (is.null(taxaRank1)) {
    taxaRank1 <- rank_names(physeq)[1]
    taxaSet1 <- NULL
  }
  stopifnot(all(c(taxaRank1, taxaRank2, fill) %in% c(rank_names(physeq), "OTU_rank")))
  ranks <- find_upper_ranks(physeq, c(taxaRank1, taxaRank2, fill))

  ## Subset at TaxaRank1
  if (!is.null(taxaSet1)) {
    physeq <- prune_taxa(tax_table(physeq)[, taxaRank1] %in% taxaSet1, physeq)
    if (ntaxa(physeq) == 0) {
      stop(paste(
        "No otu belongs to", paste(taxaSet1, collapse = ","), "\n",
        "at taxonomic level", taxaRank1
      ))
    }
  }

  ## Spread affiliations
  if (spread) {
    physeq <- tax_spread(physeq)
  }

  ## Correct taxonomy
  tax <- as(tax_table(physeq), "matrix")
  tax <- cbind(tax, OTU_rank = taxa_names(physeq))
  if (!spread) {
    tax[is.na(tax)] <- "Unknown"
    tax[grepl("unknown", tax)] <- "Unknown"
    tax[tax %in% c("", "unclassified", "Unclassified", "NA")] <- "Unknown"
  }
  tax_table(physeq) <- tax

  ## agglomerate at TaxaRank2
  if (taxaRank2 != "OTU_rank") {
    physeq <- fast_tax_glom(physeq, taxrank = taxaRank2)
  } else {
    ranks <- c(ranks, "OTU_rank")
  }

  ## Sort taxa by abundance, remove unwanted taxa and remove affiliation of least abundant taxa
  last_levels <- c("Multi-affiliation", "Unknown", "Other")
  topTaxa <- data.frame(
    abundance = taxa_sums(physeq),
    taxa = taxa_names(physeq),
    stringsAsFactors = FALSE
  ) %>%
    bind_cols(as(tax_table(physeq), "matrix")[, ranks, drop = F] %>% as_tibble()) %>%
    arrange(desc(abundance)) %>%
    mutate(
      rank = row_number(),
      status = case_when(
        row_number() < startFrom - 1 ~ "filtered out",
        row_number() < startFrom + numberOfTaxa ~ "conserved",
        TRUE ~ "aggregated"
      )
    ) %>%
    ## remove unwanted taxa
    filter(status != "filtered out") %>%
    ## propagate Other / Multi-affiliation / Unknown across all ranks in least abundant taxa
    mutate(across(all_of(ranks), ~ case_when(
      status == "conserved" ~ .x,
      status == "aggregated" & .data[[taxaRank2]] %in% last_levels ~ .data[[taxaRank2]],
      TRUE ~ "Other"
    )))


  ## Warning if no taxa remains
  if (nrow(topTaxa) == 0) {
    stop(paste0("Not enough taxa to show. Consider decreasing `startFrom` to a lower value."))
  }

  ### Check that final names of most abundant taxa are unique; if not print a warning message and sanitize them.
  res <- topTaxa %>%
    filter(status == "conserved") %>%
    pull(taxaRank2)
  if (length(res) < numberOfTaxa) {
    warning(paste0(
      "Not enough taxa to show and/or all remaining taxa have unknown affiliation at rank ",
      taxaRank2, ". Consider using a smaller value."
    ))
  }
  problematic_taxa <- duplicated(res) | res %in% last_levels
  if (any(problematic_taxa)) {
    warning(paste("Some of the most abundant taxa are unknown or have the same name at rank", taxaRank2, "but not at upper ranks.\nUsing suffix '_x' to distinguish them. See table for further details."))
    cat("Problematic taxa", sep = "\n")
    topTaxa %>%
      filter(status == "conserved", .data[[taxaRank2]] %in% res[problematic_taxa]) %>%
      select(-abundance, -status) %>%
      print()
    .f <- function(status, x) {
      if (length(x) == 1) {
        return(x)
      }
      if_else(status == "conserved", paste(x, seq_along(x), sep = "_"), x)
    }
    topTaxa <- topTaxa %>%
      group_by(.data[[taxaRank2]]) %>%
      mutate({{ taxaRank2 }} := .f(status, .data[[taxaRank2]])) %>%
      ungroup()
  }

  ## Change taxaRank2/fill to ordered levels and sort topTaxa
  taxaOrder <- match.arg(taxaOrder)
  taxa_levels <- topTaxa %>%
    filter(status == "conserved") %>%
    pull(taxaRank2)
  if (taxaOrder == "name") {
    taxa_levels <- sort(taxa_levels)
  }
  taxa_levels <- c(taxa_levels, last_levels) %>% unique()

  topTaxa <- topTaxa %>%
    mutate({{ taxaRank2 }} := factor(.data[[taxaRank2]], levels = taxa_levels)) %>%
    arrange(.data[[taxaRank2]])

  if (taxaRank2 != fill) {
    fill_levels <- topTaxa %>%
      filter(status == "conserved") %>%
      pull(fill) %>%
      unique()
    if (taxaOrder == "name") {
      fill_levels <- sort(fill_levels)
    }
    fill_levels <- c(fill_levels, last_levels)
    topTaxa <- topTaxa %>%
      mutate({{ fill }} := factor(.data[[fill]], levels = fill_levels))
  } else {
    fill_levels <- taxa_levels
  }

  ## Compact and simplify physeq object
  physeq <- prune_taxa(topTaxa$taxa, physeq)
  tax_table(physeq) <- topTaxa %>%
    tibble::remove_rownames() %>%
    select(all_of(c("taxa", ranks))) %>%
    column_to_rownames(var = "taxa") %>%
    as.matrix()
  physeq <- fast_tax_glom(physeq, taxrank = taxaRank2)

  tdf <- psmelt(physeq) %>%
    mutate(
      {{ taxaRank2 }} := factor(.data[[taxaRank2]], levels = taxa_levels),
      {{ fill }} := factor(.data[[fill]], levels = fill_levels)
    ) %>%
    arrange(desc(.data[[fill]]), Abundance)

  return(tdf)
  # ## Bind with undefined taxa and propagate Unknown / Multi-affiliation / Other across all ranks
  # topTaxa <- bind_rows(topTaxa, topTaxa_undefined) %>%
  #   mutate(across(all_of(ranks), ~ if_else(.data[[taxaRank2]] %in% last_levels, .data[[taxaRank2]], .x)))



  #
  #
  # ## TOD0: change taxaRank2 and fill directly to factors in topTaxa.
  # ## Change to character and correct taxonomic levels
  # correct_taxonomy <- function(x) {
  #   last_levels <- c("Multi-affiliation", "Unknown", "Other")
  #   res <- c(setdiff(x, last_levels), last_levels)
  #   if (any(duplicated(res))) {
  #     warning(paste("Please check upper ranks of", res[duplicated(res)], "as they may have typos.\n"))
  #     res <- unique(res)
  #   }
  #   res
  # }

  # ## Replace all levels in taxonomy of non-top / non-unknown taxa to Other
  # tax <- as(tax_table(physeq), "matrix")
  # ii <- (tax[ , taxaRank2] == "Unknown") | (taxa_names(physeq) %in% topTaxa$taxa)
  # tax[!ii, ] <- "Other"
  # tax_table(physeq) <- tax
  # archetype <- taxa_names(physeq)[!ii][1]
  # physeq <- merge_taxa(physeq, eqtaxa = which(!ii), archetype = archetype)
  # taxa_names(physeq)[taxa_names(physeq) == archetype] <- "Other"
  # # physeq <- fast_tax_glom(physeq, taxrank = taxaRank2)
  #
  # tdf <- psmelt(physeq)
  # tdf[, taxaRank2] <- factor(tdf[, taxaRank2], levels = correct_taxonomy(topTaxa[[taxaRank2]]))
  #
  # ## Create correct order for the filling variable (if different from taxaRank2)
  # if (fill != taxaRank2) {
  #   fill_order <- count(topTaxa, .data[[fill]], wt = abundance, sort = TRUE) %>%
  #     pull(fill)
  #   tdf[, fill] <- factor(tdf[, fill], levels = correct_taxonomy(fill_order)) %>%
  #     droplevels()
  # }
  #
  # tdf <- tdf %>% arrange(desc(.data[[fill]]), desc(Abundance))
  #
  # return(tdf)
}

#' Plot a distance matrix as a heatmap
#'
#' @param dist Distance matrix between the samples
#' @param order Optional. Sample order used for plotting
#' @param title Optional. Plot title.
#' @param low Optional. Default "#B1F756". Color for low values
#' @param high Optional. Default "#132B13". Color for high values
#' @param show.names Optional. Default FALSE. Show sample names in heatmap. Use only
#'                   if there are not too many samples.
#'
#' @return A ggplot2 object
#' @export
#'
#' @importFrom dplyr as_tibble mutate
#' @importFrom ggplot2 aes element_blank element_text geom_tile ggplot ggtitle scale_fill_gradient theme
#' @importFrom tidyr pivot_longer
#' @examples
#' data(food)
#' dist.bc <- phyloseq::distance(food, "bray")
#' plot_dist_as_heatmap(dist.bc)
plot_dist_as_heatmap <- function(dist, order = NULL, title = NULL,
                                 low = "#B1F756", high = "#132B13",
                                 show.names = FALSE) {
  data <- dist %>%
    as.matrix() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(x = names(.)) %>%
    pivot_longer(cols = -x, values_to = "distance", names_to = "y")
  if (!is.null(order)) {
    data$x <- factor(data$x, levels = order)
    data$y <- factor(data$y, levels = order)
  }
  p <- ggplot(data, aes(x = x, y = y, fill = distance)) +
    geom_tile()
  p <- p + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90)
  )
  if (!show.names) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
  }
  p <- p + scale_fill_gradient(low = low, high = high)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

## Wrapper around hclust to represent clustering tree
## with leaves colored according to some variables
#' Wrapper around hclust to represent clustering tree with leaves colored according to a factor.
#'
#' @param physeq phyloseq class object
#' @param dist distance matrix (dist class) or character to be used in phyloseq::distance function
#' @param method (character) linkage method used in hclust, defaults to "ward.D2"
#' @param color (character) variable name used to color tree leaves. Defaults to NULL for black leaves
#' @param label (character) variable name used to label tree leaves. Defaults to NULL for sample names.
#' @param title (character) optional. Plot title, defaults to "method" clustering tree.
#' @param palette (named color vector) optional. Manual color palette
#' @param ... (optional) additional parameters passed on to theme(axis.text.x = element_text(...)) to control label size, justification, ...
#'
#' @return Nothing, used for its plotting side effect
#' @export
#'
#' @examples
#' data(food)
#' ## Basic plot
#' plot_clust(food, dist = "unifrac")
#' plot_clust(food, dist = "unifrac", color = "EnvType")
#' ## Slightly better plot
#' plot_clust(food, dist = "unifrac", color = "EnvType", label = "EnvType", size = 8) + theme(legend.position = "none")
#' @importFrom ape as.phylo
#' @importFrom dplyr as_tibble mutate
#' @importFrom ggplot2 aes labs scale_color_manual theme
#' @importFrom ggtext element_markdown
#' @importFrom ggtree geom_tiplab geom_tippoint ggtree layout_dendrogram
#' @importFrom phyloseq distance get_variable nsamples sample_data sample_variables
#' @importFrom scales col_factor hue_pal
#' @importFrom stats hclust
plot_clust <- function(physeq, dist, method = "ward.D2", color = NULL,
                       label = NULL,
                       title = paste(method, "linkage clustering tree"),
                       palette = NULL, ...) {
  # label
  if (is.null(label) || !is.character(label) || !(label %in% sample_variables(physeq, errorIfNULL = FALSE))) {
    label <- "label"
  }
  # color
  if (is.character(color)) {
    legend.title <- NULL
    color_var <- phyloseq::get_variable(physeq, color)
  } else {
    legend.title <- NULL
    color_var <- rep("black", nsamples(physeq))
  }
  color_var <- as.factor(color_var)
  color_levels <- levels(color_var)
  ## compute distance
  if (is.character(dist)) {
    dist <- dist[1]
    dist <- phyloseq::distance(physeq, method = dist)
  }
  ## automatic color palette: one color per different sample type
  if (is.null(palette)) {
    if (length(color_levels) == 1) {
      color_palette <- "black"
    } else {
      color_palette <- hue_pal()(length(color_levels))
    }
  } else {
    color_palette <- palette[color_levels]
  }
  ## Change hclust object to phylo object and plot
  clust <- as.phylo(hclust(dist, method = method))
  ## change tip label if needed
  # if (!is.null(label)) {
  #   tip.dict <- setNames(as.character(phyloseq::get_variable(physeq, label)),
  #                        sample_names(physeq))
  #   clust$tip.label <- tip.dict[clust$tip.label]
  # }
  ## extract metadata
  meta <- phyloseq::sample_data(physeq) %>%
    dplyr::as_tibble(rownames = "label")
  if (!is.null(color)) {
    meta <- dplyr::mutate(meta, tip_color = scales::col_factor(color_palette, levels = color_levels)(.data[[color]]))
  } else {
    meta <- dplyr::mutate(meta, tip_color = "black")
  }
  ## plot clustering tree
  ggtree::`%<+%`(ggtree::ggtree(clust), meta) +
    ggtree::layout_dendrogram() +
    ggtree::geom_tippoint(if (!is.null(color)) aes(color = .data[[color]])) +
    ## as_ylab allows to print tip labels as axis.tick
    ggtree::geom_tiplab(as_ylab = TRUE, aes(label = .data[[label]])) +
    scale_color_manual(values = color_palette) +
    ## the color of which can be changed via theme()
    theme(axis.text.x = ggtext::element_markdown(color = meta$tip_color, hjust = 1, vjust = 0.5, ...)) +
    labs(title = title)
}

## Extract legend from a ggplot object
#' @importFrom ggplot2 ggplot_build ggplot_gtable
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}