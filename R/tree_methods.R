## Methods for representing differentially abundant (across conditions) otus on a grid or on a tree.
## Also useful for general sets of interesting otus. 

## Small multiple plot: each otu in set daOTUS has its own plot and each
## its (mean) relative abundance in different levels of condition variable
## is represented as a point. If replicates are available for a given level,
## SEM are represented as error bars.
daOTU_plot <- function(physeq, daOTUS, variable, title = NULL, plot = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object
  ## - daOTUS: set of OTUS to plot
  ## - variable: x axis, used for plotting and aggregating samples
  ## - title: optional, plot title
  ##
  ## Returns,
  ## - p: invisible, ggplot object
  ##
  ## check variable and otu set
  numberOfProvidedTaxa <- length(daOTUS)
  daOTUS <- intersect(daOTUS, taxa_names(physeq))
  stopifnot(variable %in% sample_variables(physeq), length(daOTUS) > 0)
  if ( numberOfProvidedTaxa != length(daOTUS) ) {
    cat(paste(length(daOTUS), "otus out of", numberOfProvidedTaxa,
              "provided were found, use only those for plotting"), sep = "\n")
  }
  ## Get count data set
  tdf <- otu_table(physeq)
  if (!taxa_are_rows(physeq)) { tdf <- t(tdf) }
  ## Normalize counts to relative abundances
  tdf <- apply(tdf, 2, function(x) x/sum(x))
  tdf <- tdf[daOTUS, ]
  ## sort OTUs by overall abundance
  otuAbun <- sort(rowSums(tdf), TRUE)
  tdf <- data.frame(otu = rownames(tdf), tdf[daOTUS, ])
  ## Get within group mean and sem
  tdf <- melt(tdf)
  colnames(tdf)[2] <- "Sample"
  tdf <- merge(tdf,
               data.frame(Sample = sample_names(physeq), Replicate = get_variable(physeq, variable)),
               by.y = "Sample")
  meandf <- ddply(tdf, .(otu, Replicate), summarize,
                  mean = mean(value), sem = sd(value)/sqrt(length(value)))
  colnames(meandf)[ colnames(meandf) == "Replicate"] <- variable
  ## Get best taxonomic annotation
  bestAnn <- function(x) {
    x <- x[!is.na(x)]
    return(ifelse(length(x), x[length(x)], ""))
  }
  taxdf <- as(tax_table(physeq), "matrix")
  ann <- apply(taxdf, 1, bestAnn)
  taxdf <- data.frame(otu = rownames(taxdf), taxdf, tax = ann)[daOTUS, ]
  taxdf$otu2 <- paste(taxdf$otu, taxdf$tax)
  mdf <- merge(meandf, taxdf, by.x = "otu")
  ## Sort results by overall abundance
  mdf$otu2 <- factor(mdf$otu2, levels = taxdf[names(otuAbun), "otu2"])
  mdf[[variable]] <- factor(mdf[[variable]], levels = levels(get_variable(physeq, variable)))
  ## Plot results
  p <- ggplot(mdf, aes_string(x = variable, y = "mean", color = variable))
  p <- p + geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1) + geom_point()
  ## Add custom title if provided
  if (is.null(title)) {
    title.text <- paste("Differentially Abundant OTUs (", variable, ")", sep = "")
  } else {
    title.text <- title
  }
  p <- p + facet_wrap(~otu2, scales = "free_y") + ggtitle(title)
  if (plot) {
      plot(p)
  }
  return(invisible(p))
}

## subsidiary function to daOTU_tree_plot
## Add thermometer to figure
addThermo <- function(lastPP, meandf, adj.x, width.thermo = ncol(meandf)*8*strwidth(" ")) {
  ## Add Thermometer
  BOTHlabels(XX = max(lastPP$xx[1:lastPP$Ntip] + adj.x) + 0.3*8*strwidth(" ") + 0.5*width.thermo,
             YY = lastPP$yy[1:lastPP$Ntip], thermo = meandf,
             horiz = TRUE, height = 0.9, adj = c(0.5, 0.5), pch = NULL,
             width = width.thermo, piecol = gg_color_hue(ncol(meandf)),
             frame = "none", pie = NULL, color = "black")
  return(invisible(width.thermo))
}

## subsidiary function to daOTU_tree_plot
## Add Histogram to figure
addHist <- function(lastPP, meandf, adj.x, width.thermo = ncol(meandf)*8*strwidth(" "), cex = CEX) {
  ## Add squares
  nSquares <- ncol(meandf)
  XX <- max(lastPP$xx[1:lastPP$Ntip] + adj.x) + 0.3*8*strwidth(" ")
  YY = lastPP$yy[1:lastPP$Ntip] - 1/2 * strheight(" ", cex = cex)
  base.size <- width.thermo/nSquares ## baseline for histogram bin width
  height.size <- 0.9
  xl <- rep(XX + (seq(nSquares) - 1)*base.size, times = nrow(meandf))     ## left x positions
  xr <- xl + base.size ## right x positions
  yb <- rep(YY, each = ncol(meandf))     ## bottom y positions
  yt <- as.vector(yb + t(height.size*meandf))     ## top y positions
  rect(xl, yb, xr, yt, col = gg_color_hue(ncol(meandf)))
  return(invisible(width.thermo))
}

## subsidiary function to daOTU_tree_plot
## Add Rectangles to figure
addRectangle <- function(lastPP, meandf, adj.x, width.thermo = ncol(meandf)*8*strwidth(" ")) {
  ## Add rectangles
  nRectangles <- ncol(meandf)
  XX <- max(lastPP$xx[1:lastPP$Ntip] + adj.x) + 0.3*8*strwidth(" ")
  YY = lastPP$yy[1:lastPP$Ntip]
  base.size <- width.thermo/nRectangles ## baseline for rectangles side length
  height.size <- 0.9
  scaling <- sqrt(meandf)
  xcenter <- rep(XX + (seq(nRectangles) - 1/2)*base.size, times = nrow(scaling))
  ycenter <- rep(YY, each = ncol(scaling))
  rectangles <- cbind(as.vector(base.size*t(scaling)), as.vector(height.size*t(scaling)))
  symbols(xcenter, ycenter, rectangles = rectangles,
          inches = FALSE, bg = gg_color_hue(ncol(scaling)), add = TRUE)
  return(invisible(widh.thermo))
}

## subsidiary function to daOTU_tree_plot
## Add circles to figure
addCircle <- function(lastPP, meandf, adj.x, width.thermo = ncol(meandf)*8*strwidth(" ")) {
  ## Add circles
  ## Requires plotrix
  nCircles <- ncol(meandf)
  XX <- max(lastPP$xx[1:lastPP$Ntip] + adj.x) + 0.3*8*strwidth(" ")
  YY = lastPP$yy[1:lastPP$Ntip]
  ## baseline for circles radii
  ratio <- diff(lastPP$y.lim)/diff(lastPP$x.lim)
  base.size <- min(width.thermo/nCircles, 0.95/(ratio * sqrt(2))) 
  scaling <- sqrt(meandf)
  xcenter <- rep(XX + (seq(nCircles) - 1/2)*base.size, times = nrow(scaling))
  ycenter <- rep(YY, each = ncol(scaling))
  ## Area proportional to frequency, sqrt(2) used to maximize space occupation while
  ## avoiding overlap (worst case scenario is two neighbor modalities with frequencies of
  ## 0.5 and 0.5)
  circles <- 0.5 * as.vector(base.size*t(scaling))
  symbols(xcenter, ycenter, circles = circles,
          inches = FALSE, bg = gg_color_hue(ncol(scaling)), add = TRUE)
  width.thermo <- min(width.thermo, nCircles * base.size)
  return(invisible(width.thermo))
}

## subsidiary function to daOTU_tree_plot
## Add squares to figure
addSquare <- function(lastPP, meandf, adj.x, width.thermo = ncol(meandf)*8*strwidth(" ")) {
  ## Add squares
  ## Requires plotrix
  nSquares <- ncol(meandf)
  XX <- max(lastPP$xx[1:lastPP$Ntip] + adj.x) + 0.3*8*strwidth(" ")
  YY = lastPP$yy[1:lastPP$Ntip]
  ## baseline for squares radii
  ratio <- diff(lastPP$y.lim)/diff(lastPP$x.lim)
  base.size <- min(width.thermo/nSquares, 0.95/ratio) 
  scaling <- sqrt(meandf)
  xcenter <- rep(XX + (seq(nSquares) - 1/2)*base.size, times = nrow(scaling))
  ycenter <- rep(YY, each = ncol(scaling))
  squares <- as.vector(base.size*t(scaling)) / sqrt(2)
  symbols(xcenter, ycenter, squares = squares,
          inches = FALSE, bg = gg_color_hue(ncol(scaling)), add = TRUE)
  width.thermo <- min(width.thermo, nSquares * base.size)
  return(invisible(width.thermo))
}


## subsidiary function to daOTU_tree_plot
## Add abundance, automatically calibrate scale to roughly match 1:4
addAbundance <- function(otuAbun, lastPP, adj.x, width.thermo,
                         abundance.type = c("circle", "bar"),
                         add.ruler = FALSE) {
  scaling <- log(otuAbun, base = 10)
  minAbun <- floor(min(scaling))
  maxAbun <- floor(max(scaling)) ## Only to construct scale and breaks
  if (maxAbun==minAbun) {
      scale <- minAbun
  } else {
      scale <- seq(minAbun, maxAbun, by = floor((maxAbun - minAbun)/4) + 1)
  }
  breaks <- 10^scale
  maxAbun <- max(scaling) ## to properly scale all quantities
  abundance.type <- match.arg(abundance.type)
  if (abundance.type == "circle") {
      ## Maximum and minimum circle radius
      ratio <- diff(lastPP$y.lim)/diff(lastPP$x.lim)
      max.size <- 0.95 /(ratio * 2 * sqrt(2))
      min.size <- max.size / 20
      ## Scaling radius
      scaling <- sqrt( (scaling - minAbun) / (maxAbun - minAbun))
      scaling <- scaling * (max.size - min.size) + min.size
      ## Construct scale
      scale <- sqrt( (scale - minAbun) / (maxAbun - minAbun))
      scale <- scale * (max.size - min.size) + min.size
      ## Centre coordinates
      adj.radii <- max(scaling)
      xcenter <- rep(max(lastPP$xx[1:lastPP$Ntip] + adj.x) + 0.5*8*strwidth(" ") + width.thermo + adj.radii,
                     lastPP$Ntip)
      ycenter <- lastPP$yy[1:lastPP$Ntip]
      ## Plot circles or segments
      symbols(xcenter, ycenter, circles = scaling,
              inches = FALSE, bg = "grey40", fg = "grey40", add = TRUE)
  } else {
      ## Min and max widths
      if (width.thermo >= 4*8*strwidth(" ")) {
          max.size <- width.thermo / 2
      } else  {
          max.size <- width.thermo 
      }
      min.size <- max.size / 20
      ## Scaling widths
      scaling <- (scaling - minAbun) / (maxAbun - minAbun) * (max.size - min.size) + min.size
      ## Construct scale
      scale <- (scale - minAbun) / (maxAbun - minAbun) * (max.size - min.size) + min.size
      XX <- rep(max(lastPP$xx[1:lastPP$Ntip] + adj.x) + 0.5*8*strwidth(" ") + width.thermo,
                lastPP$Ntip) + scaling / 2
      YY <- lastPP$yy[1:lastPP$Ntip]
      symbols(XX, YY, inches = FALSE, rectangles = cbind(scaling, 0.5),
              add = TRUE, fg = "grey40", bg = "grey40")
      if (add.ruler) {
          offset <- max(lastPP$xx[1:lastPP$Ntip] + adj.x) + 0.5 * 8 * strwidth(" ") + width.thermo
          segments(x0 = offset + scale,
                   y0 = rep(lastPP$y.lim[1] - 0.5, length(scale)),
                   x1 = offset + scale,
                   y1 = rep(lastPP$y.lim[2] + 0.5, length(scale)), col = "grey80")
          text(x = offset + scale,
               y = lastPP$y.lim[1] - 1,
               adj = c(0, 0.5),
               labels = prettyNum(breaks),
               srt = -90, cex = 0.7)
      }
  }
  return(invisible(list(scale = scale, type = abundance.type,
                        breaks = breaks)))
}

## subsidiary function to daOTU_tree_plot
## Add legend
addLegend <- function(meandf, variable, scale, type, breaks) {
  leg <- legend(x = "topleft",
                xjust = 0,
                legend = colnames(meandf),
                fill = gg_color_hue(ncol(meandf)),
                border = NA,
                title = variable,
                bty = "n")
  leg.abundance <- legend(x = leg$rect$left,
                          y = leg$rect$top - 1.1*leg$rect$h, 
                          col = "grey40",
                          legend =  prettyNum(breaks),
                          bty = "n", 
                          title = "Relative abundance",
                          plot = TRUE)
  if (type == "circle") {
      xcenter <- leg.abundance$text$x - max(scale) - strwidth(" ")
      ycenter <- leg.abundance$text$y
      symbols(xcenter, ycenter, circles = scale,
              inches = FALSE, bg = "grey40", fg = "grey40", add = TRUE)
  } else {
      xcenter <- leg.abundance$text$x - strwidth(" ") - scale / 2
      ycenter <- leg.abundance$text$y
      symbols(xcenter, ycenter, rectangles = cbind(scale, 0.5),
              inches = FALSE, bg = "grey40", fg = "grey40", add = TRUE)
  }
}


## Tree plot with thermometer: plot a tree with all otus in set daOTUS and add
## overall abundance and distribution between levels of variable at tips of tree. 
## Samples are first agreggated by level of variable to compute relative abundance
## within this level and different levels are then agreggated (and eventually weighted)
## to determine which levels contribute most to the overall relative abundance. Overall
## abundance is also plotted (in logarithmic scale) next to the thermometer
daOTU_tree_plot <- function(physeq, daOTUS, variable, colorTip = FALSE, 
                            type = c("thermo", "circle", "square", "rectangle", "hist"),
                            abundance.type = c("bar", "circle"),
                            add.ruler = TRUE, 
                            weight = NULL, ...) {
  ## Args:
  ## - physeq: phyloseq class object
  ## - daOTUS: set of OTUS to plot
  ## - variable: x axis, used for plotting and aggregating samples
  ## - type: type of graphic to plot on leaves of tree
  ## - abundance.type: type of symbol to use for abundance plotting
  ## - add.ruler: logical, should a ruler be added to abundance (only when abundance.type is bar)
  ## - colorTip: logical, should tips be colored according to dominant condition
  ##             (defaults to FALSE)
  ## - weight: optional, contribution of each factor level of variable to the whole
  ## - ...: Additional arguments passed to ‘plot.phylo’.
  ##
  ## Returns,
  ## - plot of tree, thermometer and relative abundance
  ##
  ## check variable and otu set
  numberOfProvidedTaxa <- length(daOTUS)
  daOTUS <- intersect(daOTUS, taxa_names(physeq))
  stopifnot(variable %in% sample_variables(physeq),
            length(daOTUS) > 0)
  if ( numberOfProvidedTaxa != length(daOTUS) ) {
    cat(paste(length(daOTUS), "otus out of", numberOfProvidedTaxa,
              "provided were found, use only those for plotting"), sep = "\n")
  }
  ## Get count data frame
  tdf <- otu_table(physeq)
  if (!taxa_are_rows(physeq)) { tdf <- t(tdf) }
  ## Normalize counts to relative abundances
  tdf <- apply(tdf, 2, function(x) x/sum(x))
  tdf <- data.frame(tdf[daOTUS, ])
  ## sort OTUs by overall abundance
  otuAbun <- sort(rowMeans(tdf), TRUE)
  ## Get within group mean
  meandf <- aggregate(t(tdf), by = list(variable = get_variable(physeq, variable)), FUN = mean)
  groupLevels <- rownames(meandf) <- meandf[ , 1]
  meandf <- meandf[ , -1]
  ## Normalize by otu (use weight if necessary)
  if (is.null(weight)) {
    meandf <- t(apply(meandf, 2, function(x) x/sum(x)))
  } else {
    if (!is.null(names(weight)) && all(groupLevels %in% names(weight))) {
      w <- weight[groupLevels]
    } else {
      cat("No name or non-matching names for weight, using the first weights.",
          paste("Group names are", groupLevels), sep = "\n")
      w <- weight[length(groupLevels)]
    }
    meandf <- t(apply(meandf, 2, function(x) w*x/sum(w*x)))
  }
  ## Get best taxonomic annotation
  bestAnn <- function(x) {
    x <- x[!is.na(x)]
    return(ifelse(length(x), x[length(x)], ""))
  }
  taxdf <- as(tax_table(physeq), "matrix")
  ann <- apply(taxdf, 1, bestAnn)
  taxdf <- data.frame(otu = rownames(taxdf), tax = paste(rownames(taxdf), ann))
  rownames(taxdf) <- taxdf$otu
  ## get pruned tree
  phy <- phy_tree(physeq)
  phy <- drop.tip(phy, tip = setdiff(phy$tip.label, daOTUS))
  ## Reorder taxdf, meandf, otuAbun to match phy tips and add trailing
  ## spaces at end of phy tip labels
  taxdf <- taxdf[phy$tip.label, ]
  meandf <- meandf[phy$tip.label, ]
  otuAbun <- otuAbun[phy$tip.label]
  trailing.space <- paste(rep(" ", 8*ncol(meandf)*ifelse(ncol(meandf) >= 4, 1.5, 2)), collapse = "")
  phy$tip.label <- paste(taxdf$tax, trailing.space, sep = "")
  ## Plot results
  plot.phylo(phy, tip.color = "white", ...)
  args <- list(...)
  CEX <- ifelse("cex" %in% names(args), args$cex, par("cex"))
  if (colorTip) {
      tipColor <- gg_color_hue(ncol(meandf))[ apply(meandf, 1, which.max) ]
      tiplabels(paste("", taxdf$tax), frame = "n", adj = c(0, 0.5), col = tipColor, cex = CEX)
  } else {
      tiplabels(paste("", taxdf$tax), frame = "n", adj = c(0, 0.5), cex = CEX)
  }
  width.thermo <- strwidth(paste(rep(" ", 8*ncol(meandf)), collapse = ""))
  adj.x <- strwidth(as.character(taxdf$tax), units = "user", cex = CEX)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  ## Add Graphics at leaves of tree
  type <- match.arg(type)
  width.thermo <- switch(type,
                         thermo    = addThermo(lastPP, meandf, adj.x, width.thermo),
                         hist      = addHist(lastPP, meandf, adj.x, width.thermo, CEX),
                         square    = addSquare(lastPP, meandf, adj.x, width.thermo),
                         rectangle = addRectangle(lastPP, meandf, adj.x, width.thermo),
                         circle    = addCircle(lastPP, meandf, adj.x, width.thermo)
                         )
  ## Add overall average abundance
  abundance.type <- match.arg(abundance.type)
  abunScale <- addAbundance(otuAbun, lastPP, adj.x, width.thermo, abundance.type, add.ruler)
  ## Add segments
  XX <- lastPP$xx[1:lastPP$Ntip] + strwidth(as.character(taxdf$tax), cex = CEX)
  YY <- lastPP$yy[1:lastPP$Ntip] 
  segments(x0 = XX, x1 = max(XX),
           y0 = YY, y1 = YY, col = "grey80", lty = 2)
  ## Add legend
  addLegend(meandf, variable, abunScale[[1]], abunScale[[2]], abunScale[[3]])
}



#' ggplot-like implementation of plot.phylo from the ape package for \code{phyloseq} class
#' objects. Allows more representation that ggphylo and plot_tree. Returns a data frame representing
#' the tree that is useful for ggplot graphical function. 
#'
#' @title ggtree
#' @param physeq Required. An instance of a \code{phyloseq} class object that has
#'               tree slot. 
#' @param group  Optional. Either a single character string matching a variable
#'               name in the corresponding sam_data of `physeq`, or a factor with
#'               the same length as the number of taxa in `physeq`. If provided, the samples
#'               are aggregated by factor levels before being plotted
#' @param freq    Logical. Should counts be replaced by frequencies (TRUE) or kept as such (FALSE).
#'                 Defaults to TRUE. TRUE corresponding to weighting all samples equally while FALSE
#'                 weights them with their library size.
#' @param method Optional. Edge fattening method; either 'linear' or
#'                    'logarithmic' (or an abbrevation of these). Defaults to "linear"
#' @param width.lim Optional. Numeric. Minimum and maximal edge after fattening.
#'                  Defaults to c(0.1, 4). Set to c(0, x) for true linear scaling.
#' @param base      Optional. Numeric. Base used for logarithmic scaling. Defaults to 10.
#' @param edge.group  Optional. Either a single character string matching a variable
#'               name in the corresponding tax_table of `physeq`, or a factor with
#'               the same length as the number of taxa in `physeq`. Defaults to "Phylum"
#' @param edge.method Optional. Ancestral group reconstruction method; either 'majority' or
#'                    'ace' (or an abbrevation of these). Defaults to "majority"
#' @param ...    Optional. Additional arguments passed on to plot.phylo to determine
#'               the shape of the tree. 
#' @note ggtree assumes that the tree is rooted. 
#' @return A list of two data frames properly formatted for ggplot. 
ggtree <- function(physeq, group = NULL, freq = TRUE, method = "linear", base = 10,
                              edge.group = "Phylum", edge.method = "majority",
                              width.lim = c(0.1, 4), ...) {
    ## Exception handling
    if (is.null(phy_tree(physeq, FALSE))) {
        stop("Object \"physeq\" does not have a tree slot")
    } else {
        tree <- phy_tree(physeq)
    }
    if (!is.rooted(tree)) {
        stop("Tree must be rooted, consider using midpoint rooting")
    }

    ## Potentially normalize counts to frequencies
    if (freq) {
        physeq <- transform_sample_counts(physeq, function(x) x /sum(x))
    }
    
    ## Merge samples by grouping factor.
    if (!is.null(group)) {
        physeq <- merge_samples(physeq, group)
    }

    ## Get edge group factor
    if (!is.null(edge.group)) {
        x <- color_edges(physeq, edge.group, edge.method)
        edge.group <- factor(names(x$palette)[match(x$edge, x$palette)])
        tip.group <- factor(names(x$palette)[match(x$tip, x$palette)])
    }

    ## Fatten edges and compute average 
    fattened.edges <- fattenEdges(physeq, method, width.lim, base)$edge.raw.width
    mean.fattened.edges <- colMeans(fattened.edges)
    edge.names <- colnames(fattened.edges)
    fattened.edges <- data.frame(edge.name = edge.names,
                                 t(fattened.edges))
    mean.fattened.edges <- data.frame(edge.name = edge.names,
                                      mean.width = mean.fattened.edges)
    fattened.edges <- melt(fattened.edges, id.vars = "edge.name")
    colnames(fattened.edges)[2:3] <- c("Sample", "width")
    fattened.edges <- merge(fattened.edges, mean.fattened.edges)
    
    ## Call phylo.plot
    plot.phylo(phy_tree(physeq), plot = FALSE, ...)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    ggdata <- data.frame(edge.name = edge.names,
                         edge.x1   = lastPP$xx[ lastPP$edge[, 1]],
                         edge.x2   = lastPP$xx[ lastPP$edge[, 2]],
                         edge.y1   = lastPP$yy[ lastPP$edge[, 1]],
                         edge.y2   = lastPP$yy[ lastPP$edge[, 2]])
    ggdata$edge <- lastPP$edge
    if (!is.null(edge.group)) {
        edge.group <- data.frame(edge.name = edge.names, edge.group = edge.group)
        ggdata <- merge(ggdata, edge.group)
    }
    ggdata <- merge(ggdata, fattened.edges)
    if (lastPP == "cladogram") {
        p <- ggplot(ggdata) + geom_segment(aes(x = edge.x1, y = edge.y1, xend = edge.x2, yend = edge.y2))
    }
    if (lastPP == "phylogram") {
        p <- ggplot(ggdata) + geom_segment(aes(x = edge.x1, y = edge.y1, xend = edge.x1, yend = edge.y2))
        p <- p + geom_segment(aes(x = edge.x1, y = edge.y2, xend = edge.x2, yend = edge.y2))
    }
    if (lastPP == "radial") {
        
    }
    if (lastPP == "") {}
}
