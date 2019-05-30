#' Build rarefaction curves from a count matrice
#'
#' @description Warning, this function is deprecated and will be phased out at some point.
#'
#' @param physeq A phyloseq class object
#' @param step   Step size for sample sizes in rarefaction curves. Defaults to 1.
#' @param sample Subsample size for rarefying community, either a single value or a vector.
#' @param xlab   Axis labels in plots of rarefaction curves.
#' @param ylab   Axis labels in plots of rarefaction curves.
#' @param label  Label rarefaction curves by rownames of x (logical).
#' @param col    Custom colours for line rarefactions curves.
#' @param ...    parameters passed to ordilabels
#'
#' @export
#'
#' @importFrom vegan specnumber ordilabel
#' @examples
#' data(food)
#' rarecurve2(food)
#'
rarecurve2 <- function (physeq, step = 1, sample, xlab = "Sample Size", ylab = "Number of species", label = TRUE, col = "black", ...)
{
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) {
    x <- t(x)
  }
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
