#' Build rarefaction curves from a count matrice
#'
#' @description Warning, this function is deprecated and will be phased out at some point.
#'
#' @param physeq      A phyloseq class object
#' @param step   Step size for sample sizes in rarefaction curves. Defaults to 1.
#' @param sample Subsample size for rarefying community, either a single value or a vector.
#' @param xlab,ylab   Axis labels in plots of rarefaction curves.
#' @param ylab
#' @param label       Label rarefaction curves by rownames of x (logical).
#' @param col         Custom colours for line rarefactions curves.
#' @param ...         parameters passed to ordilabels
#'
#' @return
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
