## Create ggplot style data frame for ternary plots from a phyloseq class objects
## Samples are grouped according to `group` and an error is returned is `group` has
## 5 or more levels
ternary_norm <- function(physeq, group, levelOrder = NULL, raw = FALSE, normalizeGroups = TRUE) {
    ## Args:
    ## - phyloseq class object, otus abundances are extracted from this object
    ## - group: Either the a single character string matching a
    ##          variable name in the corresponding sample_data of ‘physeq’, or a
    ##          factor with the same length as the number of samples in ‘physeq’.
    ## - raw: logical, should raw read counts be used to compute relative abudances of an
    ##        OTU among different conditions (defaults to FALSE)
    ## - levelOrder: Order along which to rearrange levels of `group`. Goes like (left, top, right) for
    ##               ternary plots and (left, top, right, bottom) for diamond plots. 
    ## - normalizeGroups: logical, only used if raw = FALSE, should all levels be given
    ##                    equal weights (TRUE, default) or weights equal to their sizes (FALSE)
    
    ## Get grouping factor 
    if (!is.null(sam_data(physeq, FALSE))) {
        if (class(group) == "character" & length(group) == 1) {
            x1 <- data.frame(sam_data(physeq))
            if (!group %in% colnames(x1)) {
                stop("group not found among sample variable names.")
            }
            group <- x1[, group]
        }
    }
    if (class(group) != "factor") {
        group <- factor(group)
    }

    ## Reorder levels of factor
    if (length(levels(group)) > 4) {
        warnings("There are 5 groups or more, the data frame will not be suitable for ternary plots.")
    }
    if (!is.null(levelOrder)) {
        if (any(! group %in% levelOrder)) {
            stop("Some levels of the factor are not included in `levelOrder`")
        } else {
            group <- factor(group, levels = levelOrder)
        }
    }
        
    ## construct relative abundances matrix
    tdf <- as(otu_table(physeq), "matrix")
    if (!taxa_are_rows(physeq)) { tdf <- t(tdf) }
    
    ## If raw, no normalisation should be done
    if (raw) {
        tdf <- t(tdf)
        abundance <- rowSums(t(tdf))/sum(tdf)
        meandf <- t(rowsum(tdf, group, reorder = TRUE))/rowSums(t(tdf))
    } else {        
        ## Construct relative abundances by sample
        tdf <- apply(tdf, 2, function(x) x/sum(x))
        if (normalizeGroups) {
            meandf <- t(rowsum(t(tdf), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(tdf)),
                                                  nrow = nrow(tdf))
            abundance <- rowSums(meandf)/sum(meandf)
            meandf <- meandf / rowSums(meandf)
        } else {
            abundance <- rowSums(tdf)/sum(tdf)
            meandf <- t(rowsum(t(tdf), group, reorder = TRUE))/rowSums(tdf)
        }
    }

    ## Construct cartesian coordinates for de Finetti's diagram
    ## (taken from wikipedia, http://en.wikipedia.org/wiki/Ternary_plot)
    if (ncol(meandf) == 3) {
        ternary.coord <- function(a,b,c) { # a = left, b = right, c = top
            return(data.frame(x = 1/2 * (2*b + c)/(a + b + c),
                              y = sqrt(3) / 2 * c / (a + b + c)))
        }
        cat(paste("(a, b, c) or (left, right, top) are (",
                  paste(colnames(meandf), collapse = ", "),
                  ")", sep = ""), sep = "\n")
        ## Data points
        df <- data.frame(x = 1/2 * (2*meandf[ , 2] + meandf[ , 3]),
                         y = sqrt(3)/2 * meandf[ , 3],
                         abundance = abundance, 
                         row.names = rownames(meandf))
        ## Extreme points
        extreme <- data.frame(ternary.coord(a = c(1, 0, 0),
                                            b = c(0, 1, 0),
                                            c = c(0, 0, 1)),
                              labels = colnames(meandf),
                              row.names = c("left", "right", "top"))
    }

    if (ncol(meandf) == 4) {
        diamond.coord <- function(a, b, c, d) {
            return(data.frame(x = (a - c) / (a + b + c + d),
                              y = (b - d) / (a + b + c + d)))
        }
        cat(paste("(a, b, c, d) or (right, top, left, bottom) are (",
                  paste(colnames(meandf), collapse = ", "),
                  ")", sep = ""), sep = "\n")
        ## data points
        df <- data.frame(x = (meandf[ , 1] - meandf[ , 3]),
                         y = (meandf[ , 2] - meandf[ , 4]),
                         abundance = abundance, 
                         row.names = rownames(meandf))
        ## extreme points
        extreme <- data.frame(diamond.coord(a = c(1, 0, 0, 0),
                                            b = c(0, 1, 0, 0),
                                            c = c(0, 0, 1, 0),
                                            d = c(0, 0, 0, 1)),
                              labels = colnames(meandf),
                              row.names = c("right", "top", "left", "bottom"))
    }

    ## Merge coordinates with taxonomix information
    df$otu <- rownames(df)
    ## Add taxonomic information
    if (!is.null(tax_table(physeq, FALSE))) {
        tax <- data.frame(otu = rownames(tax_table(physeq)),
                          tax_table(physeq))
        df <- merge(df, tax, by.x = "otu")
    }

    ## Add attributes
    attr(df, "labels") <- colnames(meandf)
    attr(df, "extreme") <- extreme
    attr(df, "type") <- c("ternary", "diamond", "other")[cut(ncol(meandf), breaks = c(0, 3, 4, Inf))]
    return(df)
}


ternary_plot <- function(physeq, group, grid = TRUE, size = "log2(abundance)",
                         color = NULL, shape = NULL, label = NULL,
                         levelOrder = NULL, plot = TRUE,
                         raw = FALSE, normalizeGroups = TRUE) {
    ## Args:
    ## - phyloseq class object, otus abundances are extracted from this object
    ## - group: Either the a single character string matching a
    ##          variable name in the corresponding sample_data of ‘physeq’, or a
    ##          factor with the same length as the number of samples in ‘physeq’.
    ## - raw: logical, should raw read counts be used to compute relative abudances of an
    ##        OTU among different conditions (defaults to FALSE)
    ## - normalizeGroups: logical, only used if raw = FALSE, should all levels be given
    ##                    equal weights (TRUE, default) or weights equal to their sizes (FALSE)
    ## - levelOrder: Order along which to rearrange levels of `group`. Goes like (left, top, right) for
    ##               ternary plots and (left, top, right, bottom) for diamond plots.
    ## - plot: logical, should the figure be plotted
    ## - grid: logical, should a grid be plotted.
    ## - size: mapping for size aesthetics, defaults to `abundance`.
    ## - shape: mapping for shape aesthetics.
    ## - color: mapping for color aesthetics.
    ## - label: Default `NULL`. Character string. The name of the variable
    ##          to map to text labels on the plot. Similar to color option
    ##          but for plotting text.
    data <- ternary_norm(physeq, group, levelOrder, raw, normalizeGroups)
    labels <- attr(data, "labels")
    extreme <- attr(data, "extreme")
    type <- attr(data, "type")

    if (type == "other") {
        stop("Ternary plots are only available for 3 or 4 levels")
    }
    
    ## borders
    borders <- data.frame(x = extreme$x,
                          y = extreme$y,
                          xend = extreme$x[c(2:nrow(extreme), 1)],
                          yend = extreme$y[c(2:nrow(extreme), 1)])
    ## grid
    ternary.coord <- function(a,b,c) { # a = left, b = right, c = top
        return(data.frame(x = 1/2 * (2*b + c)/(a + b + c),
                          y = sqrt(3) / 2 * c / (a + b + c)))
    }
    diamond.coord <- function(a, b, c, d) {
        return(data.frame(x = (a - c) / (a + b + c + d),
                          y = (b - d) / (a + b + c + d)))
    }
    x <- seq(1, 9, 1) / 10    
    
    ## Create base plot with theme_bw
    p <- ggplot() + theme_bw()
    ## Remove normal grid, axes titles and axes ticks
    p <- p + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), 
                   panel.border = element_blank(),
                   axis.ticks = element_blank(), 
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank())
    
    if (type == "ternary") {
        ## prepare levels' labels
        axes <- extreme
        axes$x <- axes$x + c(-1/2, 1/2, 0) * 0.1
        axes$y <- axes$y + c(-sqrt(3)/4, -sqrt(3)/4, sqrt(3)/4) * 0.1
        
        ## prepare ternary grid
        bottom.ticks <- ternary.coord(a = x, b = 1-x, c = 0)
        left.ticks <- ternary.coord(a = x, b = 0, c = 1-x)
        right.ticks <- ternary.coord(a = 0, b = 1 - x, c = x)
        ticks <- data.frame(bottom.ticks, left.ticks, right.ticks)
        colnames(ticks) <- c("xb", "yb", "xl", "yl", "xr", "yr")
        
        ## Add grid (optional)
        if (grid == TRUE) {
            p <- p + geom_segment(data = ticks, aes(x = xb, y = yb, xend = xl, yend = yl),
                                  size = 0.25, color = "grey40")
            p <- p + geom_segment(data = ticks, aes(x = xb, y = yb, xend = xr, yend = yr),
                                  size = 0.25, color = "grey40")
            p <- p + geom_segment(data = ticks, aes(x = rev(xl), y = rev(yl), xend = xr, yend = yr),
                                  size = 0.25, color = "grey40")
        }
    }

    if (type == "diamond") {
        ## prepare levels' labels
        axes <- extreme
        axes$x <- axes$x + c(1, 0, -1, 0) * 0.1
        axes$y <- axes$y + c(0, 1, 0, -1) * 0.1
        
        ## prepare diamond grid 
        nw.ticks <- diamond.coord(a = x, b = 1-x, c = 0, d = 0)
        ne.ticks <- diamond.coord(a = 0, b = x, c = 1-x, d = 0)
        sw.ticks <- diamond.coord(a = x, b = 0, c = 0, d = 1 - x)
        se.ticks <- diamond.coord(a = 0, b = 0, c = 1-x, d = x)
        ticks <- data.frame(nw.ticks, ne.ticks, se.ticks, sw.ticks)
        colnames(ticks) <- c("xnw", "ynw", "xne", "yne",
                             "xse", "yse", "xsw", "ysw")        
        ## Add grid (optional)
        if (grid == TRUE) {
            p <- p + geom_segment(data = ticks, aes(x = xnw, y = ynw, xend = xse, yend = yse),
                                  size = 0.25, color = "grey40")
            p <- p + geom_segment(data = ticks, aes(x = xne, y = yne, xend = xsw, yend = ysw),
                                  size = 0.25, color = "grey40")
            p <- p + geom_segment(aes(x = c(0, -1), y = c(-1, 0),
                                      xend = c(0, 1), yend = c(1, 0)),
                                  size = 0.25, color = "grey40")
        }
    }
    
    ## Add borders
    p <- p + geom_segment(data = borders, aes(x = x, y = y, xend = xend, yend = yend))
    ## Add levels' labels
    p <- p + geom_text(data = axes, aes(x = x, y = y, label = labels))
    
    ## Add, any custom-supplied plot-mapped variables
    if( length(color) > 1 ){
        data$color <- color
        names(data)[names(data)=="color"] <- deparse(substitute(color))
        color <- deparse(substitute(color))
    }
    if( length(shape) > 1 ){
        data$shape <- shape
        names(data)[names(data)=="shape"] <- deparse(substitute(shape))
        shape <- deparse(substitute(shape))
    }	
    if( length(label) > 1 ){
        data$label <- label
        names(data)[names(data)=="label"] <- deparse(substitute(label))
        label <- deparse(substitute(label))
    }
    if( length(size) > 1 ){
        data$size <- size
        names(data)[names(data)=="size"] <- deparse(substitute(size))
        size <- deparse(substitute(size))
    }

    ## Add data points
    ternary_map <- aes_string(x = "x", y = "y", color = color,
                              shape = shape, size = size, na.rm = TRUE)
    p <- p + geom_point(data = data, mapping = ternary_map)

    ## Add the text labels
    if( !is.null(label) ){
        label_map <- aes_string(x="x", y="y", label=label, na.rm=TRUE)
        p <- p + geom_text(data = data, mapping = label_map,
                           size=3, vjust=1.5, na.rm=TRUE)
    }

    if (plot) {
        plot(p)
    }
    
    invisible(p)   
}
