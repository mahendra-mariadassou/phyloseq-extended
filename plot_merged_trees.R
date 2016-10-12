#' Functions for plotting trees with colored edges and tips according to data stored in
#' tax_table and otu_table of corresponding phyloseq class object

library(scales)

#' Fatten edges of a tree monotonically with the number of reads originating from that edge
#' in a community matrix. 
#'
#' @title fattenEdges
#' @param physeq Required. An instance of a \code{phyloseq} class object with otu_table and
#'               phy_tree slots.
#' @param method Optional. Edge fattening method; either 'linear' or
#'                    'logarithmic' (or an abbrevation of these). Defaults to "linear"
#' @param width.lim Optional. Numeric. Minimum and maximal edge after fattening.
#'                  Defaults to c(0.1, 4). Set to c(0, x) for true linear scaling.
#' @param base      Optional. Numeric. Base used for logarithmic scaling. Defaults to 10.
#' @param deviation Optional. Logical. Should the frequency (or count) of a taxa in a sample
#'                  be normalized with respect to the average of that taxa over all samples
#'                  before fattening. Defaults to FALSE. If method = "linear", the values are
#'                  centered around the mean (deviation). If method = "logarithmic", they are
#'                  divided by the mean (relative deviation)
#' @note width.lim behaves differently in "linear" and "logarithmic" modes. In linear mode, it is
#'       used to scale edges widths. In logarithmic mode, widths are only scaled if they
#'       go out of the range specified by `width.lim`.
#' @return A list with components
#' \item{edge.width} A matrix of size nsamples(physeq) x number of edges in \code{phy_tree(physeq)}
#'         with fattened edge widths
#' \item{edge.raw.width} A matrix of size nsamples(physeq) x number of edges in \code{phy_tree(physeq)}
#'         with raw edge widths (untransformed counts, frequencies or deviations)
#' \item{edge.presence} A logical matrix of size nsamples(physeq) x number of edges in
#'        \code{phy_tree(physeq)}. edge.presence[i, j] is TRUE if edge j is found in community i.
#' \item{tip.presence} A logical matrix of size nsamples(physeq) x ntaxa(physeq) in
#'        tip.presence[i, j] is TRUE if taxa j is found in community i.
#' \item{pendant.edges} A logical vector indicating which edges are pendant (i.e lead to a tip)
#' \item{legend} A named vector featuring edge widths and labels for legend. 
fattenEdges <- function(physeq, method = c("linear", "logarithmic"),
                        width.lim = c(0.1, 4), base = 10, deviation = FALSE) {
    ## Exception handling
    if (deviation & nsamples(physeq) == 1) {
        stop("There is only one sample (after potential aggregation),\n
              use of deviation is meaningless")
    }

    ## Extract otu_table matrix
    x <- as(otu_table(physeq), "matrix")
    if (taxa_are_rows(physeq)) { x <- t(x) }
    phy <- phy_tree(physeq)
    
    ## Scale counts to frequencies
    x <- x/rowSums(x)
    
    ## Construct incidence matrix of the tree
    incidence <- incidenceMatrix(phy)
        
    ## Order community table according to edge order and create
    ## community phylogenetic matrix
    x <- x[ , rownames(incidence), drop = FALSE]
    cpm <- x %*% incidence

    method <- match.arg(method)
    ## If needed, transform counts to deviations
    if (deviation) {
        if (method == "linear") {
            width <- sweep(cpm, 2, colMeans(cpm), "-")
            leg <- cbreaks(range(width))
        } else {
            cpm <- sweep(cpm, 2, colMeans(cpm), "/")
            leg <- symmetric_log_breaks(base = base)(cpm[cpm > 0])
            leg <- list(labels = prettyNum(leg), breaks = leg)
            epsilon <- min(cpm[cpm > 0]) ## small offset
            width <- log(cpm / epsilon, base = base)
            width[width == -Inf] <- 0 ## Is it wise?
            leg$breaks <- log(leg$breaks / epsilon, base = base)
        }
    }  else {
        if (method == "linear") { ## Linear scaling
            leg <- cbreaks(range(cpm))
            max.cpm <- max(cpm)
            ## Manual rescaling rather than with rescale
            width <- width.lim[1] + (width.lim[2] - width.lim[1]) * (cpm / max.cpm)
            leg$breaks <-
                width.lim[1] + (width.lim[2] - width.lim[1]) * (leg$breaks / max.cpm)
        } else { ## logarithmic scaling
            leg <- log_breaks(base = base)(cpm[cpm > 0])
            leg <- list(labels = prettyNum(leg), breaks = leg)
            epsilon <- min(cpm[cpm > 0]) ## small offset
            width <- log(cpm / epsilon, base = base)
            width[width == -Inf] <- NA
            min.width <- min(width, na.rm = TRUE)
            width <- width - min.width
            width[is.na(width)] <- 0
            leg$breaks <- log(leg$breaks / epsilon, base = base) - min.width
            max.width <- max(width)
            ## rescale widths if outside of range
            if (max.width + width.lim[1] > width.lim[2]) {
                width <- width.lim[1] + (width.lim[2] - width.lim[1]) * (width / max.width)
                leg$breaks <-
                    width.lim[1] + (width.lim[2] - width.lim[1]) * (leg$breaks / max.width)
            } else {
                width <- width + width.lim[1]
                leg$breaks <- leg$breaks + width.lim[1]
            }
            ## check leg$breaks for consistency
            leg$breaks[leg$breaks < 0] <- width.lim[1]
        }
    }
    return(list(edge.width = width,
                edge.raw.width = cpm,
                edge.presence = (cpm > 0),
                tip.presence  = (x > 0),
                pendant.edges = attr(incidence, "pendant.edges"),
                legend = leg))
}


#' Color edges and leaves of a tree according to a `group` factor. `group` should give a class
#' for all leaves. The class of inner nodes is inferred using majority rule or ancestral character
#' estimation (via ape::ace function, used with ER model)
#'
#' @title color_edges
#' @param physeq Required. An instance of a \code{phyloseq} class object that has
#'               tree slot. 
#' @param group  Optional. Either a single character string matching a variable
#'               name in the corresponding tax_table of `physeq`, or a factor with
#'               the same length as the number of taxa in `physeq`. Defaults to "Phylum"
#' @param method      Optional. Ancestral group reconstruction method; either 'majority' or
#'                    'ace' (or an abbrevation of these). Defaults to "majority"
#' @param tip.only  Optional. Logical. Should colors be computed for tips only? Defaults to FALSE.
#' @note The function assumes that the tree is rooted. The color palette is constructed
#'       automatically according to ggplot discrete hue scale 
#' @return A list with components
#' \item{edge} A color vector for edges of \code{phy_tree{physeq}}
#' \item{tip} A color vector for tip labels of \code{phy_tree{physeq}}
#' \item{palette} The named color palette color, useful for the legend
color_edges <- function(physeq,
                        group = "Phylum",
                        method = c("majority", "ace"),
                        tip.only = FALSE) {
    tree <- phy_tree(physeq)
    ## Get group factor 
    if (!is.null(tax_table(physeq, FALSE))) {
        if (class(group) == "character" & length(group) == 1) {
            x1 <- as(tax_table(physeq), "matrix")
            if (!group %in% colnames(x1)) {
                stop("group not found among sample variable names.")
            }
            group <- x1[, group]
            group[is.na(group)] <- "Unassigned"
        }
    }
    if (class(group) != "factor") {
        group <- factor(group)
    }
    ## Reorder group to match tree tip labels
    names(group) <- taxa_names(physeq)
    group <- group[tree$tip.label]
    ## Create color palette
    color.palette <- gg_color_hue(length(levels(group)))
    names(color.palette) <- levels(group)
    if ("Unassigned" %in% levels(group)) {
        color.palette <- c(color.palette[ names(color.palette) != "Unassigned" ], "grey")
        names(color.palette)[length(color.palette)] <- "Unassigned"
    }
    ## Create tip color vector from group and color palette
    tip.color <- color.palette[ group ]
    ## Get group of inner nodes
    if (tip.only) {
        edge.color <- NULL
    } else {
        method <- match.arg(method) 
        if (method == "majority") { 
            node.descendants <- prop.part(tree)
            MostPresentDescendants <- function(x) {
                descendants.groups <- table(group[x])
                return(names(descendants.groups)[which.max(descendants.groups)])
            }
            group.inner.node <- unlist(lapply(node.descendants, MostPresentDescendants))
        } else {
            cat("Using ace to estimate inner node groups, be patient...\n")
            inner.node.state <- ace(group, tree, type = "discrete", model = "ER")
            inner.node.map.state <- apply(inner.node.state$lik.anc, 1, which.max)
            group.inner.node <- levels(group)[ inner.node.map.state ]
        }
        ## Regroup inner nodes and leaves
        group <- c(as.character(group), group.inner.node)
        group <- factor(group)
        ## Create edge color vector from group and color palette
        ## group[tree$edge[i , 2]] is the group of downstream node of edge i
        edge.color <- color.palette[ group[tree$edge[ , 2]] ]
    }
    ## Return results 
    return(list(edge = edge.color, tip = tip.color, palette = color.palette))
}


#' Plot tree with edges colored by `group`. `group` should give a class for all leaves.
#' The class of inner nodes is inferred using majority rule or ancestral character
#' estimation (via ape::ace function, used with ER model)
#'
#' @title tree_colored_by
#' @param physeq Required. An instance of a \code{phyloseq} class object that has
#'               tree slot. 
#' @param group  Optional. Either a single character string matching a variable
#'               name in the corresponding tax_table of `physeq`, or a factor with
#'               the same length as the number of taxa in `physeq`. Defaults to "Phylum"
#' @param legend.title Optional. Character string used for legend title. Defaults to
#'                     `group`.
#' @param legend.pos  Optional. keyword used for legend positon. Defaults to "bottomright"
#'                    and can take any value accepted by legend. "none" removes legend. 
#' @param method      Optional. Ancestral group reconstruction method; either 'majority' or
#'                    'ace' (or an abbrevation of these). Defaults to "majority"
#' @param plot        Optional. Logical. Should the tree be plotted or not. 
#' @param ...    Optional. Additional arguments passed on to plot.phylo.
#' @note The function assumes that the tree is rooted. 
#' @return The function is mainly called for its side effect of plotting a tree
#'         with colored edges and leaves. 
tree_colored_by <- function(physeq, group = "Phylum", method = c("majority", "ace"), 
                            legend.title = deparse(substitute(group)),
                            legend.pos = "bottomright", plot = TRUE, ...) {
    ## Exception handling
    if (is.null(phy_tree(physeq, FALSE))) {
        stop("Object \"physeq\" does not have a tree slot")
    } else {
        tree <- phy_tree(physeq)
    }
    if (!is.rooted(tree)) {
        stop("Tree must be rooted, consider using midpoint rooting")
    }
    color <- color_edges(physeq, group, method)
    if (plot) {
        ## plot tree
        plot(tree, edge.color = color$edge, tip.color = color$tip, ...)
        ## Add legend
        if (legend.pos != "none") {
            legend(x = legend.pos,
                   xjust = 0,
                   legend = names(color$palette),
                   fill = color$palette,
                   border = NA,
                   title = legend.title,
                   bty = "n")
        }
    }
}


#' Subsidiary function for plot_merged_trees. Transforms a list as returned by fattenEdges
#' to a list of aesthetics parameters for plot.phylo
#'
#' @title construct_aesthetics
#' @param x Required. List as returned by fattenEdges.
#' @param color Required. Colors used for edges and tips. If NULL (default), all
#'              edges and tips will be black.
#' @param fatten.edges.by Required. Aesthetics that will be reconstructed. Subset of 'size',
#'                        'alpha' and 'color'
#' @param deviation Required. Logical. Used to construct an automatic color scale. If TRUE,
#'                  a divergent color scale is used, if FALSE, a sequential one.
#' @return A list with components
#' \item{edge.size} A matrix of size nsamples(physeq) x number of edges in \code{phy_tree(physeq)}
#'         with fattened edge widths
#' \item{edge.color} A matrix of size nsamples(physeq) x number of edges in \code{phy_tree(physeq)}
#'         with colored edge widths
#' \item{legend} A named list with components \item{size}, \item{alpha} and \item{color}
#'               used for constructing a legend.
construct_aesthetics <- function(x, color, fatten.edges.by,
                                 deviation = FALSE) {
    edge <- x$edge.width
    legend <- x$legend
    ## Start with color
    if (is.null(color)) {
        if ("color" %in% fatten.edges.by) {
            if (deviation) {
                ## Choose diverging palette (brewer.pal(11, "RdYlBu"))
                ## To do, manually create color scale with different thresholds. 
                color.palette <- rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61",
                                       "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9",
                                       "#74ADD1", "#4575B4", "#313695"))
            } else {
                ## Choose sequential palette (brewer.pal(9, "YlOrRd"))
                ## color.palette <- c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026")
                ## Choose sequential palette (brewer.pal(9, "BuPu"))
                color.palette <- c("#F7FCFD", "#E0ECF4", "#BFD3E6", "#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D", "#810F7C", "#4D004B")
            }
            color.scale <- gradient_n_pal(color.palette,
                                          seq(min(legend$breaks), max(legend$breaks),
                                              length.out = length(color.palette)))
            edge.color <- color.scale(edge)
            ## Reformat as a matrix
            edge.color <- matrix(edge.color, nrow = nrow(edge), ncol = ncol(edge))
            legend$color.scale <- color.scale
            ## Recover values outside label range
            edge.color[edge > max(legend$breaks)] <- color.scale(max(legend$breaks))
            edge.color[edge < min(legend$breaks)] <- color.scale(min(legend$breaks))
        } else {
            edge.color <- matrix("black", nrow = nrow(edge), ncol = ncol(edge))
        }
    } else {
        edge.color <- matrix(rep(color$edge, each = nrow(edge)), nrow = nrow(edge))
        legend$color.scale <- color$palette
    }
    ## Continue with alpha
    if ("alpha" %in% fatten.edges.by) {
        scaled.alpha <- rescale(edge, from = c(0, max(legend$breaks, edge)))
        edge.color <- alpha(edge.color, scaled.alpha)
        ## Reformat as matrix
        edge.color <- matrix(edge.color, nrow = nrow(edge), ncol = ncol(edge))
        legend$alpha.scale <- alpha("black",
                                    rescale(legend$breaks, from = range(legend$breaks)))
        ## If also fattening by color, udpdate color scale
        if (is.function(legend$color.scale)) {
            alpha.color.palette <- alpha(color.palette,
                                         seq(0, 1, length.out = length(color.palette)))
            color.scale <- gradient_n_pal(alpha.color.palette,
                                          seq(min(legend$breaks), max(legend$breaks),
                                              length.out = length(color.palette)))
            legend$color.scale <- color.scale
        }
    }
    ## Continue with size
    if ("size" %in% fatten.edges.by) {
        edge.size <- edge
        legend$size.scale <- legend$breaks
    } else {
        edge.size <- matrix(1, nrow = nrow(edge), ncol = ncol(edge))
    }
    return(list(edge.size  = edge.size,
                edge.color = edge.color,
                legend     = legend))
}


#' Wrapper around plot.phylo to plot tree with precomputed aesthetics.
#' Only aesthetics specified in fatten.edges.by fed to plot.phylo. The goal
#' of this function is to allow the user to feed plot.phylo with custom
#' aesthetics in addition to those automatically computed by plot_merged_trees
#' and use only the custom ones or the automatically computed ones. 
#'
#' @title plot_pretty_tree
#' @param tree Required. A \code{phylo} class object
#' @param edge.color Required. A color vector potentially passed on to plot.phylo
#' @param tip.color Required. A color vector potentially passed on to plot.phylo
#' @param edge.size Required. An edge expansion vector, potentially passed on to plot.phylo
#' @param tip.size. Required. A tip expansion vector, potentially passed on to plot.phylo
#' @param fatten.edges.by  Required. Aesthetics used for edge "fattening", vector with elements in
#'                         'size', 'alpha', 'color'. Defaults to 'size' only. Determines which
#'                         of 'edge.color', 'edge.size' are passed on to
#'                         plot.phylo
#' @param fatten.tips Required. Aesthetics used for tip "fattening", vector with elements in
#'                   'size', 'alpha', 'color'. Defaults to 'size' only. Determines which
#'                   of 'tip.color', 'tip.size' are passed on to plot.phylo
#' @param color.tip Optional. Logical. Should tips be colored even if fatten.tips is set to false. 
#' @param ... Optional. Additional arguments passed on to plot.phylo. 
#' @return Nothing. This function is used for its side effect of plotting a tree.
#' @note This function is intended for use within plot_merged_tree, no argument checking is performed. 
plot_pretty_tree <- function(tree, edge.size, edge.color,
                             tip.size, tip.color, 
                             fatten.edges.by, fatten.tips, color.tip, ...) {
    ## Prepare arguments for plot.phylo
    args <- list(x = tree)
    additional.args <- list(...)
    ## If size is a fattening factor, add edge.size and tip.size to args and
    ## remove them to additional.args
    if ("size" %in% fatten.edges.by) {
        args$edge.width <- edge.size
        additional.args$edge.width <- NULL
        if (fatten.tips) {
            args$cex <- tip.size
            additional.args$cex <- NULL
        }
    }
    if (any(c("color", "alpha") %in% fatten.edges.by)) {
        args$edge.color <- edge.color
        additional.args$edge.color <- NULL
        if (fatten.tips) {
            args$tip.color <- tip.color
            additional.args$tip.color <- NULL
        }
    }
    if (!fatten.tips & color.tip) {
        args$tip.color <- alpha(tip.color, 1)
    }
    all.args <- c(args, additional.args)
    do.call("plot.phylo", all.args)
}

#' Wrapper around legend to get legend dimension. Intended for use within
#' plot_pretty_tree
#'
#' @title get_legend_dimension
#' @param leg Required. A list a produced by 'construct_aesthetics' in its
#'               'legend' component
#' @return Dimensions of the legend. 
#' @note This function is intended for use within plot_merged_tree, no argument checking is performed. 
get_legend_dimension <- function(leg) {
    ## Legend total dimensions
    leg.dim <- list(x = 0, y = 0)
    ## Get dimensions of color legend
    color.scale <- leg$color.scale
    if (is.function(color.scale)) {
        legend.color <- legend(x = "top",
                               legend = leg$labels,
                               col  = "black",
                               title = "Frequency",
                               bty = "n",
                               plot = FALSE)
        ## Keep rectangle dimensions and
        ## Add space for the title
        color.dim <- list(rect = legend.color$rect)
        color.dim$rect$h <- color.dim$rect$h + abs(diff(legend.color$text$y[1:2]))
        ## Compute color rectangle dimentions for
        ## use in plotrix::color.legend
        color.dim$box <- legend.color$rect
        color.dim$box$w <- legend.color$text$x - legend.color$rect$left
    }
    if (is.vector(color.scale)) {
        color.dim <- legend(x = "top",
                            legend = names(color.scale),
                            fill  = color.scale,
                            title = "Group",
                            bty = "n",
                            plot = FALSE)
    }
    if (is.null(color.scale)) {
        color.dim <- list(rect = list(w = 0, h = 0))
    }
    ## Update dimensions of legend
    leg.dim$x <- max(leg.dim$x, color.dim$rect$w)
    leg.dim$y <- leg.dim$y + color.dim$rect$h * 1.1 ## expansion factor
    ## Get dimensions of size legend
    size.scale <- leg$size.scale
    alpha.scale <- leg$alpha.scale
    ## If size.scale does not exist, construct additional legend
    ## only if alpha.scale exists and color.scale is not a function
    if (is.null(size.scale)) {
        if (!is.null(alpha.scale) & !is.function(color.scale)) {
            size.scale <- rep(1, length(alpha.scale))
            size.dim <- legend("top",
                               legend = leg$labels,
                               lwd = size.scale,
                               title = "Frequency",
                               bty = "n",
                               plot = FALSE)
            leg.dim$x <- max(leg.dim$x, size.dim$rect$w)
            leg.dim$y <- leg.dim$y + size.dim$rect$h
        }
    } else {
        if (is.null(alpha.scale)) {
            alpha.scale <- rep("black", length(size.scale))
        }
        size.dim <- legend("top",
                           legend = leg$labels,
                           lwd = size.scale,
                           title = "Frequency",
                           bty = "n",
                           plot = FALSE)
        leg.dim$x <- max(leg.dim$x, size.dim$rect$w)
        leg.dim$y <- leg.dim$y + size.dim$rect$h * 1.1
    }
    return(list(leg.dim = leg.dim, color.dim = color.dim,
                alpha.scale = alpha.scale,
                color.scale = color.scale,
                size.scale = size.scale))
}

#' Wrapper around legend that constructs automatic legend for the graph.
#' Intended for use within plot_merged_trees
#'
#' @title plot_pretty_legend
#' @param leg Required. A list as produced by 'construct_aesthetics' in its
#'               'legend' component
#' @param x      Required. x-coordinate around which the legend is
#'               centered.
#' @param y      Required. y-coordinate around which the legend is
#'               centered.
#' @return Nothing. This function is used for its side effect of plotting
#'         a legend.
#' @note This function is intended for use within plot_merged_tree, no argument checking is performed. 
plot_pretty_legend <- function(x, y, leg) {
    res <- get_legend_dimension(leg)
    leg.dim <- res$leg.dim
    color.dim <- res$color.dim
    alpha.scale <- res$alpha.scale
    color.scale <- res$color.scale
    size.scale <- res$size.scale
    ## Left justify legends
    if (!is.null(color.scale)) {
        if (is.vector(color.scale)) {
            legend(x = x - leg.dim$x,
                   y = y + leg.dim$y / 2,
                   legend = names(color.scale),
                   fill  = color.scale,
                   title = "Group",
                   bty = "n",
                   xpd = NA)
        } else {
            rect.col <- color.scale(seq(min(leg$breaks),
                                        max(leg$breaks),
                                        length.out = 100))
            color.legend(x = x - leg.dim$x,
                         y = y + leg.dim$y / 2,
                         h = color.dim$rect$h,
                         col.h = color.dim$box$h,
                         col.w = color.dim$box$w,
                         rect.col  = rect.col,
                         title = "Frequency",
                         legend = leg$labels,
                         values = leg$breaks,
                         xpd = NA)
        }
    }
    if (!is.null(size.scale)) {
        legend(x = x - leg.dim$x,
               y = y + leg.dim$y / 2 - color.dim$rect$h,
               legend = leg$labels,
               lwd  = size.scale,
               col  = alpha.scale,
               title = "Frequency",
               bty = "n",
               xpd = NA)        
    }
}


#' Merge samples according to group factor and plot taxa diversity for each categrory of the group. 
#'
#' @title plot_merged_trees
#' @param physeq Required. An instance of a \code{phyloseq} class object that has
#'               tree slot. 
#' @param group  Required. Either a single character string matching a variable
#'               name in the corresponding sam_data of `physeq`, or a factor with
#'               the same length as the number of taxa in `physeq`.
#' @param freq    Logical. Should counts be replaced by frequencies (TRUE) or kept as such (FALSE).
#'                 Defaults to TRUE. TRUE corresponding to weighting all samples equally while FALSE
#'                 weights them with their library size.
#' @param method Optional. Edge fattening method; either 'linear' or
#'                    'logarithmic' (or an abbrevation of these). Defaults to "linear"
#' @param width.lim Optional. Numeric. Minimum and maximal edge after fattening.
#'                  Defaults to c(0.1, 4). Set to c(0, x) for true linear scaling.
#' @param missing.color Optional. Color used for tips and edges not present in the sample.
#'                      Defaults to "gray". If NULL, nothing happens. Use "white", "transparent",
#'                      or par("bg") to remove them from plot. 
#' @param color.edge.by      Optional. A single character string matching a variable. If NULL,
#'                           nothing happens. 
#'                           name in the corresponding tax_table of `physeq`, or a factor with
#'                           the same length as the number of taxa in `physeq`.
#'                           Defaults to NULL, corresponding to no color.
#' @param color.edge.method  Optional. Ancestral group reconstruction method; either 'majority' or
#'                           'ace' (or an abbrevation of these). Defaults to "majority"
#' @param show.missing.tip   Optional. Logical. Should tips not present in the sample be still be
#'                           labelled. If FALSE, corresponding labels are replaced with ""
#' @param fatten.edges.by    Optional. Aesthetics used for edge "fattening", vector with elements in
#'                           'size', 'alpha', 'color'. Defaults to 'size' only.
#' @param fatten.tips        Optional. Logical. Should tips labels be "fattened" as edges.
#'                           Defaults to FALSE.
#' @param base      Optional. Numeric. Base used for logarithmic scaling. Defaults to 10.
#' @param legend.title Optional. Character string used for legend title. Defaults to
#'                     `group`.
#' @param deviation Optional. Logical. Should the frequency (or count) of a taxa in a sample
#'                  be normalized with respect to the average of that taxa over all samples
#'                  before fattening. Defaults to FALSE. If method = "linear", the values are
#'                  centered around the mean (deviation). If method = "logarithmic", they are
#'                  divided by the mean (relative deviation)
#' @param ...    Optional. Additional arguments passed on to plot.phylo.
#' @note The function assumes that the tree is rooted. 
#' @return The function is mainly called for its side effect of plotting a tree
#'         with colored edges and leaves. 
plot_merged_trees <- function(physeq, group, freq = TRUE, method = "linear",
                              missing.color = "gray",
                              show.missing.tip = TRUE,
                              fatten.edges.by = c("size"),
                              fatten.tips = FALSE,
                              color.edge.by = NULL, color.edge.method = "majority",
                              base = 10, width.lim = c(0.1, 4), deviation = FALSE, 
                              legend.title = deparse(substitute(group)), ...) {
    ## Exception handling
    if (is.null(phy_tree(physeq, FALSE))) {
        stop("Object \"physeq\" does not have a tree slot")
    } else {
        tree <- phy_tree(physeq)
    }
    if (!is.rooted(tree)) {
        stop("Tree must be rooted, consider using midpoint rooting")
    }

    args <- list(...)
    ## Check consistency
    if (any(! fatten.edges.by %in% c("size", "color", "alpha"))) {
        stop("Elements of fatten.edges.by must belong to size, color and alpha")
    }
    if (fatten.tips & "size" %in% fatten.edges.by & "cex" %in% names(args)) {
        stop("Can't use cex when fattening edges and tips by size")
    }
    if ("color" %in% fatten.edges.by & !is.null(color.edge.by)) {
        stop("Can't use color.edge.by when fattening edges by color")
    }
    if (any(c("size", "alpha") %in% fatten.edges.by) & deviation) {
        warning("Fattening edges by size or alpha when\n using deviation may not be very relevant")
    }
    
    ## Potentially normalize counts to frequencies
    if (freq) {
        physeq <- transform_sample_counts(physeq, function(x) {x /sum(x)} )
    }
    
    ## Merge samples by grouping factor, if NULL merge all samples together
    if (is.null(group)) {
        group <- factor(rep("All", nsamples(physeq)))
    }
    physeq <- merge_samples(physeq, group)
    
    ## Fatten edges for future use in aesthetics
    fattened.edges <- fattenEdges(physeq, method, width.lim, base, deviation)
    ## Color edges and add color to fatten.edges.by for plot_pretty_tree
    if (!is.null(color.edge.by)) {
        color <- color_edges(physeq, color.edge.by, color.edge.method)
        fatten.edges.by <- unique(c(fatten.edges.by, "color"))
    } else {
        color <- NULL
    }

    ## Construct edge and tip aesthetics
    aesthetics <- construct_aesthetics(fattened.edges, color, fatten.edges.by, deviation)
    
    ## Optionally, update missing edges and tips color and add color to
    ## fatten.edges.by for plot_pretty_tree
    if (!is.null(missing.color)) {
        aesthetics$edge.color[ !fattened.edges$edge.presence ] <- missing.color
        fatten.edges.by <- unique(c(fatten.edges.by, "color"))
    }

    ## Prepare layout for plot
    number.categories <- nsamples(physeq)
    layout.nrow <- floor(sqrt(number.categories))
    layout.ncol <- ceiling(number.categories / layout.nrow)

    par(omd = c(0, 0.9, 0, 1), mfcol = c(layout.nrow, layout.ncol), mar = c(0,0,1,0))
    ## Plot trees
    for (i in 1:number.categories) {
        ## Get edge and tip colors
        sample.edge.color <- aesthetics$edge.color[i, ]
        sample.tip.color <- sample.edge.color[fattened.edges$pendant.edges]
        ## Get edge and tip size
        sample.edge.size <- aesthetics$edge.size[i, ]
        sample.tip.size <- sample.edge.size[fattened.edges$pendant.edges]
        
        ## prepare tree labels        
        tree <- phy_tree(physeq)
        if (!show.missing.tip) {
            absent.tips <- which(! fattened.edges$tip.presence[i, ])
            tree$tip.label[ absent.tips ] <- ""
        }
        
        ## Plot tree according to aesthetics
        plot_pretty_tree(tree, sample.edge.size, sample.edge.color,
                         sample.tip.size, sample.tip.color,
                         fatten.edges.by, fatten.tips,
                         color.tip = !is.null(color.edge.by), ...)
        title(sample_names(physeq)[i])
    }

    ## Add common legends in outer margins for whichever of alpha, size and color requires it
    ## Convert normalized device coordinates
    ## in current plot ("user") coordinates
    legend.x = grconvertX(1, "ndc", "user")
    legend.y = grconvertY(0.5, "ndc", "user")
    plot_pretty_legend(x = legend.x, y = legend.y, leg = aesthetics$legend)
}



#' These functions are adapted from plotrix package and simplified here
#' for our purpose
color.legend <- function(x, y, h, col.h, col.w, title, rect.col,
                         legend, values = NULL, 
                         xpd = NA, ...) {
    if (!missing(xpd)) {
        op <- par("xpd")
        on.exit(par(xpd = op))
        par(xpd = xpd)
    }
    ## insert title
    text(x = x + col.w,
         y = y - (h - col.h) / 4,
         labels = title)
    xl <- x
    xr <- x + col.w - 0.2 * strwidth("O")
    yt <- y - (h - col.h)
    yb <- yt - col.h
    gradient.rect(xl, yb, xr, yt, col = rect.col, nslices = length(rect.col))
    ysqueeze <- (yt - yb)/(2 * length(rect.col))
    if (is.null(values)) {
        texty <- seq(yb + ysqueeze, yt - ysqueeze, length.out = length(legend))
    } else {
        f <- approxfun(range(values), c(yb + ysqueeze, yt - ysqueeze))
        texty <- f(values)
    }
    textx <- xr + 0.2 * strwidth("O")
    textadj <- c(0, 0.5)
    text(textx, texty, labels = legend, adj = textadj, ...)
}

gradient.rect <- function (xleft, ybottom, xright, ytop, col = NULL,
                          nslices = 50, border = par("fg")) {
    yinc <- (ytop - ybottom)/nslices
    ybottoms <- seq(ybottom, ytop - yinc, length = nslices)
    ytops <- ybottoms + yinc
    rect(xleft, ybottoms, xright, ytops, col = col, lty = 0)
    rect(xleft, ybottoms[1], xright, ytops[nslices], 
         border = border)
}

symmetric_log_breaks <- function (n = 5, base = 10) {
    function(x) {
        rng <- max(abs(log(range(x, na.rm = TRUE), base = base)))
        rng <- c(-rng, rng)
        min <- floor(rng[1])
        max <- ceiling(rng[2])
        if (max == min) 
            return(base^min)
        by <- floor((max - min)/n) + 1
        breaks <- base^seq(min, max, by = by)
        unique(sort(c(1, breaks)))
    }
}
