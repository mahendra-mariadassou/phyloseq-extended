estimate_prevalence <- function(physeq, group, format = c("long", "wide")) {
    ## Args:
    ## - physeq: phyloseq class object, otu abundances are extracted from this object
    ## - group:  Either the a single character string matching a
    ##           variable name in the corresponding sample_data of ‘physeq’, or a
    ##           factor with the same length as the number of samples in ‘physeq’.
    ## - format: long (default) or wide. Should the result be organised as a long data.frame
    ##           with columns "otu", "group", "prevalence" and "specificity" or a numeric matrix
    ##           with otus in rows and group specificty/group prevalence in column.
    ## 
    ## Returns;
    ## data frame with components
    ## - prevalence: observed prevalence in group
    ## - specifity : specifity (prevalence / global prevalence) in group
    ## - group: factor level
    ## - otu: otu
    ## or
    ## matrix m of size ntaxa(physeq) times 2*N where N is the number of levels (groups) in group
    ## and
    ## - m[i, g] is the prevalence of otu i in group g
    ## - m[i, N+g] is the specificity of otu i to group g
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
    tdf <- 0 + (tdf > 0)
    ## prevalence
    frac <- t(rowsum(t(tdf), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(tdf)),
                                                              nrow = nrow(tdf))
    ## specificity
    spec <- t(rowsum(t(tdf), group, reorder = TRUE)) / rowSums(tdf)
    spec[rowSums(tdf) == 0, ] <- 0 
    ## cbind prevalence and specificity
    res <- merge(melt(frac, varnames = c("otu", "group"), value.name = "prevalence"),
                 melt(spec, varnames = c("otu", "group"), value.name = "specificity"))
    recast_prevalence <- function() {
        prev <- acast(prevalence, otu ~ group, value.var = "prevalence")
        colnames(prev) <- paste("prev", colnames(prev), sep = "_")
        spec <- acast(prevalence, otu ~ group, value.var = "specificity")
        colnames(prev) <- paste("spec", colnames(spec), sep = "_")
        return(cbind(prev, spec))
    }
    format <- match.arg(format)
    if (format == "wide") {
        res <- recast_prevalence()
    }
    return(res)
}


