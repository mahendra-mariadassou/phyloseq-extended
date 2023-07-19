## Modification of phyloseq import_biom function to format Function to read the biom file produced by FROGS and load it as a phyloseq object

#' Import function to read FROGS format OTU table
#'
#' @inheritParams phyloseq::import_biom
#' @inheritDotParams phyloseq::import_biom
#' @param taxMethod Default "blast". Either blast or RDP, the method used for taxonomy affiliation.
#'
#' @return A \code{\link{phyloseq-class}} object
#' @export
#'
#' @importFrom biomformat biom_data read_biom sample_metadata
#' @importFrom methods as
#' @importFrom phyloseq otu_table read_tree sample_data tax_table
#' @importFrom purrr map map_chr map_int map_lgl pluck
#'
#' @examples
#' biom_file <- system.file("extdata", "frogs_data.biom", package = "phyloseq.extended")
#' data <- import_frogs(biom_file)
#' data
import_frogs <- function(biom,
                         treefilename = NULL,
                         refseqfilename = NULL,
                         refseqFunction = readDNAStringSet,
                         refseqArgs = NULL,
                         taxMethod = c("blast", "rdp"),
                         parallel = TRUE, ...) {
    x <- biomformat::read_biom(biom)
    ## otu table
    otuTable <- otu_table(as(biom_data(x), "matrix"), taxa_are_rows = TRUE)
    ## sample data
    if (is.null(sample_metadata(x))) {
      samdata <- NULL
    }
    else {
      samdata = sample_data(sample_metadata(x))
      # failsafe of biom_data() when there is only one sample
      if (ncol(otuTable) == 1 & nrow(samdata) == 1) {
        colnames(otuTable) <- rownames(samdata)
      }
    }
    argumentlist <- c(list(otuTable), list(samdata))
    ## tax table
    taxMethod <- match.arg(taxMethod)
    field <- paste(taxMethod, "taxonomy", sep = "_")
    taxnames <- purrr::map_chr(x$rows, "id")
    taxlist <- purrr::map(x$rows, purrr::pluck, "metadata", field)
    if (all(purrr::map_lgl(taxlist, is.null))) {
        taxtab <- NULL
        warning(paste("No taxonomy tag", field ,"was found in the biom file"))
    } else {
        n_taxa <- length(taxnames)
        n_ranks <- purrr::map_int(taxlist, length)
        max_n_ranks <- max(n_ranks)
        ranknames <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
                       paste0("Rank_", seq(from = 8, length.out = max(0, max_n_ranks - 7))))[1:max_n_ranks]
        ## Empty tax table
        taxtab <- matrix(NA_character_, nrow = n_taxa, ncol = max_n_ranks,
                         dimnames = list(taxnames, ranknames))
        ## Fill tax table using two column index array
        ii <- cbind(rep(1:n_taxa, times = n_ranks),               ## row index
                    purrr::map(taxlist, seq_along) %>% unlist()   ## col index
                    )
        taxtab[ii] <- unlist(taxlist)
    }
    # if (all(sapply(sapply(x$rows, function(i) {
    #     i$metadata[[field]]
    # }), is.null))) {
    #     taxtab <- NULL
    #     warning(paste("No taxonomy tag", field ,"was found in the biom file"))
    # } else {
    #     taxlist <- lapply(x$rows, function(i) {
    #         tmp <- i$metadata[[field]]
    #         if (is.null(tmp)) {
    #             return(matrix(NA, nrow = 1))
    #         } else {
    #             return(matrix(tmp, nrow = 1))
    #         }
    #     })
    #     taxnames <- vapply(x$rows, function(i) {
    #         i$id
    #     }, "character")
    #     taxtab <- plyr::rbind.fill.matrix(taxlist)
    #     rownames(taxtab) <- taxnames
    #     colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:ncol(taxtab)]
    #     taxtab <- tax_table(taxtab)
    # }
    argumentlist <- c(argumentlist, list(tax_table(taxtab)))
    ## tree
    if (!is.null(treefilename)) {
        if (inherits(treefilename, "phylo")) {
            tree = treefilename
        }
        else {
            tree <- read_tree(treefilename, ...)
        }
        if (is.null(tree)) {
            warning("treefilename failed import. It not included.")
        }
        else {
            argumentlist <- c(argumentlist, list(tree))
        }
    }
    ## reference sequences
    if (!is.null(refseqfilename)) {
        if (inherits(refseqfilename, "XStringSet")) {
            refseq = refseqfilename
        }
        else {
            if (!is.null(refseqArgs)) {
                refseq = do.call("refseqFunction", c(list(refseqfilename),
                  refseqArgs))
            }
            else {
                refseq = refseqFunction(refseqfilename)
            }
        }
        argumentlist <- c(argumentlist, list(refseq))
    }
    return(do.call("phyloseq", argumentlist))
}
