## Modification of phyloseq import_biom function to format Function to read the biom file produced by FROGS and load it as a phyloseq object
library(biomformat)

import_frogs <- function(biom, treefilename = NULL, refseqfilename = NULL, refseqFunction = readDNAStringSet, refseqArgs = NULL, taxMethod = c("blast", "rdp"), parallel = TRUE, ...) {
    ## Args
    ## biom: (Required). a character indicating the location of the biom file
    ## taxMethod: (Optional). Default "rdp". Either 'rdp' or 'blast' (or any unambiguous abrevation of those). Taxonomic affiliation to be imported in the phyloseq object.
    ## For all other arguments, refer to the help of import_biom
    argumentlist <- list()
    x <- read_biom(biom)
    ## otu table
    otuTable <- otu_table(as(biom_data(x), "matrix"), taxa_are_rows = TRUE)
    argumentlist <- c(argumentlist, list(otuTable))
    ## tax table
    taxMethod <- match.arg(taxMethod)
    field <- paste(taxMethod, "taxonomy", sep = "_")
    if (all(sapply(sapply(x$rows, function(i) {
        i$metadata[[field]]
    }), is.null))) {
        taxtab <- NULL
        warning(paste("No taxonomy tag", field ,"was found in the biom file"))
    } else {
        taxlist <- lapply(x$rows, function(i) {
            tmp <- i$metadata[[field]]
            if (is.null(tmp)) {
                return(matrix(NA, nrow = 1))
            } else {
                return(matrix(tmp, nrow = 1))
            }
        })
        taxnames <- vapply(x$rows, function(i) {
            i$id
        }, "character")
        taxtab <- plyr::rbind.fill.matrix(taxlist)
        rownames(taxtab) <- taxnames
        colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:ncol(taxtab)]
        taxtab <- tax_table(taxtab)
    }
    argumentlist <- c(argumentlist, list(taxtab))
    ## sample data
    if (is.null(sample_metadata(x))) {
        samdata <- NULL
    }
    else {
        samdata = sample_data(sample_metadata(x))
    }
    argumentlist <- c(argumentlist, list(samdata))
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
