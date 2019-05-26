## Modification of phyloseq import_biom function to format Function to read the biom file produced by FROGS and load it as a phyloseq object

#' Import function to read FROGS format OTU table
#'
#' @param biom A character string indicating the file location of the BIOM formatted file. This is a JSON formatted file, specific to biological datasets, as described in http://www.qiime.org/svn_documentation/documentation/biom_format.htmlthe biom-format home page.
#' @param treefilename Default value is NULL. A file representing a phylogenetic tree or a phylo object. Files can be NEXUS or Newick format.
#' @param refseqfilename Default NULL. The file path of the biological sequence file that contains at a minimum a sequence for each OTU in the dataset. Alternatively, you may provide an already-imported XStringSet object that satisfies this condition. In either case, the names of each OTU need to match exactly the taxa_names of the other components of your data.
#' @param refseqFunction Default is readDNAStringSet, which expects to read a fasta-formatted DNA sequence file. If your reference sequences for each OTU are amino acid, RNA, or something else, then you will need to specify a different function here. This is the function used to read the file connection provided as the the previous argument, refseqfilename. This argument is ignored if refseqfilename is already a XStringSet class.
#' @param refseqArgs Default NULL. Additional arguments to refseqFunction. See XStringSet-io for details about additional arguments to the standard read functions in the Biostrings package.
#' @param taxMethod Default "blast". Either blast or RDP, the method used for taxonomy affiliation.
#' @param parallel Logical. Wrapper option for .parallel parameter in plyr-package functions. If TRUE, apply parsing functions in parallel, using parallel backend provided by foreach and its supporting backend packages. One caveat, plyr-parallelization currently works most-cleanly with multicore-like backends (Mac OS X, Unix?), and may throw warnings for SNOW-like backends. See the example below for code invoking multicore-style backend within the doParallel package.
#' @param ... Additional parameters passed on to read_tree.
#'
#' @return A phyloseq class object
#' @export
#'
#' @importFrom biomformat read_biom biom_data sample_metadata
#' @importFrom plyr rbind.fill.matrix
import_frogs <- function(biom,
                         treefilename = NULL,
                         refseqfilename = NULL,
                         refseqFunction = readDNAStringSet,
                         refseqArgs = NULL,
                         taxMethod = c("blast", "rdp"),
                         parallel = TRUE, ...) {
    argumentlist <- list()
    x <- biomformat::read_biom(biom)
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
