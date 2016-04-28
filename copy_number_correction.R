#' Use information from Vetroský and Baldrian (2013) to correct for 16S rRNA copy number variation. 
#' Given a taxonomic level (Phylum or Class within Proteobacteria), returns the average copy number
#' (CN) within that level in the database constructed by the authors. 
#' 
#' @title sample_copy_number_internal
#' @param level Required. Character. Taxonomic level for which the CN is required.
#' @param dispersion. Optional. Logical. Should the returned CN be the average within
#'                    the level or sampled according to the distribution within that level.
#'                    Distribution is approximated by a normal with mean and sd set to those
#'                    observed in the sample database. Defaults to FALSE.
#' @note The copy number is truncated to the [1, 15] range observed in the database. Levels
#'       not observed in the database have CN set to 4.2 or sample from \code{rnorm(1, 4.2, 2.7)}
#'       which correspond to the global mean and sd in the database.
#' @return A numeric value corresponding to imputed CN. 
#' @references Vetroský, T. and Baldrian, P.(2013). The variability of the 16S rRNA Gene in Bacterial
#'             Genomes and Its Consequences for Bacterial Community Analyses. _Plos One_ *8*(2):e57923
#' @seealso \code{\link{sample_copy_number}}
sample_copy_number_internal <- function(level, dispersion = FALSE) {
    known.levels <- c("Acidobacteria", "Actinobacteria", "Aquificae", "Bacteroidetes", "Caldiserica",
                      "Chlamydiae", "Chlorobi", "Chloroflexi", "Cyanobacteria", "Defferibacteres",
                      "Deinoccocus-Thermus", "Dyctioglomi", "Elusimicrobia", "Fibrobacteres",
                      "Firmicutes", "Fusobacteria", "Gemmatimonadetes", "Ignavibacteria",
                      "Nitrospirae", "Planctomycetes", "Alphaproteobacteria", "Betaproteobacteria",
                      "Deltaproteobacteria", "Epsilonproteobacteria", "Gammaproteobacteria",
                      "Spirochaetes", "Synergistetes", "Tenericutes", "Thermodelsufobacteria",
                      "Thermotogae", "Verrucomicrobia")
    if (level %in% known.levels) {
        if (dispersion) {
            copy.number <- switch(level,
                                  "Acidobacteria" = rnorm(1, 1, 0),
                                  "Actinobacteria" = rnorm(1, 3.1, 1.7),
                                  "Aquificae" = rnorm(1, 2, 0.6),
                                  "Bacteroidetes" = rnorm(1, 3.5, 1.5),
                                  "Caldiserica" = rnorm(1, 1, 0),
                                  "Chlamydiae" = rnorm(1, 1.4, 0.5),
                                  "Chlorobi" = rnorm(1, 1.7, 0.7),
                                  "Chloroflexi" = rnorm(1, 2.2, 1.2),
                                  "Cyanobacteria" = rnorm(1, 2.3, 1.2),
                                  "Defferibacteres" = rnorm(1, 2, 0),
                                  "Deinoccocus-Thermus" = rnorm(1, 2.7, 1),
                                  "Dyctioglomi" = rnorm(1, 2, 0),
                                  "Elusimicrobia" = rnorm(1, 1, 0),
                                  "Fibrobacteres" = rnorm(1, 3, 0),
                                  "Firmicutes" = rnorm(1, 5.8, 2.8),
                                  "Fusobacteria" = rnorm(1, 5, 0.7),
                                  "Gemmatimonadetes" = rnorm(1, 1, 0),
                                  "Ignavibacteria" = rnorm(1, 1, 0),
                                  "Nitrospirae" = rnorm(1, 2, 1.4),
                                  "Planctomycetes" = rnorm(1, 1.7, 0.8),
                                  "Alphaproteobacteria" = rnorm(1, 2.2, 1.3),
                                  "Betaproteobacteria" = rnorm(1, 3.3, 1.6),
                                  "Deltaproteobacteria" = rnorm(1, 2.7, 1.4),
                                  "Epsilonproteobacteria" = rnorm(1, 3, 1.1),
                                  "Gammaproteobacteria" = rnorm(1, 5.8, 2.8),
                                  "Spirochaetes" = rnorm(1, 2.4, 1),
                                  "Synergistetes" = rnorm(1, 2.5, 1),
                                  "Tenericutes" = rnorm(1, 1.6, 0.5),
                                  "Thermodelsufobacteria" = rnorm(1, 2, 0),
                                  "Thermotogae" = rnorm(1, 1.8, 1),
                                  "Verrucomicrobia" = rnorm(1, 1.8, 1)
                                  )
        } else {
            copy.number <- switch(level,
                                  "Acidobacteria" = 1,
                                  "Actinobacteria" = 3.1,
                                  "Aquificae" = 2,
                                  "Bacteroidetes" = 3.5,
                                  "Caldiserica" = 1,
                                  "Chlamydiae" = 1.4,
                                  "Chlorobi" = 1.7,
                                  "Chloroflexi" = 2.2,
                                  "Cyanobacteria" = 2.3,
                                  "Defferibacteres" = 2,
                                  "Deinoccocus-Thermus" = 2.7,
                                  "Dyctioglomi" = 2,
                                  "Elusimicrobia" = 1,
                                  "Fibrobacteres" = 3,
                                  "Firmicutes" = 5.8,
                                  "Fusobacteria" = 5,
                                  "Gemmatimonadetes" = 1,
                                  "Ignavibacteria" = 1,
                                  "Nitrospirae" = 2,
                                  "Planctomycetes" = 1.7,
                                  "Alphaproteobacteria" = 2.2,
                                  "Betaproteobacteria" = 3.3,
                                  "Deltaproteobacteria" = 2.7,
                                  "Epsilonproteobacteria" = 3,
                                  "Gammaproteobacteria" = 5.8,
                                  "Spirochaetes" = 2.4,
                                  "Synergistetes" = 2.5,
                                  "Tenericutes" = 1.6,
                                  "Thermodelsufobacteria" = 2,
                                  "Thermotogae" = 1.8,
                                  "Verrucomicrobia" = 1.8
                                  )

        }
    } else {
        copy.number <- ifelse(dispersion, rnorm(1, 4.2, 2.7), 4.2)
    }
    copy.number <- min(15, copy.number)
    copy.number <- max(1, copy.number)
    return(copy.number)
}

#' Use information from Vetroský and Baldrian (2013) to correct for 16S rRNA copy number variation. 
#' Given a vector of taxonomic level (Phylum or Class within Proteobacteria), returns the average copy
#' number (CN) for each level. The average are taken from the database reconstructed by the authors. 
#' 
#' @title sample_copy_number
#' @param level Required. Character. Taxonomic level for which the CN is required.
#' @param dispersion. Optional. Logical. Should the returned CN be the average within
#'                    the level or sampled according to the distribution within that level.
#'                    Distribution is approximated by a normal with mean and sd set to those
#'                    observed in the sample database. Defaults to FALSE. 
#' @note The copy number is truncated to the [1, 15] range observed in the database. Levels
#'       not observed in the database have CN set to 4.2 or sample from \code{rnorm(1, 4.2, 2.7)}
#'       which correspond to the global mean and sd in the database.
#' @return A numeric vector corresponding to imputed CN. 
#' @references Vetroský, T. and Baldrian, P.(2013). The variability of the 16S rRNA Gene in Bacterial
#'             Genomes and Its Consequences for Bacterial Community Analyses. _Plos One_ *8*(2):e57923
#' @seealso \code{\link{sample_copy_number_internal}}, \code{\link{extract_level}}
sample_copy_number <- Vectorize(sample_copy_number_internal, "level")


#' Extract taxonomic level from a \code{\link{phyloseq}} class object for use in
#' \code{\link{sample_copy_number}}. The required levels are "Phylum" or "Class" within
#' Proteobacteria
#'
#' @title extract_level
#' @param physeq Required. \code{\link{phyloseq}} class object with \code{tax_table} slot.
#' @return A level vector, for use in \code{\link{sample_copy_number}}
#' @seealso \code{\link{sample_copy_number}}
extract_level <- function(physeq) {
    res <- tax_table(physeq)[, "Phylum"]
    res[ res %in% "Proteobacteria" ] <- tax_table(physeq)[res %in% "Proteobacteria", "Class"]
    return(res)
}

#' Scale sample counts from a \code{\link{phyloseq}} class object according to a
#' scaling vector or using the Phylum (or Class within Proteobacteria) average copy number.
#' Average copy number are taken from Vetroský and Baldrian (2013).
#'
#' @title scale_sample_counts
#' @param physeq Required. \code{\link{phyloseq}} class object
#' @param scaling Optional. Either a numeric scaling factor or any
#'                of "Phylum" or "Abundance". If "Phylum", The scaling is reconstructed
#'                as the average copy number within "Phylum" or "Class" for Proteobacteria
#'                using \code{\link{sample_copy_number}}. If "Abundance", counts are moderated
#'                within each sample to shrink them towards evenness. See Details. 
#' @param dispersion Optional. Logical. Should the returned CN be the average within
#'                   the level or sampled according to the distribution within that level.
#'                   Distribution is approximated by a normal with mean and sd set to those
#'                   observed in the sample database. Defaults to FALSE. 
#' @return a \code{\link{phyloseq}} class objet with scaled \code{otu_table} slot.
#' @note The copy number inferred under "Phylum" scaling is truncated to the [1, 15] range observed
#'       in the database. Levels not observed in the database have CN set to 4.2 or sample from
#'       \code{rnorm(1, 4.2, 2.7)} which correspond to the global mean and sd in the database.
#'       Under "Abundance" scaling, counts are moderated within each sample by assuming that the most
#'       abundant taxa has CN of 15, the least abundant has CN of 1 and linear interpolation between
#'       these extremes. This effectively evens bacterial counts (as opposed to 16S rRNA gene counts). 
#' @references Vetroský, T. and Baldrian, P.(2013). The variability of the 16S rRNA Gene in Bacterial
#'             Genomes and Its Consequences for Bacterial Community Analyses. _Plos One_ *8*(2):e57923
#' @seealso \code{\link{sample_copy_number}}, \code{\link{extract_level}}
scale_sample_counts <- function(physeq, scaling = "Phylum", dispersion = FALSE) {
    ## Test scaling and reconstruct it if needed
    if (length(scaling) == 1 && scaling == "Abundance") {
        scaling_function <- function(x) {
            index <- (x > 0)
            x[index] <- x[index] / rescale(x[index], to = c(1, 15))
            return(x)
        }
    } else {
        if (length(scaling) == 1 && scaling == "Phylum") {
            scaling <- sample_copy_number(extract_level(physeq), dispersion = dispersion)
        } else {
            stopifnot(is.numeric(scaling), length(scaling) == ntaxa(physeq))
        }
        scaling_function <- function(x, x.scaling = scaling) {
            return(x / x.scaling[1:length(x)])
        }
    }
    return(transform_sample_counts(physeq, scaling_function))
}
