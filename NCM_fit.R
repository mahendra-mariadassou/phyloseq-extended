## fitting function for neutral community model of Sloan et al. (2006).
## Assumes equivalent species, migration rate m (to fit) and an array of relative abundances.

## likelihood of the NCM using for a given migration rate and threshold. The local abundances are estimated from the data using a BLUE. 
ncm_ll <- function(abundances, migration = 0.1, threshold = 0.01, absolute.threshold = NULL) {
    Nt <- colSums(abundances)
    local.abundances <- local_abundances(abundances) ## pi
    ## compute the probability of appearance at a given threshold of an OTU given Nt, pi and the threshold
    .local <- function(pi, N, lower.tail = FALSE) {
        ## Parameters of the beta distributions
        alpha <-  migration*N*pi
        beta <-  migration*N*(1-pi)
        if (is.null(absolute.threshold)) {
            return(pbeta(threshold, shape1 = alpha, shape2 = beta,
                         lower.tail = lower.tail, log.p = TRUE))
        } else {
              return(pbeta(absolute.threshold/N, shape1 = alpha, shape2 = beta,
                           lower.tail = lower.tail, log.p = TRUE))
          }
    }
    ## Probabilities of appearance
    lprobs.presence <- outer(local.abundances, Nt, .local)
    lprobs.absence <- outer(local.abundances, Nt, .local, lower.tail = TRUE)
    llikelihood <- sum(ifelse(abundances > 0, lprobs.presence, lprobs.absence))
    return(llikelihood)
}


ncm_fit <- function(abundances, threshold, absolute.threshold = NULL) {
    optim.function <- function(m) { ncm_ll(abundances, m, threshold, absolute.threshold) }
    solution <- optimize(optim.function, interval = c(1e-10, 1), maximum = TRUE)
    return(solution)
}

ncm_curve <- function(threshold, Nm) {
    curve <- function(pi) {
        pbeta(threshold, shape1 = Nm*pi, shape2 = Nm*(1-pi), lower.tail = FALSE)
    }
    return(curve)
}

## Multi-sample ncm curve (for samples with different library sizes)
ncm_curve_multi <- function(threshold, N, m) {
    curve <- function(pi) {
        proba_sample <- rep(NA, length(N))
        for (i in seq_along(N)) {
            proba_sample[i] <- pbeta(threshold,
                                     shape1 = N[i]*m*pi,
                                     shape2 = N[i]*m*(1-pi),
                                     lower.tail = FALSE)
        }
        return(mean(proba_sample))
    }
    return(curve)
}


local_abundances <- function(abundances, migration = 0.1) {
    Nt <- colSums(abundances)
    weights <- (Nt * migration + 1) 
    weights <- weights / ( Nt * sum(weights) )
    ## estimates of local abundances, using BLUE estimates of the proportions
    weights.matrix <- matrix( rep(weights, each = nrow(abundances)) , ncol = ncol (abundances))
    local.abundances <- rowSums( abundances * weights.matrix ) ## pi
    return(local.abundances)
}
