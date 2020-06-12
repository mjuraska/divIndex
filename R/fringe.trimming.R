#' Fringe Trimming
#'
#' This function applies the fringe-trimming bias correction method (Gregori et al. 2014)
#' to sequences and read counts in an individual's sample, given confidence and filtering levels, and returns
#' a list containing the fringe-trimmed sequences and read counts. Fringe-trimming reduces potential
#' bias due to difference in the number of sequence reads across samples. Haplotypes may be
#' nucleotide or amino acid sequences.
#'
#' @param seqs A character vector containing the distinct haplotype sequences from an individual subject.
#' @param count A numeric vector of observed read counts for each haplotype.
#' @param conf A numeric value indicating the confidence level for the trimming. The default value is 0.9.
#' @param filter A numeric value indicating level of filtering. The default value is 0.01.
#'
#' @details Sample estimates of diversity and richness tend to be increasing functions
#' of sample size, with larger samples yielding higher observed diversity. This leads to bias when
#' sample read counts vary drastically among individuals. Traditional methods of correcting for this
#' sample size bias include the Chao-1 and Chao-2 estimators, the abundance-based coverage estimator,
#' and the incidence-based coverage estimator (Colwell and Coddington, 1994; Gotelli and Caho, 2013).
#' However, these estimators are calculated using observed singletons and doubletons in samples,
#' which represent rare haplotypes.
#' In data sequenced using NGS methods, singletons and doubletons are eliminated during a filtration
#' process that filters out haplotypes below a noise threshold level in an effort to eliminate artificially
#' created haplotypes arising from sequencing errors. As a result, traditional metrics for correcting
#' sample size bias are not available. This function applies a fringe-trimming bias correction proposed by
#' Gregori et al. (2014) to correct sample size bias in NGS filtered data.
#'
#' For a given sample's sequences and read counts, we exclude sequence reads representing haplotypes for which
#' \deqn{P(N \le n_i | n, p=0.01) < 0.9,} where \eqn{N} is a binomial random variable with
#' parameters \eqn{n}, the total number of sequence reads for an individual sample, {p},
#' the true haplotype frequency, and \eqn{n_i}, the number of sequence reads for the \eqn{i}-th
#' haplotype. The default filtration cut-off is 1\% and the default confidence level is 90\%,
#' based on analyses by Gregori et al. (2014). Thus, haplotypes for which the probability of
#' observing as few as \eqn{n_i} reads is less than 90\%, given a total read count of \eqn{n} and true abundance
#' of 1\%, will be excluded.
#'
#' @return Returns a list containing the modified \code{seqs} and \code{count} vectors with
#' fringe haplotypes removed.
#'
#' @references
#' Colwell, R. K., & Coddington, J. A. (1994). Estimating terrestrial biodiversity through extrapolation.
#' \emph{Phil. Trans. R. Soc. Lond. B, 345}(1311), 101-118.
#'
#' Gotelli, N. J., & Chao, A. (2013). Measuring and estimating species richness, species diversity,
#' and biotic similarity from sampling data.
#'
#' Gregori, J., Salicru, M., Domingo, E., Sanchez, A., Esteban, J. I., Rodriguez-Frias, F., & Quer, J. (2014).
#' Inference with viral quasispecies diversity indices: clonal and NGS approaches. \emph{Bioinformatics, 30}(8),
#' 1104-1111.
#'
#' @examples
#' seqs <- c("ACDEFGHI","KLMNPQRS","TVWYACDE","SMWLWCAN")
#' count <- c(80,200,366,14)
#' conf <- 0.9
#' filter <- 0.1
#' fringe.trimming(seqs, count, conf, filter)
#'
#' @importFrom stats qbinom
#'
#' @export
fringe.trimming <- function(seqs, count, conf = 0.9, filter = 0.01){
  cutoff <- qbinom(conf, size=sum(count), prob = filter, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seqs[j] <- NA
      count[j] <- NA
    }
  }
  seqs <- seqs[!is.na(seqs)]
  count <- count[!is.na(count)]
  list(seqs, count)
}

