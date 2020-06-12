#' Corrected Number of Segregating Sites
#'
#' This function calculates the corrected number of segregating sites, S/a1, given a vector of an individual's
#' haplotype sequences and a vector of the observed counts for each haplotype. Segregating sites are
#' also known as polymorphic sites.
#'
#' @param seqs A character vector containing the distinct haplotype sequences from an individual subject.
#' @param count A numeric vector of observed read counts for each haplotype.
#' @param amino A logical value indicating whether haplotypes are amino acid sequences or nucleotide
#' sequences. If \code{FALSE}, nucleotide sequences are assumed.
#'
#' @details The corrected number of segregating sites is a genetic diversity index that estimates the total
#' number of polymorphic sites among the haplotypes within an individual. Also known as Watterson's
#' estimator (Watterson 1975), the corrected number of segregating sites is a measure of mutation rate
#' and is given by the equation \eqn{S/a1,} where \eqn{S} is the total number of segregating sites and
#' \eqn{a1= \sum\lim_{i=1}^{n-1} 1/i} is the correction factor, with \eqn{n} being the number of read counts.
#' Because larger numbers of reads will increase the number of segregating sites, the correction factor,
#' \eqn{a1}, is used to correct for this sample size bias that may arise when read counts vary drastically
#' between individuals.
#'
#' The corrected number of segregating sites is often used to calculate Tajima's D (Tajima 1989),
#' a population genetic test statistic used to detect selection in the evolution of a DNA sequence.
#'
#' @return Returns the corrected number of segregating sites as a numeric value.
#'
#' @references
#' Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism.
#' \emph{Genetics}, 123(3), 585-595.
#'
#' Watterson, G. A. (1975). On the number of segregating sites in genetical models without recombination.
#' \emph{Theoretical population biology}, 7(2), 256-276.
#'
#' @examples
#' # amino acid sequences
#' seqs <- c("ACDEFGHI","KLMNPQRS","TVWYACDE","SMWLWCAN")
#' count <- c(10,200,366,14)
#' amino <- TRUE
#' cor.seg.sites(seqs, count, amino)
#'
#' # nucleotide sequences
#' seqs <- c("TACCTGGCG","TACTAAGGG", "TACGATGAC")
#' count <- c(200,50,36)
#' amino <- FALSE
#' cor.seg.sites(seqs, count, amino)
#'
#' @seealso See \code{\link{polymorphic.sites}} for the function that calculates the number of polymorphic sites.
#' See \code{\link[ape]{seg.sites}} in \strong{ape} for another function that calculates the indices of segregating
#' (polymorphic) sites in a sample of DNA sequences.
#'
#' @import Biostrings
#'
#' @export
cor.seg.sites <- function(seqs, count, amino){
  if(length(count)<2) return(0)
  S <- polymorphic.sites(seqs, amino)$psites
  h <- length(count) #number of haplotypes
  a1 <- numeric(h-1)
  for(i in 1:(h-1)){
    a1[i] <- 1/i
  }
  a1 <- sum(a1)
  cor.seg.sites <- S/a1
  cor.seg.sites
}
