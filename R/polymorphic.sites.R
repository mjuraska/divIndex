#' Polymorphic Sites
#'
#' This function computes the number of polymorphic sites, the number of unique mutations, and
#' the locations of polymorphic sites, given a vector of an individual's
#' haplotype sequences. Haplotypes may be nucleotide or amino acid sequences.
#'
#' @param seqs A character vector containing the distinct haplotype sequences from an individual subject.
#' @param amino A logical value indicating whether haplotypes are amino acid sequences or nucleotide
#' sequences. If \code{FALSE}, nucleotide sequences are assumed.
#'
#' @details The number of polymorphic sites is used in the calculation of the corrected number of
#' segregating sites, a genetic diversity index that estimates the total number of polymorphic sites among
#' the haplotypes within an individual.
#'
#' @return Returns a list with the following items:
#'   \item{\code{psites}}{Number of polymorphic sites, defined as the number of positions with at least
#'   one base or amino acid.}
#'   \item{\code{nmuts}}{Number of unique mutations, assumed with reference to a consensus sequence. For
#'   each position, there is assumed to be one reference base or amino acid. All other bases or amino acids
#'   appearing in that position are counted as unique mutations.}
#'   \item{\code{poly.pos}}{Numeric vector of the numbered positions of polymorphic sites.}
#'
#' @references
#' Gregori, J., Perales, C., Rodriguez-Frias, F., Esteban, J. I., Quer, J., & Domingo, E. (2016).
#' Viral quasispecies complexity measures. \emph{Virology}, 493, 227-237.
#'
#' @examples
#' seqs <- c("ACDEFGAI","KLMNPQAS","TVWYACAE","SMWLWCAN")
#' polymorphic.sites(seqs, amino=TRUE)
#'
#' seqs <- c("TACCTGGCG","TACTAAGGG", "TACGATGAC")
#' polymorphic.sites(seqs, amino=FALSE)
#'
#' @seealso See \code{\link[ape]{seg.sites}} in \strong{ape} for another function that calculates the
#' indices of segregating (polymorphic) sites in a sample of DNA sequences. See \code{\link{cor.seg.sites}}
#' for a function that estimates the number of polymorphic sites, correcting for sample size bias.
#'
#' @import Biostrings
#'
#' @export
polymorphic.sites <- function(seqs, amino){
  nsq <- length(seqs)
  if(amino==TRUE){
    frqmat <- consensusMatrix(AAStringSet(seqs)) # frequency matrix
  } else{
    frqmat <- consensusMatrix(DNAStringSet(seqs))[c("A","C","G","T"),]
  }
  psites <- apply(frqmat,2,function(mut) sum(mut>0)) # number of bases or amino acids per site
  poly.pos <- which(psites>1) # polymorphic sites (more than one base or amino acid)
  psites <- psites[poly.pos] # number of bases or amino acids at each polymorphic site
  tbl <- table(psites)
  nmuts <- sum(psites-1) # number of unique mutations
  psites <- length(psites) # number of polymorphic sites
  list(psites=psites,nmuts=as.integer(nmuts),poly.pos=poly.pos)
}
