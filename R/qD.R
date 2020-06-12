#' Hill Numbers
#'
#' This function calculates the Hill number (Hill 1973) of order q, given the exponent q and
#' an individual's vector of observed read counts for each haplotype. Haplotypes may be nucleotide or
#' amino acid sequences.
#'
#' @param count A numeric vector of observed read counts for each haplotype.
#' @param q A number denoting the order of the Hill number to be returned. Any value from -\code{Inf}
#' to \code{Inf} is allowed.
#'
#' @details Hill numbers are diversity indices that give the 'effective number' of haplotypes for a
#' certain order, q, which specifies the contribution of different haplotypes to the diversity calculation.
#' As q increases, more weight is given to abundant haplotypes and the diversity calculated becomes
#' less sensitive to rare haplotypes.
#' Hill numbers are given by the equation \deqn{qD=(\sum_{i=1}^H p_i^q)^{1/(1-q)},}
#' where \eqn{H} is the number of haplotypes and \eqn{p_i} is the proportional abundance of the
#' ith haplotype.
#'
#' When q=0, the Hill number gives the number of haplotypes.
#' When q=1, each haplotype is weighted by its proportional abundance, and the Hill number
#' approaches the exponential of the Shannon entropy. This can be interpreted as the number of
#' common or typical haplotypes.
#' When q=2, the Hill number is the inverse of the Simpson index, which is the probability that
#' two sequences drawn at random from a given sample belong to the same haplotype. The Hill
#' number with q=2 can be interpreted as the number of very abundant haplotypes.
#' Hill numbers are preferred over the traditional Shannon entropy and Simpson index diversity indices
#' because Hill numbers obey the replication principle (Gregori et al. 2016).
#'
#' @return Returns the Hill number of order q as a numeric value.
#'
#' @references
#' Gregori, J., Perales, C., Rodriguez-Frias, F., Esteban, J. I., Quer, J., & Domingo, E. (2016).
#' Viral quasispecies complexity measures. \emph{Virology}, 493, 227-237.
#'
#' Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences.
#' \emph{Ecology}, 54(2), 427-432.
#'
#' @examples
#' count <- c(10,200,366,14)
#' # Hill number with q=0
#' qD(count, 0)
#' # Hill number with q=1
#' qD(count, 1)
#'
#' @seealso See \code{\link[vegan]{renyi}} in \strong{vegan} for another function that calculates Hill numbers.
#' See \code{\link{Shannon.entropy}} to calculate Shannon entropy.
#'
#' @export
qD <- function(count,q){
  if(q==0){return(length(count))}
  if(q==1){return(exp(Shannon.entropy(count)))}
  p <- count/sum(count)
  if(q==Inf){return(1/max(p))}
  if(q==-Inf){return(1/min(p))}
  res <- sum(p^q)^(1/(1-q))
  if (res==Inf & q>0){return(1/max(p))}
  if(res==Inf & q<0){return(1/min(p))}
  res
}

