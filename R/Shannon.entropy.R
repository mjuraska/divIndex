#' Shannon Entropy
#'
#' This function calculates the Shannon entropy diversity index, given a vector of the observed
#' read counts for each haplotype. Haplotypes may be nucleotide or amino acid sequences.
#'
#' @param count A numeric vector of read counts for each haplotype.
#'
#' @details The shannon entropy is a diversity index that measures the uncertainty in predicting
#' the specific haplotype of a randomly sampled sequence in an individual. Maximum uncertainty occurs
#' when all haplotypes are equally abundant. Minimum uncertainty occurs when there is only one haplotype;
#' in this case, the Shannon entropy has a value of 0. Each haplotype is weighed exactly by its frequency.
#' Shannon entropy is given by the equation \deqn{H' = -\sum p_i ln(p_i),}
#' where \eqn{p_(i)} is the proportional abundance of the ith haplotype.
#'
#' When the number of sequence counts varies drastically from individual to individual, sample size
#' bias may occur, where individuals with larger sequence counts appear to have higher diversity.
#' The Shannon entropy is sensitive to this sample size bias, although the bias can be partially
#' corrected by fringe-trimming (Gregori et al., 2014) and by a Taylor series expansion:
#' \eqn{H'_{MLE} = H' - \frac{H-1}{2n} + \frac{1-\sum p_i^{-1}}{12n^2} + \dots}.
#'
#' @return Returns the Shannon entropy as a numeric value.
#'
#' @examples
#' count <- c(10,200,366,14)
#' Shannon.entropy(count)
#'
#' @seealso See \code{\link[vegan]{diversity}} in \strong{vegan} for another function that calculates Shannon entropy.
#'
#' @export
Shannon.entropy <- function(count){
  h <- length(count)
  if(h<2){
    return(0)
  }
  p <- count/sum(count)
  lnp <- ifelse(count==0,0,log(p))
  S <- -sum(p*lnp)
  ###  Correct bias only if vector of counts given
  if(all(count >= 1)){
    S <- S-(h-1)/(2*sum(count))
  }
  if(S>log(h)){
    S <- log(h)
  }
  return(S)
}
