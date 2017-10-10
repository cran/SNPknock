#' Sample knockoff copies of a discrete Markov chain.
#' 
#' This function constructs knockoff copies of variables distributed as a discrete Markov chain.
#' 
#' @param X a matrix of size n-by-p containing the original variables.
#' @param pInit an array of length K, containing the marginal distribution of the states for the first variable.
#' @param Q an array of size (p-1,K,K), containing a list of p-1 transition matrices between the K states of the Markov chain.
#' @param seed an integer random seed (default: 123).
#' @return A matrix of size n-by-p containing the knockoff variables.
#' 
#' @family knockoffs
#' 
#' @details
#' Each element of the matrix X should be an integer value between 0 and K-1.
#' The transition matrices contained in Q are defined such that \eqn{P[X_{j+1}=k|X_{j}=l]=Q[j,l,k]}.
#' 
#' @references 
#'   Sesia et al., Gene Hunting with Knockoffs for Hidden Markov Models,
#'   arXiv:1706.04677 (2017).
#'   \href{https://statweb.stanford.edu/~candes/papers/HMM_Knockoffs.pdf}{https://statweb.stanford.edu/~candes/papers/HMM_Knockoffs.pdf}
#' 
#' @examples
#' p=10; K=5;
#' pInit = rep(1/K,K)
#' Q = array(stats::runif((p-1)*K*K),c(p-1,K,K))
#' for(j in 1:(p-1)) { Q[j,,] = Q[j,,] / rowSums(Q[j,,]) }
#' X = SNPknock.models.sampleDMC(pInit, Q, n=20)
#' Xk = SNPknock.knockoffDMC(X, pInit, Q)
#' 
#' @export
SNPknock.knockoffDMC <- function(X, pInit, Q, seed=123) {
  # Verify dimensions are compatible
  stopifnot(dim(X)[2]==dim(Q)[1]+1)
  stopifnot(length(pInit)==dim(Q)[2])
  stopifnot(dim(Q)[2]==dim(Q)[3])
  
  # Verify contents are compatible
  stopifnot(is.numeric(X))
  stopifnot(is.numeric(pInit))
  stopifnot(is.numeric(Q))
  stopifnot(is.numeric(seed))
  stopifnot(seed==floor(seed))
  stopifnot( all(X == floor(X)) )
  stopifnot( min(X)>=0 )
  stopifnot( max(X)<length(pInit) )
  
  # Extract dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  K = length(pInit)
  
  # Call the CPP wrapper
  Xk = knockoffDMC_wrapper(X, pInit, Q, n, p, K, seed)
}

#' Sample knockoff copies of a hidden Markov model.
#' 
#' This function constructs knockoff copies of variables distributed as a hidden Markov model.
#' 
#' @param X a matrix of size n-by-p containing the original variables.
#' @param pInit an array of length K, containing the marginal distribution of the states for the first variable.
#' @param Q an array of size (p-1,K,K), containing a list of p-1 transition matrices between the K states of the Markov chain.
#' @param pEmit an array of size (p,M,K), containing the emission probabilities for each of the M possible emission states, 
#' from each of the K hidden states and the p variables.
#' @param seed an integer random seed (default: 123).
#' @return A matrix of size n-by-p containing the knockoff variables.
#' 
#' @family knockoffs
#' 
#' @details
#' Each element of the matrix X should be an integer value between 0 and M-1.
#' The transition matrices contained in Q are defined with the same convention as in \link{SNPknock.knockoffDMC}.
#' The emission propability matrices contained in pEmit are defined such that \eqn{P[X_{j}=k|H_{j}=l]=\mathrm{pEmit}[j,k,l]},
#' where \eqn{H_j} is the latent variable associated to \eqn{X_j}.
#' 
#' @references 
#'   Sesia et al., Gene Hunting with Knockoffs for Hidden Markov Models,
#'   arXiv:1706.04677 (2017).
#'   \href{https://statweb.stanford.edu/~candes/papers/HMM_Knockoffs.pdf}{https://statweb.stanford.edu/~candes/papers/HMM_Knockoffs.pdf}
#'   
#' @examples
#' p=10; K=5; M=3;
#' pInit = rep(1/K,K)
#' Q = array(stats::runif((p-1)*K*K),c(p-1,K,K))
#' for(j in 1:(p-1)) { Q[j,,] = Q[j,,] / rowSums(Q[j,,]) }
#' pEmit = array(stats::runif(p*M*K),c(p,M,K))
#' for(j in 1:p) { pEmit[j,,] = pEmit[j,,] / rowSums(pEmit[j,,]) }
#' X = SNPknock.models.sampleHMM(pInit, Q, pEmit, n=20)
#' Xk = SNPknock.knockoffHMM(X, pInit, Q, pEmit)
#'    
#' @export
SNPknock.knockoffHMM <- function(X, pInit, Q, pEmit, seed=123) {
  # Verify dimensions are compatible
  stopifnot(dim(X)[2]==dim(Q)[1]+1)
  stopifnot(length(pInit)==dim(Q)[2])
  stopifnot(dim(pEmit)[3]==dim(Q)[2])
  stopifnot(dim(pEmit)[1]==dim(Q)[1]+1)
  stopifnot(dim(Q)[2]==dim(Q)[3])
  
  # Verify contents are compatible
  stopifnot(is.numeric(X))
  stopifnot(is.numeric(pInit))
  stopifnot(is.numeric(Q))
  stopifnot(is.numeric(pEmit))
  stopifnot(is.numeric(seed))
  stopifnot(seed==floor(seed))
  stopifnot( all(X == floor(X)) )
  stopifnot( min(X)>=0 )
  stopifnot( max(X)<length(pInit) )
  
  # Extract dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  K = length(pInit)
  M = dim(pEmit)[2]
  
  # Call the CPP wrapper
  Xk = knockoffHMM_wrapper(X, pInit, Q, pEmit, n, p, K, M, seed)
}