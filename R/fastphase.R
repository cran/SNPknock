#' Calls fastPhase to fit an HMM to genotype data.
#' 
#' This function provides a wrapper for the fastPhase executable in order to fit an HMM to the
#' genotype data. The software fastPhase will fit the HMM  to the genotype data and write the 
#' corresponding parameter estimates in three separate files. 
#' Since fastPhase is not an R package, this executable must be downloaded separately by the 
#' user. Visit \url{http://scheet.org/software.html} for more information on how to obtain fastPhase.
#' 
#' @param fp_path a string with the path to the directory with the fastPhase executable.
#' @param X_file a string with the path of the genotype input file containing X in fastPhase  
#'               format (as created by \link{SNPknock.fp.writeX}).
#' @param out_path a string with the path of the directory in which the parameter estimates 
#'                 will be saved (default: NULL). If this is equal to NULL, a temporary file 
#'                 in the R temporary directory will be used.
#' @param K the number of hidden states for each haplotype sequence (default: 12).
#' @param numit the number of EM iterations (default: 25).
#' @param seed the random seed for the EM algorithm (default: 1).
#' @return A string containing the path of the directory in which the parameter estimates 
#'         were saved. This is useful to find the data when the default option for `out_path` 
#'         is used and the output is written in an R temporary directory.
#' 
#' @family fastphase
#'
#' @details
#' The software fastPhase saves the parameter estimates in three separate files whose names
#' begin with the string contained in 'out_path' and end with:
#' \itemize{
#'   \item{"_rhat.txt"}
#'   \item{"_alphahat.txt"}
#'   \item{"_thetahat.txt"}
#' }
#' 
#' The HMM for the genotype data can then be loaded from these files by calling
#' \link{SNPknock.fp.loadFit}.
#' 
#' @references 
#'   Scheet and Stephens,  A fast and flexible statistical model for large-scale population genotype data,
#'   Am J Hum Genet (2006).
#'   \href{http://www.sciencedirect.com/science/article/pii/S000292970763701X}{http://www.sciencedirect.com/science/article/pii/S000292970763701X}
#' 
#' @examples
#' fp_path  = "~/.local/bin/fastPHASE" # Path to the fastPHASE executable
#' # Specify the path to the genotype input file in ".inp" format.
#' # An example file can be found in the package installation folder.
#' X_file = system.file("extdata", "X.inp", package = "SNPknock")
#' fp_outPath = SNPknock.fp.runFastPhase(fp_path, X_file)
#' 
#' @export
SNPknock.fp.runFastPhase <- function(fp_path, X_file, out_path=NULL, K=12, numit=25, seed=1) {
  # Verify that the fastPhase executable can be found
  if(!file.exists(fp_path)) {
    message(paste("SNPknock could find the fastPhase executable: '",fp_path,"' does not exist.
If you have not downloaded it yet, you can obtain fastPhase from: http://scheet.org/software.html", sep=""))
    return(NULL)  
  }
  
  # Write to temporary directory unless specified otherwise
  if(is.null(out_path)) {
    out_path = tempfile(pattern="file", tmpdir=tempdir(), fileext = "")
  }
  
  # Make out_path absolute
  out_path_dirname = tools::file_path_as_absolute(dirname(out_path))
  out_path_basename = basename(out_path)
  out_path_abs = paste(out_path_dirname, out_path_basename, sep="/")
  
  # Prepare arguments for fastPhase
  command = fp_path
  command = paste(command, " -Pp -T1 -K", K, sep="")
  command = paste(command, " -g -H-4 -C", numit, sep="")
  command = paste(command, " -S", seed, sep="")
  command = paste(command, " -o'", out_path_abs, "' ", X_file, sep="")
  
  # Run the fastPhase executable
  tryCatch(system(command), error=function(e) 1)
  
  return(out_path)
}

#' Convert a genetic matrix X into the fastPhase input format.
#' 
#' This function convert a genetic matrix X into the fastPhase input format and saves
#' it to a user-specified file. Then, an HMM can be fitted by calling fastPhase with
#' \link{SNPknock.fp.runFastPhase}.
#' 
#' @param X a matrix of size n-by-p containing the original variables.
#' @param out_file a string containing the path of the output file onto which X will be written (default: NULL). 
#' If this is equal to NULL, a temporary file in the R temporary directory will be used.
#' @return A string containing the path of the output file onto which X was written. This is useful to find the data 
#' when the default option for `out_file` is used and X is written onto a temporary file in the R temporary directory.
#' 
#' @family fastphase
#' 
#' @references 
#'   Scheet and Stephens,  A fast and flexible statistical model for large-scale population genotype data,
#'   Am J Hum Genet (2006).
#'   \href{http://www.sciencedirect.com/science/article/pii/S000292970763701X}{http://www.sciencedirect.com/science/article/pii/S000292970763701X}
#' 
#' @examples
#' # Load an example data matrix X from the package installation directory.
#' X_file = system.file("extdata", "X.RData", package = "SNPknock")
#' load(X_file)
#' Xinp_file = SNPknock.fp.writeX(X) # Write X in a temporary file
#' 
#' @export
SNPknock.fp.writeX <- function(X, out_file=NULL) {
  # Write to temporary file unless specified otherwise
  if(is.null(out_file)) {
    out_file = tempfile(pattern="file", tmpdir=tempdir(), fileext = ".inp")
  }
  
  n = dim(X)[1]
  p = dim(X)[2]
  X = t(X)
  # Phase X (randomly)
  v1 = array(c(0,1,1))
  v2 = array(c(0,0,1))
  Xp1 = array(v1[X+1], dim(X))
  Xp2 = array(v2[X+1], dim(X))
  
  # Write to file
  con = file(out_file, "w")
  writeLines(toString(n), con = con, sep = "\n", useBytes = FALSE)
  writeLines(toString(p), con = con, sep = "\n", useBytes = FALSE)
  for(m in 1:n) {
    text = paste(c("#id",toString(m-1)), collapse = '')
    writeLines(text, con = con, sep = "\n", useBytes = FALSE)
    text = paste(Xp1[,m], collapse = '')
    writeLines(text, con = con, sep = "\n", useBytes = FALSE)
    text = paste(Xp2[,m], collapse = '')
    writeLines(text, con = con, sep = "\n", useBytes = FALSE)
  }
  close(con)
  
  return(out_file)
}

#' Load the parameter estimates obtained by fastPhase and assembles the HMM model for the genotype data.
#' 
#' This function loads the parameter estimates obtained by fastPhase (see \link{SNPknock.fp.runFastPhase})
#' and assembles the HMM model for the genotype data, in the format required by the knockoff generation function
#' \link{SNPknock.knockoffHMM}.
#' 
#' @param r_file a string with the path of the "_rhat.txt" file produced by fastPhase.
#' @param theta_file a string with the path of the "_thetahat.txt" file produced by fastPhase.
#' @param alpha_file a string with the path of the "_alphahat.txt" file produced by fastPhase.
#' @param x a numpy array of length p, where p is the number of SNPs, containing the genotype sequence of the first individual in the dataset.
#' 
#' @return A structure describing the HMM fitted by fastPhase.
#' 
#' @family fastphase
#' 
#' @details
#' In the description of the parameter "x", the first individual is intended 
#' in the same order as provided to fastPhase. This is needed in order to correctly 
#' intepret the emission parameters estimated by fastPhase.
#'
#' This function returns a structure with three fields: 
#' \itemize{
#'   \item{"pInit": an array of length K, containing the marginal distribution of the hidden states for the first SNP.}
#'   \item{"Q": an array of size (p-1,K,K), containing a list of p-1 transition matrices between the K latent states of the HMM.}
#'   \item{"pEmit": an array of size (p,K,3), containing the emission probabilities of the hidden states for each of the p SNPs.}
#'  }
#'   
#' @references 
#'   Scheet and Stephens,  A fast and flexible statistical model for large-scale population genotype data,
#'   Am J Hum Genet (2006).
#'   \href{http://www.sciencedirect.com/science/article/pii/S000292970763701X}{http://www.sciencedirect.com/science/article/pii/S000292970763701X}
#' 
#' @examples
#' # Load an example data matrix X from the package installation directory.
#' X_file = system.file("extdata", "X.RData", package = "SNPknock")
#' load(X_file)
#' 
#' # Specify the location of the fastPhase output files containing the parameter estimates.
#' # Example files can be found in the package installation directory.
#' r_file = system.file("extdata", "X_rhat.txt", package = "SNPknock")
#' theta_file = system.file("extdata", "X_thetahat.txt", package = "SNPknock")
#' alpha_file = system.file("extdata", "X_alphahat.txt", package = "SNPknock")
#' 
#' # Read the parameter files and build the HMM
#' hmm = SNPknock.fp.loadFit(r_file, theta_file, alpha_file, X[1,])
#' 
#' @export
SNPknock.fp.loadFit <- function(r_file, theta_file, alpha_file, x) {
  # Load (r,theta,alpha) paramters from fastPhase fit
  r = loadEMParameters(r_file)
  alpha = loadEMParameters(alpha_file)
  theta = loadEMParameters(theta_file)

  # Flip theta
  X_chr_flip = x>0
  theta[X_chr_flip,] = 1-theta[X_chr_flip,]
  
  # Assemble transition matrices
  Q1 = compute_Q1(r, alpha)
  Q = assemble_Q(Q1)
  pEmit = assemble_pEmit(theta)
  pInit = assemble_pInit(alpha)
  
  hmm = NULL
  hmm$pInit = pInit
  hmm$Q = Q
  hmm$pEmit = pEmit
  return(hmm) 
}

#' Read the files produced by fastPhase
#'  
#' @rdname loadEMParameters
#' @keywords internal
loadEMParameters <- function(data_path) {
  # Load the entire file
  lines = readLines(data_path)
  # Remove first comment row
  lines = lines[2:length(lines)]
  # Keep only first EM start
  last_row = suppressWarnings(min(grep("^\\s*>", lines)))
  if (is.finite(last_row)) {
    lines = lines[1:(last_row-1)]
  }
  # Read the relevant rows as a csv dile
  con = textConnection(lines)
  return(as.matrix(utils::read.table(con)))
}

#' Compute the haplotype transition matrices based on the fastPhase HMM
#'  
#' @rdname compute_Q1
#' @keywords internal
compute_Q1 <- function(r, alpha) {
  p = dim(alpha)[1]
  K = dim(alpha)[2]
  Q = array(rep(0, (p-1)*K*K), c(p-1,K,K))
  rExp = exp(-r)
  for(j in 2:p) {
    v = (1-rExp[j])*alpha[j,]
    Q[j-1,,] = t(matrix(rep(v,K), ncol = K)) + diag(rep(rExp[j],K))
  }
  return(Q)
}

#' Compute the genotype transition matrices based on the fastPhase HMM
#'  
#' @rdname assemble_Q
#' @keywords internal
assemble_Q <- function(Q1) {
  p = dim(Q1)[1]+1
  K = dim(Q1)[2]
  Keff = K*(K+1)/2
  Q = array(rep(0, (p-1)*Keff*Keff), c(p-1,Keff,Keff))
  for (k1 in 1:K) {
    for (k2 in 1:k1) {
      i = ((k1-1)*k1)/2+k2
      for (k1p in 1:K) {
        for (k2p in 1:k1p) {
          j = ((k1p-1)*k1p)/2+k2p
          Q[,i,j] = Q1[,k1,k1p] * Q1[,k2,k2p]
          if (k1p != k2p) {# Normalization seems correct, but different from the paper
            Q[,i,j] = Q[,i,j] + Q1[,k1,k2p] * Q1[,k2,k1p]
          }
        }
      }
    }
  }
  for(m in 1:(p-1)) {
    # Enforce normalization
    QSums = rowSums(Q[m,,])
    #stopifnot(abs(QSums-1)<1e-6)
    Q[m,,] = sweep(Q[m,,],1,QSums,`/`)
  }
  return(Q)
}

#' Compute the genotype emission distributions based on the fastPhase HMM
#'  
#' @rdname assemble_pEmit
#' @keywords internal
assemble_pEmit<- function(theta) {
  p = dim(theta)[1]
  K = dim(theta)[2]
  Keff = K*(K+1)/2
  pEmit = array(rep(0,p*3*Keff), c(p,3,Keff))
  for(m in 1:p) {
    pEmit1 = array(rep(0,3*Keff), c(3,Keff))
    for (k1 in 1:K) {
      for (k2 in 1:k1) {
        i = ((k1-1)*k1)/2+k2
        pEmit1[1,i] = (1-theta[m,k1])*(1-theta[m,k2])
        pEmit1[2,i] = theta[m,k1]*(1-theta[m,k2]) + theta[m,k2]*(1-theta[m,k1])
        pEmit1[3,i] = theta[m,k1]*theta[m,k2] 
      }
    }
    # Enforce normalization
    pEmit1Sums = colSums(pEmit1)
    # stopifnot(abs(pEmit1Sums-1)<1e-6)
    pEmit[m,,] = sweep(pEmit1,2,pEmit1Sums,`/`)
  }
  return(pEmit)
}

#' Compute the genotype initial distributions based on the fastPhase HMM
#'  
#' @rdname assemble_pInit
#' @keywords internal
assemble_pInit <- function(alpha) {
  p = dim(alpha)[1]
  K = dim(alpha)[2]
  Keff = K*(K+1)/2
  pInit = rep(0,Keff)
  for (k1 in 1:K) {
    for (k2 in 1:k1) {
      i = ((k1-1)*k1)/2+k2
      if (k1==k2)
        pInit[i] = alpha[1,k1]*alpha[1,k1]
      else
          pInit[i] = 2*alpha[1,k1]*alpha[1,k2]
    }
  }
  # Enforce normalization
  pInitSum = sum(pInit)
  # stopifnot(abs(pInitSum-1)<1e-4)
  pInit = pInit/pInitSum
  
  return(pInit)
}