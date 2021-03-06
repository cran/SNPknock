% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastphase.R
\name{runFastPhase}
\alias{runFastPhase}
\title{Fit an HMM to genetic data using fastPHASE}
\usage{
runFastPhase(fp_path, X_file, out_path = NULL, K = 12, numit = 25,
  phased = FALSE, seed = 1)
}
\arguments{
\item{fp_path}{a string with the path to the directory with the fastPHASE executable.}

\item{X_file}{a string with the path of the genotype input file containing X in fastPHASE
format (as created by \link{writeXtoInp}).}

\item{out_path}{a string with the path of the directory in which the parameter estimates
will be saved (default: NULL). If this is equal to NULL, a temporary file
in the R temporary directory will be used.}

\item{K}{the number of hidden states for each haplotype sequence (default: 12).}

\item{numit}{the number of EM iterations (default: 25).}

\item{phased}{whether the data are already phased (default: FALSE).}

\item{seed}{the random seed for the EM algorithm (default: 1).}
}
\value{
A string containing the path of the directory in which the parameter estimates
        were saved. This is useful to find the data when the default option for `out_path`
        is used and the output is written in an R temporary directory.
}
\description{
This function provides a wrapper for the fastPHASE executable in order to fit an HMM to either
unphased genotype data or phased haplotype data.
The software fastPHASE will fit the HMM  to the genotype data and write the
corresponding parameter estimates in four separate files.
Since fastPHASE is not an R package, this executable must be downloaded separately by the
user. Visit \url{http://scheet.org/software.html} for more information on how to obtain fastPHASE.
}
\details{
The software fastPHASE saves the parameter estimates in four separate files whose names
begin with the string contained in 'out_path' and end with:
\itemize{
  \item{"_rhat.txt"}
  \item{"_alphahat.txt"}
  \item{"_thetahat.txt"}
  \item{"_origchars"}
}

The HMM for the genotype data can then be loaded from these files by calling
\link{loadHMM}.
}
\examples{
fp_path  = "~/bin/fastPHASE" # Path to the fastPHASE executable

# Run fastPHASE on unphased genotypes
# Specify the path to the genotype input file in ".inp" format.
# An example file containing unphased genotypes can be found in the package installation folder.
X_file = system.file("extdata", "genotypes.inp", package = "SNPknock")
fp_outPath = runFastPhase(fp_path, X_file)

# Run fastPHASE on phased haplotypes
# An example file containing phased haplotypes can be found in the package installation folder.
H_file = system.file("extdata", "haplotypes.inp", package = "SNPknock")
fp_outPath = runFastPhase(fp_path, H_file, phased=TRUE)

}
\references{
\insertRef{scheet2006}{SNPknock}
}
\seealso{
Other fastPHASE: \code{\link{loadHMM}},
  \code{\link{writeXtoInp}}
}
\concept{fastPHASE}
