## ------------------------------------------------------------------------
library(SNPknock)
X_file = system.file("extdata", "X.RData", package = "SNPknock")
load(X_file)
table(X)

## ------------------------------------------------------------------------
# Convert X into the suitable fastPhase input format, write it into a temporary file
# and return the path to that file.
Xinp_file = SNPknock.fp.writeX(X)

## ------------------------------------------------------------------------
fp_path  = "~/bin/fastPHASE" # Path to the fastPHASE executable
# Call fastPhase and return the path to the parameter estimate files
fp_outPath = SNPknock.fp.runFastPhase(fp_path, Xinp_file)

## ----eval=FALSE----------------------------------------------------------
#  r_file = paste(fp_outPath, "_rhat.txt", sep="")
#  theta_file = paste(fp_outPath, "_thetahat.txt", sep="")
#  alpha_file = paste(fp_outPath, "_alphahat.txt", sep="")

## ------------------------------------------------------------------------
r_file = system.file("extdata", "X_rhat.txt", package = "SNPknock")
theta_file = system.file("extdata", "X_thetahat.txt", package = "SNPknock")
alpha_file = system.file("extdata", "X_alphahat.txt", package = "SNPknock")

## ------------------------------------------------------------------------
hmm = SNPknock.fp.loadFit(r_file, theta_file, alpha_file, X[1,])

## ------------------------------------------------------------------------
Xk = SNPknock.knockoffHMM(X, hmm$pInit, hmm$Q, hmm$pEmit)
table(Xk)

