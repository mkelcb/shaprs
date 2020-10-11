##############################################################
# Summary stats pre-processor script:
# it performs:
# (1) modified Cochran's Q-test which optionally adjusts for overlapping controls
# (2) lFDR estimation on the p-values from the above Cochran's Q-test
# (3) lists SNPs that are fail the heterogeneity test at specified thresholds (optional)
# INPUT:
# (1) inputDataLoc: summary statistics file with the following structure (requires header):
#       SNP	CHR	BP	Beta_A	SE_A	Beta_B	SE_B
# rs4040617   1  779322 -0.0017630 0.008608 -0.010990 0.008592
# (2) rho: overlap between studies. 0 for no overlap and 1 for complete overlap
# (3) outputLoc: location of the output files
# (4) thresholds: space separated list of thresholds to be used to create list of SNPs (optional)
# OUTPUT:
# (1) lFDR table, a 2 column file with the following signature
#     SNPID lFDR
# (2) a file to hold the list of SNPids that failed the heterogeneity test at each threshold (optional)


##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 3) {stop("not enough Arguments received")} 

inputDataLoc = args[1]
rho = as.numeric(args[2])
outputLoc = args[3]
thresholds = vector(length = (length(args) -3) )
counter = 1
if (length(args) > 3) {
  for (i in 4:length(args)) { # loop through the rest of the arguments, where each is supposed to be a threshold
    thresholds[counter] = as.numeric(args[i])
    counter = counter +1
    print(thresholds[counter])
  }
}


##############################################################
#                       lFDR FUNCTION                        #
##############################################################
# stand-alone version of the one implemented by https://github.com/StoreyLab/qvalue
lfdr = function (p, pi0 = NULL, trunc = TRUE, monotone = TRUE, transf = c("probit", 
                                                                   "logit"), adj = 1.5, eps = 10^-8, ...) 
{
  lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("P-values not in valid range [0,1].")
  }
  else if (is.null(pi0)) {
    pi0 <- pi0est(p, ...)$pi0
  }
  n <- length(p)
  transf <- match.arg(transf)
  if (transf == "probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1 - eps)
    x <- qnorm(p)
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    lfdr <- pi0 * dnorm(x)/y
  }
  else {
    x <- log((p + eps)/(1 - p + eps))
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x)/(1 + exp(x))^2
    lfdr <- (pi0 * dx)/y
  }
  if (trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if (monotone) {
    o <- order(p, decreasing = FALSE)
    ro <- order(o)
    lfdr <- cummax(lfdr[o])[ro]
  }
  lfdr_out[rm_na] <- lfdr
  return(lfdr_out)
}

pi0est = function (p, lambda = seq(0.05, 0.95, 0.05), pi0.method = c("smoother", 
                                                            "bootstrap"), smooth.df = 3, smooth.log.pi0 = FALSE, 
          ...) 
{
  rm_na <- !is.na(p)
  p <- p[rm_na]
  pi0.method = match.arg(pi0.method)
  m <- length(p)
  lambda <- sort(lambda)
  ll <- length(lambda)
  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  }
  else if (ll > 1 && ll < 4) {
    stop(sprintf(paste("ERROR:", paste("length(lambda)=", 
                                       ll, ".", sep = ""), "If length of lambda greater than 1,", 
                       "you need at least 4 values.")))
  }
  else if (min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }
  if (max(p) < max(lambda)) {
    stop("ERROR: maximum p-value is smaller than lambda range. Change the range of lambda or use qvalue_truncp() for truncated p-values.")
  }
  if (ll == 1) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0.lambda <- pi0
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL
  }
  else {
    ind <- length(lambda):1
    pi0 <- cumsum(tabulate(findInterval(p, vec = lambda))[ind])/(length(p) * 
                                                                   (1 - lambda[ind]))
    pi0 <- pi0[ind]
    pi0.lambda <- pi0
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      }
      else {
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    }
    else if (pi0.method == "bootstrap") {
      minpi0 <- quantile(pi0, prob = 0.1)
      W <- sapply(lambda, function(l) sum(p >= l))
      mse <- (W/(m^2 * (1 - lambda)^2)) * (1 - W/m) + (pi0 - 
                                                         minpi0)^2
      pi0 <- min(pi0[mse == min(mse)], 1)
      pi0Smooth <- NULL
    }
    else {
      stop("ERROR: pi0.method must be one of \"smoother\" or \"bootstrap\".")
    }
  }
  if (pi0 <= 0) {
    warning("The estimated pi0 <= 0. Setting the pi0 estimate to be 1. Check that you have valid p-values or use a different range of lambda.")
    pi0 <- pi0.lambda <- 1
    pi0Smooth <- lambda <- 0
  }
  return(list(pi0 = pi0, pi0.lambda = pi0.lambda, lambda = lambda, 
              pi0.smooth = pi0Smooth))
}


estimate_lFDR = function (p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, 
                     pi0 = NULL, ...) 
{
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  }
  else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 
                                   1)) {
    stop("'fdr.level' must be in (0, 1].")
  }
  if (is.null(pi0)) {
    pi0s <- pi0est(p, ...)
  }
  else {
    if (pi0 > 0 && pi0 <= 1) {
      pi0s = list()
      pi0s$pi0 = pi0
    }
    else {
      stop("pi0 is not (0,1]")
    }
  }
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m/(i * (1 - 
                                                        (1 - p[o])^m))))[ro]
  }
  else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m/i))[ro]
  }
  qvals_out[rm_na] <- qvals
  if (lfdr.out) {
    lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  }
  else {
    lfdr_out <- NULL
  }
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out, 
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level, 
                   significant = (qvals <= fdr.level), pi0.lambda = pi0s$pi0.lambda, 
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out, 
                   pvalues = p_in, lfdr = lfdr_out, pi0.lambda = pi0s$pi0.lambda, 
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}

##############################################################
#                    lFDR FUNCTION  END                      #
##############################################################


##############################################################
#                         MAIN SCRIPT                        #
##############################################################

# 1. load data
inputData= read.table(inputDataLoc, header = T) 

# 2. Cochran's Q-test formula: from 'Meta-Analysis of Genome-wide Association Studies with Overlapping Subjects' ncbi.nlm.nih.gov/pmc/articles/PMC2790578/
Vhat = (inputData$SE_A^2 + inputData$SE_B^2 - 2 * rho * inputData$SE_A * inputData$SE_B)
Q_adjusted = (inputData$Beta_A - inputData$Beta_B)^2/ Vhat
df=2-1 # degrees of freedom, 2 studies -1
Q_adjusted_pvals = pchisq(Q_adjusted, df = df, lower.tail = F)

# 3. lFDR estimation
lfdr_obj = estimate_lFDR(p = Q_adjusted_pvals)
lfdr_qvals <- lfdr_obj$lfdr

# 4. write out a table of lFDR values for each SNP
lFDRTable <- data.frame(inputData$SNP, lfdr_qvals)
filename = paste(outputLoc, "_SNP_lFDR" , sep="")
write.table(lFDRTable, filename, row.names = F, col.names = F, quote = FALSE) 
print(paste("written lFDRs for SNPs to",filename))

# 5. Write out SNPs that fail the heterogeneity test at the specified thresholds
for (i in 1:length(thresholds)) {
  currentThreshold=thresholds[i]
  sigSNPs = lFDRTable[which(lFDRTable$lfdr_qvals <currentThreshold),]
  filename = paste(outputLoc, "_SNPs_",currentThreshold , sep="")
  write.table(sigSNPs$inputData.SNP, filename, row.names = F, col.names = F, quote = FALSE) 
  print(paste("written ", nrow(sigSNPs), " SNPs to:", filename, sep=""))
}

