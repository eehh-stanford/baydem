rm(list = ls())
library(dplyr)
library(RColorBrewer)
library(baydem)
library(evd)
library(doParallel)
library(foreach)
registerDoParallel(detectCores())

source("analysis/bayesian_radiocarbon_functions.R")

set.seed(75372) # from random.org between 1 and 1,000,000

# The simulation distribution: a Gaussian mixture with ordering pi1, pi2, mu1,
# mu2, sig1, sig2
th_sim <-
  c(
    pi1 = 0.2,
    pi2 = 0.8,
    mu1 = 775,
    mu2 = 1000,
    sig1 = 35,
    sig2 = 45
  )

hp <-
  list(
    # Class of fit (Gaussian mixture)
    fitType = "gaussmix",
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 10,
    # The gamma distribution rate parameter for sigma, yielding a mode of 100
    alpha_r = (10 - 1) / 50,
    # Minimum calendar date (years BC/AD)
    taumin = 600,
    # Maximum calendar date (years BC/AD)
    taumax = 1300,
    # Spacing for the measurement matrix (years)
    dtau = 1,
    # Number of mixtures
    K = 2
  )

# Locations for calendar date grid (spacing of 1 year)
tauVect <- seq(hp$taumin, hp$taumax, by = hp$dtau)
G <- length(tauVect)

S <- 8
TH <- matrix(NA, S, 6)

for (s in 1:S) {
  TH[s, ] <- sample_gm(hp)
}
Fmat <- bd_calc_gauss_mix_pdf_mat(TH, tauVect)

taumin <- hp$taumin
taumax <- hp$taumax

calibDf <- bd_load_calib_curve("intcal13")
tau_curve <- 1950 - calibDf$yearBP
phi_curve <- exp(-calibDf$uncalYearBP / 8033)

# Calculate the calibration curve fraction modern at the locations of the calendar grid
phiInterp <- approx(tau_curve, phi_curve, tauVect)
phiInterp <- phiInterp$y

# Use a measurement error for the fraction modern of 1e-3
measError <- 1e-3

# Determine the range of values for the fraction modern, then increase the
# range by four times the measurement error on both edges of the range
phiMin <- min(phiInterp)
phiMax <- max(phiInterp)
phiMin <- phiMin - measError * 4
phiMax <- phiMax + measError * 4
phiVect <- seq(phiMin, phiMax, len = G * 4)

# Calculate the fraction modern without adding calibration uncertainty
M <- bd_calc_meas_matrix(tauVect, phiVect, rep(measError, length(phiVect)), calibDf, addCalibUnc = F)

# Plot some sample curves (including their fraction modern probablity density)
# along with a visualization of the calibration curve
equifData <- bd_assess_calib_curve_equif(calibDf)
canInvert <- equifData$canInvert
invSpanList <- equifData$invSpanList

phiMinPlot <- 0
phiMaxPlot <- 0.08
col_vector <- mapply(brewer.pal, S, "Set1")
pdf("figures/FigS2_gm_example.pdf", width = 20, height = 12)
par(mfrow = c(3, 1))

bd_vis_calib_curve(taumin, taumax, calibDf, xlab = "Calendar Date [AD]", ylab = "Fraction Modern", invertCol = "gray80")


plot(NULL, type = "n", xlim = c(taumin, taumax), ylim = c(0, max(Fmat)), xlab = "Calendar Date [AD]", ylab = "Probability Density")

for (ii in 1:length(invSpanList)) {
  invReg <- invSpanList[[ii]]
  if (between(invReg$tau_left, taumin, taumax) || between(invReg$tau_right, taumin, taumax)) {
    rect(invReg$tau_left, 0, invReg$tau_right, max(Fmat), border = NA, col = "gray80")
  }
}

# Highlight regions 116 and 118 in all the graphs
rect(invSpanList[[116]]$tau_right, max(Fmat) - .001, invSpanList[[118]]$tau_left, max(Fmat), border = NA, col = "indianred1")

for (s in 1:S) {
  ths <- TH[s, ]
  fs <- bd_calc_gauss_mix_pdf(ths, tauVect)
  lines(tauVect, fs, col = col_vector[s], lwd = 4)
  P <- calcPerturbMatGaussMix(tauVect, ths, taumin, taumax)
  N <- MASS::Null(t(M %*% P))
  if (ncol(N) != 0) {
    stop(paste("Identifiability problem with sample", s))
  }
}

# Create matrix of fraction modern data for plotting
phiPdfMat <- matrix(NA, S, nrow(M))
for (s in 1:S) {
  phiPdfMat[s, ] <- M %*% as.matrix(Fmat[s, ])
}

pdfMin <- 0
pdfMax <- max(phiPdfMat)
plot(NULL, type = "n", xlim = c(phiMin, phiMax), ylim = c(pdfMin, pdfMax), xlab = "Fraction Modern", ylab = "Probability Density")


for (ii in 1:length(invSpanList)) {
  invReg <- invSpanList[[ii]]
  if (between(invReg$phi_left, phiMin, phiMax) || between(invReg$phi_right, phiMin, phiMax)) {
    rect(invReg$phi_left, pdfMin, invReg$phi_right, pdfMax, border = NA, col = "gray80")
  }
}


rect(invSpanList[[116]]$phi_right, max(phiPdfMat) - 10, invSpanList[[118]]$phi_left, max(phiPdfMat), border = NA, col = "indianred1")

for (s in 1:S) {
  lines(phiVect, phiPdfMat[s, ], col = col_vector[s], lwd = 4)
}

dev.off()

# Check local identifiability of simulation parameter vector
if (!is_identified(th_sim, M, tauVect, taumin, taumax)) {
  stop("Simulation parameter vector is not identified")
} else {
  print(paste("Simulation parameter vector is identified, with a measurement error of", measError))
}

# Check local identifiability for a large number of random samples
S <- 100000

identified <- rep(F, S)
print(paste("Checking local identifiability for", S, "samples"))
TH_local <- matrix(NA, S, 6)
for (s in 1:S) {
  TH_local[s, ] <- sample_gm(hp)
}

identified <- foreach(s = 1:S, .combine = cbind) %dopar% {
  output <- is_identified(TH_local[s, ], M, tauVect, taumin, taumax)
}

numBad <- sum(!identified)
print(paste0(numBad, " (out of ", S, ") non-identifiable parameter vectors found"))

# If non-identifiable parameters are found, determine whether it is caused by P
# being non-identified. If not, throw an error.
if (sum(!identified) > 0) {
  indBad <- which(!identified)
  for (ii in 1:length(indBad)) {
    s <- indBad[ii]
    print(paste("Sample", s, "is not identified"))
    ths <- TH_local[s, ]
    print(ths)
    P <- calcPerturbMatGaussMix(tauVect, ths, taumin, taumax)
    N_P <- MASS::Null(t(P))
    print(paste("The null size of P is", ncol(N_P)))
    if (ncol(N_P) == 0) {
      stop("Sample is not identified even though the null size of P zero")
    }
  }
}

numChecks <- 100000
badLocList <- list()
print(paste("Checking for non-identifiable pairs with 2 mixtures"))
start_time <- Sys.time()
print("----")
numLoc_phi <- 0
numLoc_f <- 0
numLoc_phi_and_f <- 0
numLoc_phi_not_f <- 0
numPair <- 0
relTol <- 1e-6 # relative tolerance for checking fraction modern equality
for (cc in 1:numChecks) {
  th_a <- sample_gm(hp)
  th_b <- sample_gm(hp)
  f_a <- bd_calc_gauss_mix_pdf(th_a, tauVect, taumin, taumax)
  f_b <- bd_calc_gauss_mix_pdf(th_b, tauVect, taumin, taumax)

  phiPdf_a <- M %*% f_a
  phiPdf_b <- M %*% f_b
  equalTol <- mean(c(phiPdf_a, phiPdf_b)) * relTol
  if (all(abs(phiPdf_b - phiPdf_a) <= equalTol)) {
    numPair <- numPair + 1
    print("----")
    print(p)
    print(th_a)
    print(th_b)
  }

  # Check both parameterizations for local identifiability
  Pa <- calcPerturbMatGaussMix(tauVect, th_a, taumin, taumax)
  Pb <- calcPerturbMatGaussMix(tauVect, th_a, taumin, taumax)
  Na <- MASS::Null(t(M %*% Pa))
  Na_P <- MASS::Null(t(Pa))
  Nb <- MASS::Null(t(M %*% Pb))
  Nb_P <- MASS::Null(t(Pb))

  phiBad_a <- ncol(Na) != 0
  fBad_a <- ncol(Na_P) != 0
  if (phiBad_a) {
    numLoc_phi <- numLoc_phi + 1
  }

  if (fBad_a) {
    numLoc_f <- numLoc_f + 1
  }

  if (phiBad_a && fBad_a) {
    numLoc_phi_and_f <- numLoc_phi_and_f + 1
  }

  if (phiBad_a && (!fBad_a)) {
    numLoc_phi_not_f <- numLoc_phi_not_f + 1
    badLocList[[length(badLocList) + 1]] <- th_a
  }

  phiBad_b <- ncol(Nb) != 0
  fBad_b <- ncol(Nb_P) != 0
  if (phiBad_b) {
    numLoc_phi <- numLoc_phi + 1
  }

  if (fBad_b) {
    numLoc_f <- numLoc_f + 1
  }

  if (phiBad_b && fBad_b) {
    numLoc_phi_and_f <- numLoc_phi_and_f + 1
  }

  if (phiBad_b && (!fBad_b)) {
    numLoc_phi_not_f <- numLoc_phi_not_f + 1
    badLocList[[length(badLocList) + 1]] <- th_b
  }
  end_time <- Sys.time()
}
print(paste("Finished checking in", as.character(end_time - start_time)))
print(paste("Number of pairs checked is", as.character(numChecks)))
print(paste(as.character(numPair), "observationally equivalent pairs found"))
print(paste(as.character(numLoc_phi_not_f), "parameterizations fail for phi but not f"))
print(paste(as.character(numLoc_f), "parameterizations fail for f"))
print(paste(as.character(numLoc_phi), "parameterizations fail for phi"))
print(paste(as.character(numLoc_phi_and_f), "parameterizations fail for phi and f"))
