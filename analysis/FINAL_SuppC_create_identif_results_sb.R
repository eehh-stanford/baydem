library(baydem)
library(RColorBrewer)

source(here::here("analysis/bayesian_radiocarbon_functions.R"))

# ****  null_analysis.pdf  ****
# Plot the calibration curve and null vectors on the interval AD 600 to 1300
taumin <- 600 # AD. Minimum calendar date
taumax <- 1300 # AD. Maximum calendar date
dtau <- 5 # Spacing in years of sampling grid

# To assess only the influence of the shape of the calibration curve, use a
# very small error and (below) do not add calibration curve uncertainty to the
# measurement matrix
measError <- .0001
errorSpec <- list(min = measError, max = measError)

# Load the calibration curve and call the code to assess equifinality
calibDf <- bd_load_calib_curve("intcal13")
tau_curve <- 1950 - calibDf$yearBP # AD
equiList <- bd_calc_calib_curve_equif_dates(calibDf)
phi_curve <- bd_calc_calib_curve_frac_modern(calibDf)
temp <- bd_assess_calib_curve_equif(calibDf)
canInvert <- temp$canInvert
invSpanList <- temp$invSpanList

# Determine which value of tau and phi to include in the measurement matrix. All
# the phi values on the interval taumin to taumax in the calibration curve are
# used. These values correspond to one or more calendar dates as summarized in
# the list of invertible regions, invSpanList. Use all these tau-values. Hence,
# there are more tau-values than phi-values: length(tauVect) > length(phiVect).
ind <- (tau_curve >= taumin) & (tau_curve <= taumax)
phiVect <- unique(phi_curve[ind])

tauVect <- tau_curve[ind]
for (kk in 1:length(equiList)) {
  equiEntry <- equiList[[kk]]
  if (equiEntry$indBase %in% which(ind)) {
    tauVect <- unique(c(tauVect, equiEntry$tauBase))
    tauVect <- unique(c(tauVect, equiEntry$tauEqui))
  }
}

tauVect <- sort(tauVect) # Order the tau-values

tauVect <- tauVect[tauVect >= taumin]
tauVect <- tauVect[tauVect <= taumax]

# Calculate the measurement matrix. To assess, for now, only the influence of
# the shape of the calibration curve, do not add calibration curve uncertainty
# in the measurement matrix calculation. calc_meas_matrix2 is identical to
# bd_calc_meas_matrix, except that no error checking is done to ensure that the
# grid spacing is even (this is because sometimes it is not exactly even to
# floating point accuracy).
M <- bd_calc_meas_matrix(tauVect, phiVect, rep(measError, length(phiVect)), calibDf, F, useTrapez = T)
# Calculate the null space of the measurement matrix.
N <- MASS::Null(t(M))
rankM <- Matrix::rankMatrix(M)
print("Dimension of null matrix is")
print(dim(N))
print("Rank of measurement matrix is")
print(rankM)
print("Dimension of measurement matrix is")
print(dim(M))

# Make sure the rank and null dimensions are consistent
if (rankM + ncol(N) != ncol(M)) {
  stop("rank + null dimension not equal to number of columns in M")
}

# Save named pairs of key information for reporting in the manuscript
outlog <- list()
outlog$`taumin` <- taumin
outlog$`taumax` <- taumax
outlog$`nullSize` <- ncol(N)
outlog$`numRow_M` <- nrow(M)
outlog$`numCol_M` <- ncol(M)
outlog$`rankM` <- rankM
outlog$`numPoints_tau` <- length(tauVect)
outlog$`numPoints_phi` <- length(phiVect)

# Plot the calibration curve
tauplot <- tau_curve[ind]
canInvert <- canInvert[ind]
pdf(here::here("analysis/figures/FigS2_null_analysis.pdf"), width = 20, height = 12)
par(mfrow = c(2, 1))
bd_vis_calib_curve(taumin, taumax, calibDf, invertCol = "gray80", xlab = "Calendar Date [AD]", ylab = "Fraction Modern")

# In the top plot, draw horizontal dashed lines for cluster 118 in invSpanList
invReg <- invSpanList[[118]]
lines(c(500, invReg$tau_left), c(1, 1) * invReg$phi_left, col = "black", lwd = 2)
lines(c(500, tau_curve[invReg$ii_next]), c(1, 1) * invReg$phi_right, col = "black", lwd = 2)

# Plot the null vectors
plot(NULL, type = "n", xlim = c(taumin, taumax), ylim = c(-1, 1), xlab = "Calendar Date [AD]", ylab = "Null Vector Value")

# First, add grey bands for the null vector plot to mark invertible regions of
# the calibration curve
for (ii in 1:length(invSpanList)) {
  invReg <- invSpanList[[ii]]
  if (dplyr::between(invReg$tau_left, taumin, taumax) || dplyr::between(invReg$tau_right, taumin, taumax)) {
    rect(invReg$tau_left, -1, invReg$tau_right, 1, border = NA, col = "gray80")
  }
}

# Second, plot the null vectors
pointSize <- .75
for (i in 2:ncol(N)) { # All but the first
  points(tauVect, N[, i], col = "black", cex = pointSize, pch = 19)
}
# Plot the first null vector in red
points(tauVect, N[, 1], col = "red", pch = 19, cex = pointSize)
dev.off()

# ****  standard_basis_exp.pdf  ****
# Add draws from the null space of the standard basis to an exponential growth
# curve for a part of the calibration curve with one non-invertible region.
taumin2 <- 970
taumax2 <- 1035

tau_curve2 <- 1950 - calibDf$yearBP # AD
equiList2 <- bd_calc_calib_curve_equif_dates(calibDf)
phi_curve2 <- bd_calc_calib_curve_frac_modern(calibDf)
temp2 <- bd_assess_calib_curve_equif(calibDf)
canInvert2 <- temp2$canInvert
invSpanList2 <- temp2$invSpanList

# Determine which value of tau and phi to include in the measurement matrix. All
# the phi values on the interval taumin to taumax in the calibration curve are
# used. These values correspond to one or more calendar dates as summarized in
# the list of invertible regions, invSpanList. Use all these tau-values. Hence,
# there are more tau-values than phi-values: length(tauVect) > length(phiVect).
ind2 <- (tau_curve2 >= taumin2) & (tau_curve2 <= taumax2)
phiVect2 <- unique(phi_curve2[ind2])


# A function to find invertible regions in the span taumin to taumax
findInvReg <- function(invSpanList, taumin, taumax) {
  count <- 0
  reg <- list()
  for (ii in 1:length(invSpanList)) {
    invReg <- invSpanList[[ii]]
    if (dplyr::between(invReg$tau_left, taumin, taumax) || dplyr::between(invReg$tau_right, taumin, taumax)) {
      count <- count + 1
      reg[[count]] <- invReg
    }
  }
  return(reg)
}

# Find invertible regions (used below to make plot) and make sure there are two
reg <- findInvReg(invSpanList2, taumin2, taumax2)
if (length(reg) != 2) {
  stop("There should be two invertible regions for this example")
}

tauVect2 <- tau_curve2[ind2]
counts <- 0
for (kk in 1:length(equiList2)) {
  equiEntry2 <- equiList2[[kk]]
  if (equiEntry2$indBase %in% which(ind2)) {
    tauVect2 <- unique(c(tauVect2, equiEntry2$tauBase))
    tauVect2 <- unique(c(tauVect2, equiEntry2$tauEqui))
  }
}

tauVect2 <- sort(tauVect2) # Order the tau-values

# Calculate the measurement matrix. To assess, for now, only the influence of
# the shape of the calibration curve, do not add calibration curve uncertainty
# in the measurement matrix calculation. calc_meas_matrix2 is identical to
# bd_calc_meas_matrix, except that no error checking is done to ensure that the
# grid spacing is even (this is because sometimes it is not exactly even to
# floating point accuracy).
M2 <- bd_calc_meas_matrix(tauVect2, phiVect2, rep(measError, length(phiVect2)), calibDf, F, useTrapez = T)
# Calculate the null space of the measurement matrix.
N2 <- MASS::Null(t(M2))
rankM2 <- Matrix::rankMatrix(M2)
print("Dimension of null matrix is")
print(dim(N2))
print("Rank of measurement matrix is")
print(rankM2)
print("Dimension of measurement matrix is")
print(dim(M2))

# Make sure the rank and null dimensions are consistent
if (rankM2 + ncol(N2) != ncol(M2)) {
  stop("rank + null dimension not equal to number of columns in M")
}

dtauVect2 <- bd_calc_trapez_weights(tauVect2)
r <- log(1.02)
v <- exp(r * tauVect2)
# v <- v / (sum(v*dtauVect2))
v <- calcExpPdf(tauVect2, r, taumin2, taumax2)
pdf(here::here("analysis/figures/FigS3_standard_basis_exp.pdf"), width = 10, height = 6)
plot(NULL, type = "n", xlim = c(taumin2, taumax2), ylim = c(0, .03), xlab = "Calendar Date [AD]", ylab = "Probability Density")

# Add grey identifiable regions
for (ii in 1:length(invSpanList2)) {
  invReg <- invSpanList2[[ii]]
  if (dplyr::between(invReg$tau_left, taumin2, taumax2) || dplyr::between(invReg$tau_right, taumin2, taumax2)) {
    rect(invReg$tau_left, 0, invReg$tau_right, .03, border = NA, col = "gray80")
  }
}

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- mapply(brewer.pal, ncol(N2), "Set1")
for (ii in 1:ncol(N2)) {
  eta <- N2[, ii]
  adjFact <- .015 / max(abs(N2))
  lines(tauVect2, v + adjFact * eta, col = col_vector[ii], lwd = 3)
}
# Plot the exponential in the two invertible regions only
tau_span1 <- seq(min(tauVect2), tauVect2[which(reg[[1]]$tau_right == tauVect2) - 1], len = 100)
v_span1 <- calcExpPdf(tau_span1, r, taumin2, taumax2)
tau_span2 <- seq(tauVect2[which(reg[[2]]$tau_left == tauVect2) + 1], max(tauVect2), len = 100)
v_span2 <- calcExpPdf(tau_span2, r, taumin2, taumax2)
lines(tau_span1, v_span1, col = "black", lwd = 3)
lines(tau_span2, v_span2, col = "black", lwd = 3)
dev.off()

outlog$`taumin2` <- taumin2
outlog$`taumax2` <- taumax2
outlog$`nullSize2` <- ncol(N2)
outlog$`numRow_M2` <- nrow(M2)
outlog$`numCol_M2` <- ncol(M2)
outlog$`rankM2` <- rankM2
outlog$`numPoints_tau2` <- length(tauVect2)
outlog$`numPoints_phi2` <- length(phiVect2)


print("Key parameter values:")
tibble::tibble(
  Parameter = names(outlog),
  Value = unlist(outlog)
) %T>%
  readr::write_csv(here::here("analysis/logs/SuppC_sb.csv"))
