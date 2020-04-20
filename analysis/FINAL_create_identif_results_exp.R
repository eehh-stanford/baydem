library(baydem)

source(here::here("analysis/bayesian_radiocarbon_functions.R"))

set.seed(183450) # from random.org between 1 and 1,000,000

# Load the calibration curve and call the code to assess equifinality
calibDf <- bd_load_calib_curve("intcal13")
tau_curve <- 1950 - calibDf$yearBP # AD
equiList <- bd_calc_calib_curve_equif_dates(calibDf)
phi_curve <- bd_calc_calib_curve_frac_modern(calibDf)
equifData <- bd_assess_calib_curve_equif(calibDf)
canInvert <- equifData$canInvert
invSpanList <- equifData$invSpanList


taumin2 <- 970
taumax2 <- 1035
measError2 <- 0.00001
ind2 <- (tau_curve >= taumin2) & (tau_curve <= taumax2)
phiMin2 <- min(phi_curve[ind2])
phiMax2 <- max(phi_curve[ind2])
phiVect2 <- unique(phi_curve[ind2])
tauVect2 <- tau_curve[ind2]

for (kk in 1:length(equiList)) {
  equiEntry <- equiList[[kk]]
  if (equiEntry$indBase %in% which(ind2)) {
    tauVect2 <- unique(c(tauVect2, equiEntry$tauBase))
    tauVect2 <- unique(c(tauVect2, equiEntry$tauEqui))
  }
}

tauVect2 <- sort(tauVect2)

# M2 <- calcMeasMatrix2(tauVect2,phiVect2,rep(measError2,length(phiVect2)),calibDf,T,c(phiMin2,phiMax2))
# N2 <- MASS::Null(t(M2))

pdf(here::here("analysis/figures/FigS1_exp_example.pdf"), width = 20, height = 18)
# There are three rows of plots. The top row has a single, long plot. The next
# two rows rows have three plots each.
layout(matrix(c(1, 1, 1, 2, 3, 4, 5, 6, 7), 3, 3, byrow = TRUE))
taumina <- 700
taumaxa <- 950
tauminb <- 800
taumaxb <- 850
tauminc <- 900
taumaxc <- 950
col0 <- "indianred1"
measError1 <- 1e-2
measError2 <- 1e-3

# Save named pairs of key information for reporting in the manuscript
outlog <- list()
outlog$`taumina` <- taumina
outlog$`taumaxa` <- taumaxa
outlog$`tauminb` <- tauminb
outlog$`taumaxb` <- taumaxb
outlog$`tauminc` <- tauminc
outlog$`taumaxc` <- taumaxc
outlog$`measError1` <- measError1
outlog$`measError2` <- measError2

kappa <- 8033 # Reference decay rate of carbon-14

# Calculate and save the uncertainty in uncalibrated calendar years for the two
# measurement error settings
measErrorYears1_700 <- kappa * measError1 / (exp(-(1950 - 700) / kappa))
measErrorYears2_700 <- kappa * measError2 / (exp(-(1950 - 700) / kappa))
measErrorYears1_950 <- kappa * measError1 / (exp(-(1950 - 950) / kappa))
measErrorYears2_950 <- kappa * measError2 / (exp(-(1950 - 950) / kappa))

outlog$`measErrorYears1_700` <- measErrorYears1_700
outlog$`measErrorYears2_700` <- measErrorYears2_700
outlog$`measErrorYears1_950` <- measErrorYears1_950
outlog$`measErrorYears2_950` <- measErrorYears2_950

# The vector of growth rates plot, -4% to 4% per annum
rVect <- log(1 + seq(-.04, .04, by = .01))

# Show the calibration curve in the top row
bd_vis_calib_curve(700, 950, calibDf, xlab = "Calendar Date [AD]", ylab = "Fraction Modern", invertCol = "gray80")
# Add red rectangles to the calibration curve plot to mark the time periods
# used to generate the fraction modern probability densities
rect(taumina, .865, taumaxa, .866, border = NA, col = col0)
text(mean(c(taumina, taumaxa)), .867, "Span 1", cex = 2)
rect(tauminb, .855, taumaxb, .856, border = NA, col = col0)
text(mean(c(tauminb, taumaxb)), .8535, "Span 2", cex = 2)
rect(tauminc, .855, taumaxc, .856, border = NA, col = col0)
text(mean(c(tauminc, taumaxc)), .8535, "Span 3", cex = 2)

# Add plots for the six cases (three timespans by two measurement error
# settings).
#
# If there is an identifiability problem for any cases, an error is thrown by
# visExpEquif.
visExpEquif(rVect, taumina, taumaxa, calibDf, measError1, 1000)
visExpEquif(rVect, tauminb, taumaxb, calibDf, measError1, 1000)
visExpEquif(rVect, tauminc, taumaxc, calibDf, measError1, 1000)
visExpEquif(rVect, taumina, taumaxa, calibDf, measError2, 1000)
visExpEquif(rVect, tauminb, taumaxb, calibDf, measError2, 1000)
visExpEquif(rVect, tauminc, taumaxc, calibDf, measError2, 1000)
dev.off()

tibble::tibble(
  Parameter = names(outlog),
  Value = unlist(outlog)
) %T>%
  readr::write_csv(here::here("analysis/logs/SuppB_exp.csv"))
