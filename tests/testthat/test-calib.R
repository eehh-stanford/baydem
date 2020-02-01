# Run tests related to the calibration dataframe, including functions related
# to equifinality / identifiability.

# Check that loading the intcal13 data does not raise an error
expect_error(
  calibDf <- bd_load_calib_curve("intcal13")
,NA
)
# Check that calling bd_assess_calib_curve_equif does not raise an error and
# also make sure that the return is a list with invSpanList and canIvert. This
# also checks the helper function bd_phi2tau called by
# bd_assess_calib_curve_equif.
expect_error(
  equifResult <- bd_assess_calib_curve_equif(calibDf)
,NA
)

expect_equal(
  names(equifResult), c('invSpanList','canInvert')
)

# Check that calling bd_assess_calib_curve_equif does not raise an error and
# also make check that bd_assess_calib_curve_equif works if equiList is input
expect_error(
  equifDates <- bd_calc_calib_curve_equif_dates(calibDf)
,NA
)

expect_error(
  equifResult2 <- bd_assess_calib_curve_equif(calibDf,equifDates)
,NA
)

expect_equal(
  names(equifResult2), c('invSpanList','canInvert')
)

# Check that calling bd_calc_calib_curve_frac_modern does not raise an error.
# Also check that if tau is not input the function returns exp(-calibDf$uncalYearBP / 8033)
expect_error(
  phi <- bd_calc_calib_curve_frac_modern(calibDf)
,NA
)

expect_equal(
  phi, exp(-calibDf$uncalYearBP / 8033)
)

# Check that calling bd_calc_calib_curve_frac_modern does not raise an error
# when tau is input. Also check that the output length is at least the same
# length as the input
tau2 <- c(600,602,805.89)
expect_error(
  phi2 <- bd_calc_calib_curve_frac_modern(calibDf,tau2)
,NA
)

expect_equal(
  length(phi2), length(tau2)
)

# Check that calling bd_vis_calib_curve does not raise an error
expect_error(
  bd_vis_calib_curve(600,700,calibDf)
,NA
)
