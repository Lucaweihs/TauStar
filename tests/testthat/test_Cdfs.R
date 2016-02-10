library(TauStar)
context("Testing the CDFs.")

test_that("pHoeffInd function returns correct values", {
  # Comparing against the values given by the BKR paper
  bkrToHoeff = function(x) {
    return(as.numeric(pHoeffInd(72/pi^4 * x - 1, 10^-6)))
  }

  expect_equal(.00000, bkrToHoeff(.3), tolerance = 10^-5)
  expect_true(abs(.04867 - bkrToHoeff(.6)) <= 10^-5)
  expect_equal(.29652, bkrToHoeff(.9), tolerance = 10^-5)
  expect_equal(.54354, bkrToHoeff(1.2), tolerance = 10^-5)
  expect_equal(.70763, bkrToHoeff(1.5), tolerance = 10^-5)
  expect_equal(.80922, bkrToHoeff(1.8), tolerance = 10^-5)
  expect_equal(.86406, bkrToHoeff(2.05), tolerance = 10^-5)
  expect_equal(.98546, bkrToHoeff(3.9), tolerance = 10^-5)
  expect_equal(.98965, bkrToHoeff(4.2), tolerance = 10^-5)
  expect_equal(.99261, bkrToHoeff(4.5), tolerance = 10^-5)
  expect_equal(.99471, bkrToHoeff(4.8), tolerance = 10^-5)
  expect_equal(.99994, bkrToHoeff(9), tolerance = 10^-5)

  expect_equal(.95 - bkrToHoeff(2.844), tolerance = 10^-5)
  expect_equal(.99 - bkrToHoeff(4.230), tolerance = 10^-5)
  expect_equal(.995 - bkrToHoeff(4.851), tolerance = 10^-5)

  # Just some sanity checks
  expect_equal(pHoeffInd(-2), 0)
  expect_equal(pHoeffInd(10), 1, tolerance = 10^-6)
})

test_that("dHoeffInd function returns correct values", {
  # Comparing against the values given by the BKR paper
  xVals = c(-0.8256510, -0.7465366, -0.5883078, -0.2718503, 0.3610648,
            1.6268950)
  empTruth = c(0.0004382527, 0.0181979401, 0.5434361861, 1.1925501925,
                     0.3440154196, 0.0349563799)

  for (i in 1:length(xVals)) {
    asymVal = dHoeffInd(xVals[i])
    absDiff = abs(empTruth[i] - asymVal)
    expect_true(abs(empTruth[i] - asymVal) <= 10^-3)
  }
})

test_that("pDisHoeffInd function returns correct values", {
})

