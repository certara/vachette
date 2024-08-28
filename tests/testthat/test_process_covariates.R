test_that("covariates are correctly extracted if missing reference values", {

  #DF
  obs.data <- read.csv(system.file(package = "vachette", "examples", "iv-obs.csv"))
  covariates <- vachette:::.process_covariates(covariates = c("WT"), obs.data)
  testthat::expect_equal(covariates, c("WT" = "30"))

  #Test Tibble
  mtcars$cyl <- paste0("cyl_", mtcars$cyl)
  covariates <- vachette:::.process_covariates(covariates = c("mpg", "cyl"), mtcars)
  testthat::expect_equal(covariates, c("mpg" = "19.2", "cyl" = "cyl_8"))

  covariates <- vachette:::.process_covariates(covariates = c("mpg" = 18.7, "cyl" = "cyl_6"), mtcars)
  testthat::expect_equal(covariates, c("mpg" = "18.7", "cyl" = "cyl_6"))

  covariates <- vachette:::.process_covariates(covariates = c("cyl", "hp" = 123), mtcars)
  testthat::expect_equal(covariates, c("cyl" = "cyl_8", "hp" = "123"))
})

