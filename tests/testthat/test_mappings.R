test_that("user mappings provided give correct errors", {

  obs.data <- read.csv(system.file(package = "vachette", "examples", "oral-two-cov-obs.csv"))
  typ.data <- read.csv(system.file(package = "vachette", "examples", "oral-two-cov-typ.csv"))

  colnames(obs.data) <- tolower(colnames(obs.data))
  colnames(typ.data) <- tolower(colnames(typ.data))

  #Leave out  OBS
  vachette_data(
    obs.data,
    typ.data,
    covariates =  c("wt" = 70, "age" = 30),
    ref.dosenr = 1,
    mappings = c(ID="id", IPRED="ipred", PRED="pred")
  ) %>%
    testthat::expect_error()

  #Typo in column name 'ip' instead of 'id'
  vachette_data(
    obs.data,
    typ.data,
    covariates =  c("wt" = 70, "age" = 30),
    ref.dosenr = 1,
    mappings = c(ID="ip", IPRED="ipred", PRED="pred", OBS = "obs")
  ) %>%
    testthat::expect_error()

  colnames(obs.data) <- toupper(colnames(obs.data))
  colnames(typ.data) <- toupper(colnames(typ.data))

  #Missing mappings
  vachette_data(
    obs.data,
    typ.data,
    covariates =  c("WT" = 70, "AGE" = 30),
    ref.dosenr = 1
  ) %>%
    testthat::expect_error()

  #incomplete mappings
  vachette_data(
    obs.data,
    typ.data,
    covariates =  c("WT" = 70, "AGE" = 30),
    ref.dosenr = 1,
    mappings = c(dosenr = "DOSENR")
  ) %>%
    testthat::expect_error()

})

