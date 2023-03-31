test_that("user mappings provided give correct errors", {

  indivsam.obs <- read.csv(system.file(package = "vachette", "examples", "oral-two-cov-obs.csv"))
  output.typ <- read.csv(system.file(package = "vachette", "examples", "oral-two-cov-typ.csv"))

  colnames(indivsam.obs) <- tolower(colnames(indivsam.obs))
  colnames(output.typ) <- tolower(colnames(output.typ))

  #Leave out  OBS
  vachette_data(
    indivsam.obs,
    output.typ,
    vachette.covs =  c("wt" = 70, "age" = 30),
    ref.dose = 1,
    mappings = c(ID="id", IPRED="ipred", PRED="pred")
  ) %>%
    testthat::expect_error()

  #Typo in column name 'ip' instead of 'id'
  vachette_data(
    indivsam.obs,
    output.typ,
    vachette.covs =  c("wt" = 70, "age" = 30),
    ref.dose = 1,
    mappings = c(ID="ip", IPRED="ipred", PRED="pred", OBS = "obs")
  ) %>%
    testthat::expect_error()

  colnames(indivsam.obs) <- toupper(colnames(indivsam.obs))
  colnames(output.typ) <- toupper(colnames(output.typ))

  #Missing mappings
  vachette_data(
    indivsam.obs,
    output.typ,
    vachette.covs =  c("WT" = 70, "AGE" = 30),
    ref.dose = 1
  ) %>%
    testthat::expect_error()

  #incomplete mappings
  vachette_data(
    indivsam.obs,
    output.typ,
    vachette.covs =  c("WT" = 70, "AGE" = 30),
    ref.dose = 1,
    mappings = c(dosenr = "DOSENR")
  ) %>%
    testthat::expect_error()

})

