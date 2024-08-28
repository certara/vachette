test_that("user mappings provided give correct errors", {

  obs.data <- read.csv(system.file(package = "vachette", "examples", "iv-obs.csv"))
  typ.data <- read.csv(system.file(package = "vachette", "examples", "iv-typ-minmax.csv"))

  colnames(obs.data) <- tolower(colnames(obs.data))
  colnames(typ.data) <- tolower(colnames(typ.data))


  #Typo in column name 'ip' instead of 'id'
  testthat::expect_error(
    vachette_data(
      obs.data,
      typ.data,
      covariates =  c("wt" = 70),
      ref.dosenr = 1,
      mappings = c(
        x = "time",
        ID = "ip",
        PRED = "pred",
        OBS = "dv"
      )
    ),
    regexp = "ip column\\(s\\) provided in mappings argument are not found in obs.data."
  )

  colnames(obs.data) <- toupper(colnames(obs.data))
  colnames(typ.data) <- toupper(colnames(typ.data))

  #Missing mappings
  testthat::expect_error(
    vachette_data(
    obs.data,
    typ.data,
    covariates =  c("WT" = 70),
    ref.dosenr = 1
  ),
  regexp =  "OBS x column\\(s\\) not found in obs.data")


})

