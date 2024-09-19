test_that("vachette transformations are accurate for iv", {
  skip_on_os(c("mac", "linux"))
  set.seed(121)

  vd <- suppressMessages(
    vachette_data(
      obs.data = read.csv(system.file(package = "vachette", "examples", "iv-obs.csv")),
      typ.data = read.csv(
        system.file(package = "vachette", "examples", "iv-typ-minmax.csv")
      ),
      covariates = c(WT = 70),
      mappings = c(OBS = "DV",
                   x = "time"),
      ref.dosenr = 1,
      model.name = "iv-example"
    )
  )

  log_file <- tempfile()
  vd <- apply_transformations(vd, log_file = log_file)

  # test output against reference
  # ref_log <- readLines(log_file)
  # writeLines(ref_log, "inst/test_files/iv_log.txt" )
  ref_log <-
    readLines(system.file(package = "vachette", "test_files", "iv_log.txt"))
  test_log <- readLines(log_file)
  unlink(log_file)

  # test tranformed data against reference
  #write.csv(vd$obs.all, "inst/test_files/iv_obs_all.csv", row.names = FALSE)
  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "iv_obs_all.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_equal(ref_log, test_log)
  testthat::expect_equal(vd$obs.all, tidyr::as_tibble(reference))
  testthat::expect_true(ggplot2::is.ggplot(p.vachette(vd)))
})


test_that("vachette transformations are accurate for sigmoid", {
  skip_on_os(c("mac", "linux"))
  set.seed(121)

  vd <- suppressMessages(suppressWarnings(
    vachette_data(
      obs.data = read.csv(
        system.file(package = "vachette", "examples", "sigmoid-obs.csv")
      ),
      typ.data = read.csv(
        system.file(package = "vachette", "examples", "sigmoid-typ-minmax.csv")
      ),
      covariates = c(WT = 70),
      mappings = c(OBS = "DV",
                   x = "bmx"),
      ref.dosenr = 1,
      log.x = TRUE
    )
  ))

  vd <-
    suppressMessages(
      apply_transformations(
        vd,
        asymptote_right = TRUE,
        asymptote_left = TRUE,
        zero_asymptote_right = FALSE,
        zero_asymptote_left = TRUE
      )
    )

  #write.csv(vd$obs.all, "inst/test_files/sigmoid_obs_all.csv", row.names = FALSE)

  reference <-
    read.csv(
      system.file(package = "vachette", "test_files", "sigmoid_obs_all.csv"),
      colClasses = c(COV = "character", ucov = "double"),
      numerals = "no.loss"
    )

  testthat::expect_equal(vd$obs.all, tidyr::as_tibble(reference))

})

