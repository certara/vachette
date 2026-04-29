test_that(".warn_if_sim_obs_join_empty is silent when join has rows", {
  sim_join <- data.frame(REP = 1L, ID = 1L, x = 1, y = 1)
  obs_join <- data.frame(ID = 1L, x = 1)
  out <- dplyr::inner_join(sim_join, obs_join, by = c("ID", "x"))
  expect_equal(nrow(out), 1L)
  expect_silent(vachette:::.warn_if_sim_obs_join_empty(out, sim_join, obs_join))
})

test_that(".warn_if_sim_obs_join_empty warns on near-equal x (floating-point)", {
  sim_join <- data.frame(REP = 1L, ID = 1L, x = 1 + .Machine$double.eps, y = 1)
  obs_join <- data.frame(ID = 1L, x = 1)
  out <- dplyr::inner_join(sim_join, obs_join, by = c("ID", "x"))
  expect_equal(nrow(out), 0L)
  expect_warning(
    vachette:::.warn_if_sim_obs_join_empty(out, sim_join, obs_join),
    regexp = "returned 0 rows|floating-point precision"
  )
})

test_that(".warn_if_sim_obs_join_empty warns when IDs do not overlap", {
  sim_join <- data.frame(REP = 1L, ID = 1L, x = 1, y = 1)
  obs_join <- data.frame(ID = 99L, x = 1)
  out <- dplyr::inner_join(sim_join, obs_join, by = c("ID", "x"))
  expect_equal(nrow(out), 0L)
  expect_warning(
    vachette:::.warn_if_sim_obs_join_empty(out, sim_join, obs_join),
    regexp = "No overlapping"
  )
})
