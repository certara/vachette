obs <- read.csv(system.file(package = "vachette", "examples", "iv-obs.csv"))
typ <- read.csv(system.file(package = "vachette", "examples", "iv-typ.csv"))
sim <- read.csv(system.file(package = "vachette", "examples", "iv-sim.csv"))

vd <-
  vachette_data(
    obs.data = obs,
    typ.data = typ,
    sim.data = sim,
    covariates = c(WT=70),
    mappings = c(x = "time",
                 OBS = "DV"),
    model.name  = "Intravenous (i.v)"
  ) |>
  apply_transformations(log_file = "log.txt")
