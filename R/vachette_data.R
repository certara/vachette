
#' Initialize vachette object with required data
#'
#' @param indivsam.obs Observed \code{data.frame}
#' @param output.typ Typ \code{data.frame}
#' @param indivsam.vpc VPC \code{data.frame}
#' @param vachette.covs chracter vector of covariate names
#'
#' @return `vachette_data`
#' @export
#'
#' @examples
#' indivsam.obs <- read.csv(system.file(package = "vachette", "examples", "iv-obs.csv"))
#' output.typ <- read.csv(system.file(package = "vachette", "examples", "iv-typ.csv"))
#' indivsam.vpc <- read.csv(system.file(package = "vachette", "examples", "iv-vpc.csv"))
#'
#'
vachette_data <- function(indivsam.obs, output.typ, indivsam.vpc = NULL, vachette.covs) {

  # *Enhancement support tidy evaluation of vachette.covs

  # data QC checks to perform on obs, typ, vpc data
  stopifnot(vachette.covs %in% names(indivsam.obs))
  # "Dummy" observations replicate number
  indivsam.obs$isim      <- 1

  stopifnot(vachette.covs %in% names(output.typ))

  if (!is.null(indivsam.vpc)) {
    VVPC <- TRUE
    stopifnot(vachette.covs %in% names(indivsam.vpc))
  } else {
    VVPC <- FALSE
    # "Dummy" vpc simulated observations dataset
    indivsam.vpc <- indivsam.obs
  }


  # Retain required columns
  indivsam.obs  <-
    indivsam.obs %>% dplyr::select(isim, ID, x, PRED, IPRED, OBS, dplyr::all_of(vachette.covs), dosenr)
  indivsam.vpc  <-
    indivsam.vpc %>% dplyr::select(isim, ID, x, PRED, IPRED, OBS, dplyr::all_of(vachette.covs), dosenr)

  # Extract last observed x (look for max x/time in obs and sim data)
  xstop <- max(indivsam.obs$x,indivsam.vpc$x)
  # Define unique covariate combination as new 'COV' column.
  obs.orig <- indivsam.obs %>%
    mutate(COV = paste(!!!syms(vachette.covs)))
  sim.orig <- indivsam.vpc %>%
    mutate(COV = paste(!!!syms(vachette.covs)))
  output.typ <- output.typ %>%
    mutate(COV = paste(!!!syms(vachette.covs)))

  vachette_data_env <- environment()

  .define_and_enumerate_regions(output.typ, obs.orig, sim.orig, vachette.covs) %>%
    list2env(envir = vachette_data_env)


  list(
    indivsam.obs = indivsam.obs,
    output.typ = output.typ,
    indivsam.vpc = indivsam.vpc,
    vachette.covs = vachette.covs,
    VVPC = VVPC,
    output.typ = output.typ,
    obs.orig = obs.orig,
    sim.orig = sim.orig,
    tab.ucov = tab.ucov
  ) %>%
    structure(class = "vachette_data")
}


.define_and_enumerate_regions <- function(output.typ, obs.orig, sim.orig, vachette.covs) {
  # ----- Define regions and provide a number to all combined covariate effects ----
  # Region number
  output.typ$region  <- NA
  obs.orig$region    <- NA
  sim.orig$region    <- NA

  # Region type (open/closed)
  output.typ$region.type  <- NA
  obs.orig$region.type    <- NA
  sim.orig$region.type    <- NA

  # Covariate combinations number
  output.typ$ucov  <- NA
  obs.orig$ucov    <- NA
  sim.orig$ucov    <- NA

  n.ucov    <- 0
  tab.ucov  <- NULL # Table with ucov properties

  # Get all covariate combinations
  comb.ucov <- output.typ %>% tidyr::expand(!!!syms(vachette.covs))

  # Collect unique covariate/dose combinations and define region.type
  # Add to info to data frames
  if(length(vachette.covs)==1)
  {
    for(i in c(1:dim(comb.ucov)[1])) #e.g., unique values of covariates nrow(comb.ucov)
    {
      # Get number of doses for each covariate combination:
      z <- output.typ %>% filter(!!sym(vachette.covs) == comb.ucov[vachette.covs][[1]][i])
      ndose <- length(unique(z$dosenr))
      for(idose in c(1:ndose))
      {
        if (idose == ndose) region.type <- 'open'   # Last dose
        if (idose <  ndose) region.type <- 'closed' # Any dose before last dose

        # New unique combination
        n.ucov <- n.ucov + 1
        add <- data.frame(ucov = n.ucov) %>%
          mutate(!!vachette.covs := comb.ucov[vachette.covs][[1]][i],
                 region = idose,
                 region.type = region.type)
        tab.ucov <- rbind(tab.ucov, add)

        # Selection of new unique combination
        sel.typ  <- output.typ[[vachette.covs]] == comb.ucov[vachette.covs][[1]][i] &
          output.typ$dosenr == idose
        sel.obs  <- obs.orig[[vachette.covs]] == comb.ucov[vachette.covs][[1]][i] &
          obs.orig$dosenr == idose
        sel.sim  <- sim.orig[[vachette.covs]] == comb.ucov[vachette.covs][[1]][i] &
          sim.orig$dosenr == idose

        # Add new unique combination info
        output.typ$region[sel.typ]   <- idose
        obs.orig$region[sel.obs]     <- idose
        sim.orig$region[sel.sim]     <- idose

        output.typ$ucov[sel.typ]   <- n.ucov
        obs.orig$ucov[sel.obs]     <- n.ucov
        sim.orig$ucov[sel.sim]     <- n.ucov

        output.typ$region.type[sel.typ]   <- region.type
        obs.orig$region.type[sel.obs]     <- region.type
        sim.orig$region.type[sel.sim]     <- region.type
      }
    }
  }

  # To do, improve to account for n vachette.covs
  # if(length(vachette.covs)==2)
  # {
  #   for(i in c(1:dim(comb.ucov)[1]))
  #   {
  #     # Get number of doses for covariate combination:
  #     z <- output.typ %>% filter(vachette.cov1 == comb.ucov$vachette.cov1[i] &
  #                                  vachette.cov2 == comb.ucov$vachette.cov2[i])
  #     ndose <- length(unique(z$dosenr))
  #
  #     for(idose in c(1:ndose))
  #     {
  #       if (idose == ndose) region.type <- 'open'
  #       if (idose <  ndose) region.type <- 'closed'
  #
  #       # New unique combination
  #       n.ucov <- n.ucov + 1
  #       add <- data.frame(ucov = n.ucov,
  #                         vachette.cov1 = comb.ucov$vachette.cov1[i],
  #                         vachette.cov2 = comb.ucov$vachette.cov2[i],
  #                         region = idose,
  #                         region.type = region.type)
  #       tab.ucov <- rbind(tab.ucov, add)
  #
  #       # Selection of new unique combination
  #       sel.typ  <- output.typ$vachette.cov1 == comb.ucov$vachette.cov1[i] &
  #         output.typ$vachette.cov2 == comb.ucov$vachette.cov2[i] &
  #         output.typ$dosenr == idose
  #       sel.obs  <- obs.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
  #         obs.orig$vachette.cov2 == comb.ucov$vachette.cov2[i] &
  #         obs.orig$dosenr == idose
  #       sel.sim  <- sim.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
  #         sim.orig$vachette.cov2 == comb.ucov$vachette.cov2[i] &
  #         sim.orig$dosenr == idose
  #
  #       # Add new unique combination info
  #       output.typ$region[sel.typ]   <- idose
  #       obs.orig$region[sel.obs]     <- idose
  #       sim.orig$region[sel.sim]     <- idose
  #
  #       output.typ$ucov[sel.typ]   <- n.ucov
  #       obs.orig$ucov[sel.obs]     <- n.ucov
  #       sim.orig$ucov[sel.sim]     <- n.ucov
  #
  #       output.typ$region.type[sel.typ]   <- region.type
  #       obs.orig$region.type[sel.obs]     <- region.type
  #       sim.orig$region.type[sel.sim]     <- region.type
  #     }
  #   }
  # }
  return(
    list(
      output.typ = output.typ,
      obs.orig = obs.orig,
      sim.orig = sim.orig,
      tab.ucov = tab.ucov
    )
  )
}

