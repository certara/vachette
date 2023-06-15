
#' Initialize vachette object with required data
#'
#' @param obs.data data.frame; Observed data
#' @param typ.data data.frame; Typical (population) curves
#' @param sim.data data.frame; Simulated (VPC) data
#' @param covariates named character vector; Covariate names with reference values in vachette transformation
#' @param ref.dosenr integer; Dose number to use as the reference dose, corresponding to value in "dosenr" column in input data
#' @param iiv.correction logical; Apply inter-individual variability correction. Default \code{FALSE}
#' @param error.model character; Applied error model, \code{"proportional"} or \code{"additive"}. Default \code{"proportional"}
#' @param model.name character; Optional model name for plot output
#' @param mappings named character vector;  Optional mappings to be included if column names in input \code{data.frame} differ from required column names.
#' See Required Columns section:
#' @section Required columns \code{obs.data}:
#' \itemize{
#'  \item{"ID"}{ - Subject ID}
#'  \item{"x"}{ - Typically time}
#'  \item{"PRED"}{ - Population prediction, required if \code{iiv.correction = TRUE}}
#'  \item{"IPRED"}{ - Individual prediction, required if \code{iiv.correction = TRUE}}
#'  \item{"OBS"}{ - DV}
#'  \item{"dosenr"}{ - Dose number; unique dose number for ID/time point.
#'  }
#' }
#'
#' @section Required columns \code{typ.data}:
#' \itemize{
#'  \item{"ID"}{ - Subject ID}
#'  \item{"x"}{ - Typically time}
#'  \item{"PRED"}{ - Population prediction}
#'  \item{"dosenr"}{ - Dose number; unique dose number for ID/time point}
#' }
#'
#' @section Required columns \code{sim.data}:
#' \itemize{
#'  \item{"ID"}{ - Subject ID}
#'  \item{"x"}{ - Typically time}
#'  \item{"PRED"}{ - Population prediction, required if \code{iiv.correction = TRUE}}
#'  \item{"IPRED"}{ - Individual prediction, required if \code{iiv.correction = TRUE}}
#'  \item{"REP"}{ - Replicate number}
#' }
#'
#' @details If "dosenr" column is missing it will be automatically calculated using the priority of available columns:
#'    \itemize{
#'      \item{"EVID": If available in data, "dosenr" will be calculated using \code{cumsum(EVID==1)}}
#'      \item{"ADDL"/"II": If "ADDL" and "II" are available in data, "dosenr" will be calculated given additional dose number and interval}
#'      \item{"AMT": If only "AMT" column exists in data, "dosenr" will be calculated using \code{cumsum(AMT!=0)}}
#'    }
#'  If \code{sim.data} argument is used, "dosenr" is not required in corresponding \code{sim.data} \code{data.frame}; "dosenr" is extracted from
#'  \code{obs.data} and merged into \code{sim.data} by "ID" and "x" key.
#'
#' @return \code{vachette_data}
#' @export
#'
#' @examples
#' obs <- read.csv(system.file(package = "vachette", "examples", "iv-obs.csv"))
#' typ <- read.csv(system.file(package = "vachette", "examples", "iv-typ.csv"))
#'
#' vd <- vachette_data(
#'   obs.data = obs,
#'   typ.data = typ,
#'   covariates = c(WT=70),
#'   model.name  = "IV"
#'  ) |>
#'  apply_transformations()
#'
#'
vachette_data <-
  function(obs.data,
           typ.data,
           sim.data = NULL,
           covariates,
           ref.dosenr,
           iiv.correction = FALSE,
           error.model = c("proportional", "additive"),
           model.name = NULL,
           mappings = NULL) {

  vachette_data_env <- environment()
  # Column Validation Check
  .validate_columns(mappings, obs.data, typ.data, sim.data, iiv.correction) %>%
    list2env(envir = vachette_data_env)

  # Process/assign covariates
  covariates <- .process_covariates(covariates, obs.data)

  obs.data <- .calculate_dose_number(obs.data, data_type = "obs.data")
  typ.data <- .calculate_dose_number(typ.data, data_type = "typ.data")

  number_of_doses <- unique(obs.data[["dosenr"]])
  if (missing(ref.dosenr)) {
    if (length(number_of_doses) > 1) {
      warning("ref.dosenr argument is missing and more than one dose number found in data,
              setting ref.dosenr to 1", call. = FALSE)
    }
    ref.dosenr <- 1
  }
  # ensure value supplied to ref.dosenr exists inside dosenr column
  if (ref.dosenr %notin% number_of_doses) {
    stop("ref.dosenr value of ", ref.dosenr, " not found in dosenr column in ", data_type, call. = FALSE)
  }

  # Region is the Vachette terminology for the time between two dose administrations
  ref.region <- ref.dosenr
  # "Dummy" observations replicate number
  obs.data$REP <- 1
  error <- match.arg(error.model)
  # assign error flags
  if (error == "proportional") {
    ADD_TR <- FALSE
    PROP_TR <- TRUE
  } else {
    ADD_TR <- TRUE
    PROP_TR <- FALSE
  }

  stopifnot(names(covariates) %in% names(typ.data))

  if (!is.null(sim.data)) {
    VVPC <- TRUE
    stopifnot(names(covariates) %in% names(sim.data))
    #Join dosenr from obs.data by ID and x key if dosenr missing
    if ("dosenr" %notin% colnames(sim.data)) {
      sim.data <- left_join(sim.data, obs.data %>% select(ID, x, dosenr), by = c("ID", "x"))
    }
  } else {
    VVPC <- FALSE
    # "Dummy" vpc simulated observations dataset
    sim.data <- obs.data
  }

  # Retain required columns
  if (iiv.correction) {
    obs_cols_to_select <- c("REP", "ID", "x", "PRED", "IPRED", "OBS", "dosenr")
  } else {
    obs_cols_to_select <- c("REP", "ID", "x", "OBS", "dosenr")
  }
  typ_cols_to_select <- c("ID", "x", "PRED", "dosenr")
  obs.data  <-
    obs.data %>% dplyr::select(dplyr::all_of(obs_cols_to_select), dplyr::all_of(names(covariates)))
  sim.data  <-
    sim.data %>% dplyr::select(dplyr::all_of(obs_cols_to_select), dplyr::all_of(names(covariates)))
  typ.data  <-
    typ.data %>% dplyr::select(dplyr::all_of(typ_cols_to_select), dplyr::all_of(names(covariates)))

  # Extract last observed x (look for max x/time in obs and sim data)
  xstop <- max(obs.data$x,sim.data$x)
  # Define unique covariate combination as new 'COV' column.
  obs.orig <- obs.data %>%
    mutate(COV = paste(!!!syms(names(covariates))))
  sim.orig <- sim.data %>%
    mutate(COV = paste(!!!syms(names(covariates))))
  typ.orig <- typ.data %>%
    mutate(COV = paste(!!!syms(names(covariates))))

  .define_and_enumerate_regions(typ.orig, obs.orig, sim.orig, covariates, ref.dosenr, ref.region) %>%
    list2env(envir = vachette_data_env)

  # Apply IIV Correction
  typ.orig <- typ.orig %>%
    mutate(y=PRED)
  # to apply
  if(iiv.correction) {
    obs.orig$y <- obs.orig$OBS  - obs.orig$IPRED + obs.orig$PRED
    sim.orig$y <- sim.orig$OBS  - sim.orig$IPRED + sim.orig$PRED
  } else {
    obs.orig$y <- obs.orig$OBS
    sim.orig$y <- sim.orig$OBS
  }

  list(
    model.name = model.name,
    covariates = covariates,
    VVPC = VVPC,
    ADD_TR = ADD_TR,
    PROP_TR = PROP_TR,
    obs.orig = obs.orig,
    sim.orig = sim.orig,
    typ.orig = typ.orig,
    tab.ucov = tab.ucov,
    ref.region = ref.region,
    ref.dosenr = ref.dosenr,
    xstop = xstop,
    n.ucov = n.ucov
  ) %>%
    structure(class = "vachette_data")
}

#' @export
update.vachette_data <- function(vachette_data, ...) {
  args <- list(...)
  for (i in names(args)) {
    vachette_data[[i]] <- args[[i, exact=TRUE]]
  }
  vachette_data
}

.define_and_enumerate_regions <- function(typ.orig, obs.orig, sim.orig, covariates, ref.dosenr, ref.region) {
  # ----- Define regions and provide a number to all combined covariate effects ----
  # Region number
  typ.orig$region  <- NA
  obs.orig$region    <- NA
  sim.orig$region    <- NA

  # Region type (open/closed)
  typ.orig$region.type  <- NA
  obs.orig$region.type    <- NA
  sim.orig$region.type    <- NA

  # Covariate combinations number
  typ.orig$ucov  <- NA
  obs.orig$ucov    <- NA
  sim.orig$ucov    <- NA

  n.ucov    <- 0
  tab.ucov  <- NULL # Table with ucov properties

  # Get all covariate combinations
  comb.ucov <- typ.orig %>% tidyr::expand(!!!syms(names(covariates)))

  # Collect unique covariate/dose combinations and define region.type
  # Add to info to data frames
  if(length(covariates)==1)
  {
    for(i in c(1:dim(comb.ucov)[1])) #e.g., unique values of covariates nrow(comb.ucov)
    {
      # Get number of doses for each covariate combination:
      z <- typ.orig %>% filter(!!sym(names(covariates)) == comb.ucov[names(covariates)][[1]][i])
      ndose <- length(unique(z$dosenr))
      for(idose in c(1:ndose))
      {
        if (idose == ndose) region.type <- 'open'   # Last dose
        if (idose <  ndose) region.type <- 'closed' # Any dose before last dose

        # New unique combination
        n.ucov <- n.ucov + 1
        add <- data.frame(ucov = n.ucov) %>%
          mutate(!!names(covariates) := comb.ucov[names(covariates)][[1]][i],
                 region = idose,
                 region.type = region.type)
        tab.ucov <- rbind(tab.ucov, add)

        # Selection of new unique combination
        sel.typ  <- typ.orig[[names(covariates)]] == comb.ucov[names(covariates)][[1]][i] &
          typ.orig$dosenr == idose
        sel.obs  <- obs.orig[[names(covariates)]] == comb.ucov[names(covariates)][[1]][i] &
          obs.orig$dosenr == idose
        sel.sim  <- sim.orig[[names(covariates)]] == comb.ucov[names(covariates)][[1]][i] &
          sim.orig$dosenr == idose

        # Add new unique combination info
        typ.orig$region[sel.typ]   <- idose
        obs.orig$region[sel.obs]     <- idose
        sim.orig$region[sel.sim]     <- idose

        typ.orig$ucov[sel.typ]   <- n.ucov
        obs.orig$ucov[sel.obs]     <- n.ucov
        sim.orig$ucov[sel.sim]     <- n.ucov

        typ.orig$region.type[sel.typ]   <- region.type
        obs.orig$region.type[sel.obs]     <- region.type
        sim.orig$region.type[sel.sim]     <- region.type
      }
    }
  }

  # To do, improve to account for n covariates
  if(length(covariates) > 1)
  {
    vars <- names(covariates)

    for(i in c(1:dim(comb.ucov)[1]))
    {
      filter_query <- paste0(vars, "==", "comb.ucov$", vars, "[i]", collapse = " & ")
      # Get number of doses for covariate combination:
      z <- typ.orig %>% filter(eval(parse(text = filter_query)))

      ndose <- length(unique(z$dosenr))

      for(idose in c(1:ndose))
      {
        if (idose == ndose) region.type <- 'open'
        if (idose <  ndose) region.type <- 'closed'

        # New unique combination
        n.ucov <- n.ucov + 1

        add <- data.frame(ucov = n.ucov)
        vals <- comb.ucov %>% slice(i)
        add[vars] <- vals
        add <- mutate(add, region = idose,
                      region.type = region.type)

        tab.ucov <- rbind(tab.ucov, add)
        filter_query_typ <- paste0(
          paste0("typ.orig$", vars, "==", "comb.ucov$", vars, "[i]", collapse = " & "),
          " & typ.orig$dosenr == idose"
        )
        filter_query_obs <- paste0(
          paste0("obs.orig$", vars, "==", "comb.ucov$", vars, "[i]", collapse = " & "),
          " & obs.orig$dosenr == idose"
        )
        filter_query_sim <- paste0(
          paste0("sim.orig$", vars, "==", "comb.ucov$", vars, "[i]", collapse = " & "),
          " & sim.orig$dosenr == idose"
        )
        # Selection of new unique combination
        sel.typ <- eval(parse(text=filter_query_typ))
        sel.obs <- eval(parse(text=filter_query_obs))
        sel.sim <- eval(parse(text=filter_query_sim))
        # sel.typ  <- typ.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        #   typ.orig$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        #   typ.orig$dosenr == idose
        # sel.obs  <- obs.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        #   obs.orig$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        #   obs.orig$dosenr == idose
        # sel.sim  <- sim.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        #   sim.orig$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        #   sim.orig$dosenr == idose

        # Add new unique combination info
        typ.orig$region[sel.typ]   <- idose
        obs.orig$region[sel.obs]     <- idose
        sim.orig$region[sel.sim]     <- idose

        typ.orig$ucov[sel.typ]   <- n.ucov
        obs.orig$ucov[sel.obs]     <- n.ucov
        sim.orig$ucov[sel.sim]     <- n.ucov

        typ.orig$region.type[sel.typ]   <- region.type
        obs.orig$region.type[sel.obs]     <- region.type
        sim.orig$region.type[sel.sim]     <- region.type
      }
    }
  }

  tab.ucov$ref <- "No"

  for(i.ucov in c(1:dim(tab.ucov)[1]))
  {
    # if(length(covariates)==1)
    #   if(tab.ucov[i.ucov, covariates] == ref.cov1 && tab.ucov[i.ucov, "region"] == ref.region)
    #     tab.ucov[i.ucov, "ref"] <- "Yes"
    # # where both covs equal ref values in table, set yes
    # if(length(covariates)>1)
      condition <- paste0(
      paste0("as.character(tab.ucov$", names(covariates), "[i.ucov])", "==", "'", covariates, "'", collapse = " && "),
      " & tab.ucov$region[i.ucov] == ref.region"
    )
      if(eval(parse(text = condition)))
        tab.ucov$ref[i.ucov] <- "Yes"
  }

  # Add reference flag to typ.orig
  typ.orig <- typ.orig %>%
    dplyr::left_join(tab.ucov[,c('ucov','ref')],by='ucov')

  return(
    list(
      typ.orig = typ.orig,
      obs.orig = obs.orig,
      sim.orig = sim.orig,
      tab.ucov = tab.ucov,
      n.ucov = n.ucov
    )
  )
}


#' Apply vachette transformations
#'
#' @param vachette_data object of class \code{vachette_data}
#' @param lm.refine logical; Refine landmark x-position. Default \code{FALSE}
#' @param tol.end numeric; Relative tolerance to determine last x open end reference
#' @param tol.noise numeric; Relative tolerance for landmark determination typical curves
#' @param step.x.factor numeric; x-axis extension factor to search for last x, i.e., to determine where close enough to asymptote
#' @param ngrid.fit numeric; number of grid points in last query segment for matching last reference segment
#' @param window integer; size (gridpoints) of Savitzky Golay smoothing window for landmark position determination
#' @param window.d1.refine integer; size (gridpoints) of Savitzky Golay smoothing window for refinement of first derivative landmark position
#' @param window.d2.refine integer; size (gridpoints) of Savitzky Golay smoothing window for refinement of second derivative landmark position.
#'
#' @name apply_transformations
#' @export
apply_transformations <- function(vachette_data, ...) UseMethod("apply_transformations")

#' @rdname apply_transformations
#' @export
apply_transformations.vachette_data <-
  function(vachette_data,
           lm.refine = FALSE,
           tol.end = 0.001,
           tol.noise = 1e-8,
           step.x.factor = 1.5,
           ngrid.fit = 100,
           window = 17,     # Savitzky Golay smoothing - initial landmarks
           window.d1.refine = 7,     # Savitzky Golay smoothing - refine landmarks - first derivative
           window.d2.refine = 5) {

  stopifnot(inherits(vachette_data, "vachette_data"))

  vachette_transformed_data <- .calculate_transformations(vachette_data,
                                                         lm.refine,
                                                         tol.end,
                                                         tol.noise,
                                                         step.x.factor,
                                                         ngrid.fit,
                                                         window,
                                                         window.d1.refine,
                                                         window.d2.refine)
  if (vachette_data$VVPC) {
    vachette_transformed_data_sim <- .calculate_transformations(vachette_data,
                                                           lm.refine,
                                                           tol.end,
                                                           tol.noise,
                                                           step.x.factor,
                                                           ngrid.fit,
                                                           window,
                                                           window.d1.refine,
                                                           window.d2.refine,
                                                           run_sim = TRUE)
    sim.all <- vachette_transformed_data_sim$obs.all
  } else {
    sim.all <- NULL
  }

    update(vachette_data,
           obs.all = vachette_transformed_data$obs.all,
           sim.all = sim.all,
           curves.all = vachette_transformed_data$curves.all,
           curves.scaled.all = vachette_transformed_data$curves.scaled.all,
           ref.extensions.all = vachette_transformed_data$ref.extensions.all,
           lm.all = vachette_transformed_data$lm.all,
           nseg = vachette_transformed_data$nseg)
  }

