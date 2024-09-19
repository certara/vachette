
#' Initialize vachette object with required data
#'
#' @param obs.data data.frame; Observed data
#' @param typ.data data.frame; Typical (population) curves
#' @param sim.data data.frame; Simulated (VPC) data
#' @param covariates named character vector; Covariate names with reference values in vachette transformation
#' @param ref.dosenr integer; Dose number to use as the reference dose, corresponding to value in "dosenr" column in input data
#' @param iiv.correction logical; Apply inter-individual variability correction. Default \code{FALSE}
#' @param log.x logical; Apply log(x) conversion. Default \code{FALSE}
#' @param error.model character; Applied error model, \code{"proportional"} or \code{"additive"}. Default \code{"proportional"}
#' @param model.name character; Optional model name for plot output
#' @param mappings named character vector;  Optional mappings to be included if column names in input \code{data.frame} differ from required column names.
#' See Required Columns section:
#' @section Required columns \code{obs.data}:
#' \itemize{
#'  \item \code{"ID"} - Subject ID
#'  \item \code{"x"} - Typically time
#'  \item \code{"PRED"} - Population prediction, required if \code{iiv.correction = TRUE}
#'  \item \code{"IPRED"} - Individual prediction, required if \code{iiv.correction = TRUE}
#'  \item \code{"OBS"} - DV
#'  \item \code{"dosenr"} - Dose number; unique dose number for ID/time point
#' }
#'
#' @section Required columns \code{typ.data}:
#' \itemize{
#'  \item \code{"ID"} - Subject ID
#'  \item \code{"x"} - Typically time
#'  \item \code{"PRED"} - Population prediction
#'  \item \code{"dosenr"} - Dose number; unique dose number for ID/time point
#' }
#'
#' @section Required columns \code{sim.data}:
#' \itemize{
#'  \item \code{"ID"} - Subject ID
#'  \item \code{"x"} - Typically time
#'  \item \code{"PRED"} - Population prediction, required if \code{iiv.correction = TRUE}
#'  \item \code{"IPRED"} - Individual prediction, required if \code{iiv.correction = TRUE}
#'  \item \code{"REP"} - Replicate number
#' }
#'
#' @details If "dosenr" column is missing it will be automatically calculated using the priority of available columns:
#'    \itemize{
#'      \item "EVID": If available in data, "dosenr" will be calculated using \code{cumsum(EVID==1)}
#'      \item "ADDL"/"II": If "ADDL" and "II" are available in data, "dosenr" will be calculated given additional dose number and interval
#'      \item "AMT": If only "AMT" column exists in data, "dosenr" will be calculated using \code{cumsum(AMT!=0)}
#'    }
#'
#' @return \code{vachette_data}
#' @export
#'
#' @examples
#' obs <- read.csv(system.file(package = "vachette", "examples", "iv-obs.csv"))
#' typ <- read.csv(system.file(package = "vachette", "examples", "iv-typ-minmax.csv"))
#'
#' vd <- vachette_data(
#'   obs.data = obs,
#'   typ.data = typ,
#'   covariates = c(WT = 70),
#'   mappings = c(OBS = "DV", x = "time"),
#'   model.name  = "IV"
#'   )
#'
#'
vachette_data <-
  function(obs.data,
           typ.data,
           sim.data = NULL,
           covariates,
           ref.dosenr,
           log.x = FALSE,
           iiv.correction = FALSE,
           error.model = c("proportional", "additive"),
           model.name = NULL,
           mappings = NULL) {

    data_type <- PRED <- tab.ucov <- n.ucov <- NULL
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
        #sim.data <- left_join(sim.data, obs.data %>% select(ID, x, dosenr), by = c("ID", "x"))
        sim.data <- .calculate_dose_number(sim.data, data_type = "sim.data")
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

    if(log.x) {
      message("Applying log(x) conversion")
      if(sum(typ.orig$x)<=0) stop("Error, no log(x) conversion possible for x <= 0 (typical)")
      if(sum(obs.orig$x)<=0) stop("Error, no log(x) conversion possible for x <= 0 (observations)")
      if(sum(sim.orig$x)<=0) stop("Error, no log(x) conversion possible for x <= 0 (simulation)")
      typ.orig$x <- log(typ.orig$x)
      obs.orig$x <- log(obs.orig$x)
      sim.orig$x <- log(sim.orig$x)

      # if(xstart<=0 | xstop<=0) stop("Error, no log(x) conversion possible for observed/simulated x <= 0")
      # xstart <- log(xstart)
      # xstop  <- log(xstop)
    }

    # Check max x in typical curves is greater than max x in obs, if not, warn
    ucov <- unique(obs.orig$COV)

    for (cov in ucov){
      xmax_obs <- max(obs.orig[obs.orig$COV==cov,]$x)
      xmax_typ <- max(typ.orig[typ.orig$COV==cov,]$x)
      if (xmax_typ < xmax_obs) {
        warning("Maximum x value for ID = ", typ.orig[typ.orig$COV==cov,]$ID[1], " in `typ.data` is less than the maximum x value in `obs.data`, it is recommended to simulate longer:\n",
                "Unique covariate(s):\t", paste0(names(covariates), "=", strsplit(cov, split = " ")[[1]], collapse = "\t"), call. = FALSE)
      }
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
      log.x = log.x,
      # xstart = xstart,
      # xstop = xstop,
      n.ucov = n.ucov
    ) %>%
      structure(class = "vachette_data")
  }

update <- function(vachette_data, ...) {
  args <- list(...)
  for (i in names(args)) {
    vachette_data[[i]] <- args[[i, exact=TRUE]]
  }
  vachette_data
}

.define_and_enumerate_regions <-
  function(typ.orig,
           obs.orig,
           sim.orig,
           covariates,
           ref.dosenr,
           ref.region) {
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
    vars <- names(covariates)

    # Get all covariate combinations and ensure they exist in obs.orig
    comb.ucov <- typ.orig %>% tidyr::expand(!!!syms(vars))
    comb.ucov <- dplyr::inner_join(comb.ucov, obs.orig, by = vars) %>%
      dplyr::select(dplyr::all_of(vars)) %>%
      dplyr::distinct()

    for (i in c(1:nrow(comb.ucov))) {
      filter_query <-
        purrr::map2(vars, comb.ucov[i, vars], function(var_name, value) {
          rlang::expr(!!rlang::sym(var_name) == !!value)
        })
      # Get number of doses for covariate combination:
      z <- typ.orig %>% filter(!!!filter_query)

      ndose <- length(unique(z$dosenr))

      for (idose in c(1:ndose))
      {
        if (idose == ndose)
          region.type <- 'open'
        if (idose <  ndose)
          region.type <- 'closed'

        # New unique combination
        n.ucov <- n.ucov + 1

        add <- data.frame(ucov = n.ucov)
        vals <- comb.ucov %>% slice(i)
        add[vars] <- vals
        add <- mutate(add, region = idose,
                      region.type = region.type)

        tab.ucov <- rbind(tab.ucov, add)
        filter_query_typ <- paste0(
          paste0(
            "typ.orig$",
            vars,
            "==",
            "comb.ucov$",
            vars,
            "[i]",
            collapse = " & "
          ),
          " & typ.orig$dosenr == idose"
        )

        filter_query_obs <- paste0(
          paste0(
            "obs.orig$",
            vars,
            "==",
            "comb.ucov$",
            vars,
            "[i]",
            collapse = " & "
          ),
          " & obs.orig$dosenr == idose"
        )
        filter_query_sim <- paste0(
          paste0(
            "sim.orig$",
            vars,
            "==",
            "comb.ucov$",
            vars,
            "[i]",
            collapse = " & "
          ),
          " & sim.orig$dosenr == idose"
        )
        # Selection of new unique combination
        sel.typ <- eval(parse(text = filter_query_typ))
        sel.obs <- eval(parse(text = filter_query_obs))
        sel.sim <- eval(parse(text = filter_query_sim))

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

    tab.ucov$ref <- "No"

    for (i.ucov in c(1:nrow(tab.ucov)))
    {
      condition <- paste0(
        paste0(
          "as.character(tab.ucov$",
          names(covariates),
          "[i.ucov])",
          "==",
          "'",
          covariates,
          "'",
          collapse = " && "
        ),
        " & tab.ucov$region[i.ucov] == ref.region"
      )
      if (eval(parse(text = condition)))
        tab.ucov$ref[i.ucov] <- "Yes"
    }

    # Add reference flag to typ.orig
    typ.orig <- typ.orig %>%
      dplyr::left_join(tab.ucov[, c('ucov', 'ref')], by = 'ucov')

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
#' @param tol.end numeric; Relative tolerance to determine last x open end reference
#' @param tol.noise numeric; Relative tolerance for landmark determination typical curves
#' @param step.x.factor numeric; x-axis extension factor to search for last x, i.e., to determine where close enough to asymptote
#' @param ngrid.fit numeric; number of grid points in last query segment for matching last reference segment
#' @param window integer; size (gridpoints) of Savitzky Golay smoothing window for landmark position determination
#' @param log_file character; File path to direct console output e.g., \code{"log.txt"}
#' @param ... Additional arguments
#' @name apply_transformations
#' @return \code{vachette_data} object containing a list of vachette-transformed \code{data.frame}s
#' @examples
#' obs <- read.csv(system.file(package = "vachette", "examples", "iv-obs.csv"))
#' typ <- read.csv(system.file(package = "vachette", "examples", "iv-typ-minmax.csv"))
#'
#' vd <- vachette_data(
#'   obs.data = obs,
#'   typ.data = typ,
#'   covariates = c(WT = 70),
#'   mappings = c(OBS = "DV", x = "time"),
#'   model.name  = "IV"
#'  )
#'
#'  vd <- apply_transformations(vd)
#'
#' @export
apply_transformations <- function(vachette_data, ...) UseMethod("apply_transformations")

#' @rdname apply_transformations
#' @export
apply_transformations.vachette_data <-
  function(vachette_data,
           tol.end = 0.001,
           tol.noise = 1e-8,
           step.x.factor = 1.5,
           ngrid.fit = 100,
           window = 17,
           log_file = NULL,
           ...) {

    stopifnot(inherits(vachette_data, "vachette_data"))

    args <- list(...)

    if (!is.null(args$asymptote_right)) {
      asymptote_right <- args$asymptote_right
    } else {
      asymptote_right <- TRUE
    }


    if (!is.null(args$asymptote_left)) {
      asymptote_left <- args$asymptote_left
    } else {
      asymptote_left <- FALSE
    }

    if (!is.null(args$zero_asymptote_right)) {
      zero_asymptote_right <- args$zero_asymptote_right
    } else {
      zero_asymptote_right <- TRUE
    }

    if (!is.null(args$zero_asymptote_left)) {
      zero_asymptote_left <- args$zero_asymptote_left
    } else {
      zero_asymptote_left <- FALSE
    }

    if (!is.null(log_file)) {
      if (file.exists(log_file))  {
        res <- askYesNo(msg = paste0(log_file, " exists and will be overwritten"))
        if (res) {
          unlink(log_file)
        } else {
          stop("Cannot delete log file. Set argument `log_file`=NULL to ignore.")
        }
      }
      conn <- file(log_file, open = "at")
      assign("vachette_log_file", value = conn, envir = vachette_env)
      on.exit(close(conn))
      on.exit(assign("vachette_log_file", value = NULL, envir = vachette_env), add = TRUE)
    }

    vachette_transformed_data <-
      tryCatch({
        .calculate_transformations(
          vachette_data,
          tol.end,
          tol.noise,
          step.x.factor,
          ngrid.fit,
          window,
          asymptote_right = asymptote_right,
          asymptote_left  = asymptote_left,
          zero_asymptote_right = zero_asymptote_right,
          zero_asymptote_left  = zero_asymptote_left
        )
      }, error = function(e) {
        warning(e$message)
        # vachette_summary <- get("summary_out", envir = vachette_env)
        vachette_env$summary_out$values[["asymptote_right"]] <- asymptote_right
        vachette_env$summary_out$values[["asymptote_left"]] <- asymptote_left
        vachette_env$summary_out$values[["zero_asymptote_right"]] <- zero_asymptote_right
        vachette_env$summary_out$values[["zero_asymptote_left"]] <- zero_asymptote_left
        return(update(vachette_data, summary = vachette_env$summary_out))
      })

    vachette_env$summary_out$values[["asymptote_right"]] <- asymptote_right
    vachette_env$summary_out$values[["asymptote_left"]] <- asymptote_left
    vachette_env$summary_out$values[["zero_asymptote_right"]] <- zero_asymptote_right
    vachette_env$summary_out$values[["zero_asymptote_left"]] <- zero_asymptote_left

    return(update(vachette_data,
             obs.all = vachette_transformed_data$obs.all,
             obs.excluded = vachette_transformed_data$obs.excluded,
             sim.all = vachette_transformed_data$sim.all,
             curves.all = vachette_transformed_data$curves.all,
             curves.scaled.all = vachette_transformed_data$curves.scaled.all,
             ref.extensions.all = vachette_transformed_data$ref.extensions.all,
             lm.all = vachette_transformed_data$lm.all,
             curvature.all = vachette_transformed_data$curvature.all,
             nseg = vachette_transformed_data$nseg,
             summary = vachette_env$summary_out))

}

