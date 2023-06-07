
#' Initialize vachette object with required data
#'
#' @param obs.data data.frame; Observed data
#' @param typ.data data.frame; Typical (population) curves
#' @param sim.data data.frame; Simulated (VPC) data
#' @param covariates named character vector; Covariate names with reference values in vachette transformation
#' @param ref.dosenr integer; Number of doses given up through start of reference region
#' @param iiv.correction logical; Apply inter-individual variability correction. Default \code{FALSE}
#' @param error.model character; Applied error model, \code{"proportional"} or \code{"additive"}. Default \code{"proportional"}
#' @param model.name character; Optional model name for plot output
#' @param mappings named character vector;  Optional mappings to be included if column names in input \code{data.frame} differ from required column names
#' @details Required column names in input \code{data.frame} are:
#' \itemize{
#'  \item{"ID"}{ - Subject ID}
#'  \item{"x"}{ - Typically time}
#'  \item{"PRED"}{ - Population prediction, required if \code{iiv.correction = TRUE}}
#'  \item{"IPRED"}{ - Individual prediction, required if \code{iiv.correction = TRUE}}
#'  \item{"OBS"}{ - DV}
#' }
#'
#'
#' @return \code{vachette_data}
#' @export
#'
#' @examples
#' obs.data <- read.csv(system.file(package = "vachette", "examples", "iv-obs.csv"))
#' typ.data <- read.csv(system.file(package = "vachette", "examples", "iv-typ.csv"))
#' sim.data <- read.csv(system.file(package = "vachette", "examples", "iv-vpc.csv"))
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

  # If ref.dosenr is missing, provide warning then set to 1
  if (missing(ref.dosenr)) {
    warning("ref.dosenr argument not given, setting ref.dosenr to 1", call. = FALSE)
    ref.dosenr <- 1
  }

  obs.data <- .calculate_dose_number(obs.data, ref.dosenr = ref.dosenr, data_type = "obs.data")
  typ.data <- .calculate_dose_number(typ.data, ref.dosenr = ref.dosenr, data_type = "typ.data")

  # Region is the Vachette terminology for the time between two dose administrations
  ref.region <- ref.dosenr
  # "Dummy" observations replicate number
  obs.data$isim <- 1
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
    if (iiv.correction) {
      warning("sim.data argument used with iiv.correction = TRUE, setting iiv.correction = FALSE")
      iiv.correction <- FALSE
    }
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
    cols_to_select <- c("isim", "ID", "x", "PRED", "IPRED", "OBS", "dosenr")
  } else {
    cols_to_select <- c("isim", "ID", "x", "OBS", "dosenr")
  }
  obs.data  <-
    obs.data %>% dplyr::select(dplyr::all_of(cols_to_select), dplyr::all_of(names(covariates)))
  sim.data  <-
    sim.data %>% dplyr::select(dplyr::all_of(cols_to_select), dplyr::all_of(names(covariates)))

  # Extract last observed x (look for max x/time in obs and sim data)
  xstop <- max(obs.data$x,sim.data$x)
  # Define unique covariate combination as new 'COV' column.
  obs.orig <- obs.data %>%
    mutate(COV = paste(!!!syms(names(covariates))))
  sim.orig <- sim.data %>%
    mutate(COV = paste(!!!syms(names(covariates))))
  typ.data <- typ.data %>%
    mutate(COV = paste(!!!syms(names(covariates))))

  .define_and_enumerate_regions(typ.data, obs.orig, sim.orig, covariates, ref.dosenr, ref.region) %>%
    list2env(envir = vachette_data_env)

  # Apply IIV Correction
  typ.data <- typ.data %>%
    mutate(y=PRED)
  # to apply
  if(iiv.correction) obs.orig$y <- obs.orig$OBS  - obs.orig$IPRED + obs.orig$PRED
  if(iiv.correction) sim.orig$y <- sim.orig$OBS  - sim.orig$IPRED + sim.orig$PRED

  # to omit
  if(!iiv.correction) obs.orig$y <- obs.orig$OBS
  if(!iiv.correction) sim.orig$y <- sim.orig$OBS

  list(
    model.name = model.name,
    obs.data = obs.data,
    typ.data = typ.data,
    sim.data = sim.data,
    covariates = covariates,
    VVPC = VVPC,
    ADD_TR = ADD_TR,
    PROP_TR = PROP_TR,
    typ.data = typ.data,
    obs.orig = obs.orig,
    sim.orig = sim.orig,
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

.define_and_enumerate_regions <- function(typ.data, obs.orig, sim.orig, covariates, ref.dosenr, ref.region) {
  # ----- Define regions and provide a number to all combined covariate effects ----
  # Region number
  typ.data$region  <- NA
  obs.orig$region    <- NA
  sim.orig$region    <- NA

  # Region type (open/closed)
  typ.data$region.type  <- NA
  obs.orig$region.type    <- NA
  sim.orig$region.type    <- NA

  # Covariate combinations number
  typ.data$ucov  <- NA
  obs.orig$ucov    <- NA
  sim.orig$ucov    <- NA

  n.ucov    <- 0
  tab.ucov  <- NULL # Table with ucov properties

  # Get all covariate combinations
  comb.ucov <- typ.data %>% tidyr::expand(!!!syms(names(covariates)))

  # Collect unique covariate/dose combinations and define region.type
  # Add to info to data frames
  if(length(covariates)==1)
  {
    for(i in c(1:dim(comb.ucov)[1])) #e.g., unique values of covariates nrow(comb.ucov)
    {
      # Get number of doses for each covariate combination:
      z <- typ.data %>% filter(!!sym(names(covariates)) == comb.ucov[names(covariates)][[1]][i])
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
        sel.typ  <- typ.data[[names(covariates)]] == comb.ucov[names(covariates)][[1]][i] &
          typ.data$dosenr == idose
        sel.obs  <- obs.orig[[names(covariates)]] == comb.ucov[names(covariates)][[1]][i] &
          obs.orig$dosenr == idose
        sel.sim  <- sim.orig[[names(covariates)]] == comb.ucov[names(covariates)][[1]][i] &
          sim.orig$dosenr == idose

        # Add new unique combination info
        typ.data$region[sel.typ]   <- idose
        obs.orig$region[sel.obs]     <- idose
        sim.orig$region[sel.sim]     <- idose

        typ.data$ucov[sel.typ]   <- n.ucov
        obs.orig$ucov[sel.obs]     <- n.ucov
        sim.orig$ucov[sel.sim]     <- n.ucov

        typ.data$region.type[sel.typ]   <- region.type
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
      z <- typ.data %>% filter(eval(parse(text = filter_query)))

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
          paste0("typ.data$", vars, "==", "comb.ucov$", vars, "[i]", collapse = " & "),
          " & typ.data$dosenr == idose"
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
        # sel.typ  <- typ.data$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        #   typ.data$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        #   typ.data$dosenr == idose
        # sel.obs  <- obs.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        #   obs.orig$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        #   obs.orig$dosenr == idose
        # sel.sim  <- sim.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        #   sim.orig$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        #   sim.orig$dosenr == idose

        # Add new unique combination info
        typ.data$region[sel.typ]   <- idose
        obs.orig$region[sel.obs]     <- idose
        sim.orig$region[sel.sim]     <- idose

        typ.data$ucov[sel.typ]   <- n.ucov
        obs.orig$ucov[sel.obs]     <- n.ucov
        sim.orig$ucov[sel.sim]     <- n.ucov

        typ.data$region.type[sel.typ]   <- region.type
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

  # Add reference flag to typ.data
  typ.data <- typ.data %>%
    dplyr::left_join(tab.ucov[,c('ucov','ref')],by='ucov')

  return(
    list(
      typ.data = typ.data,
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
    # Collect all Vachette query curves and original/transformed observations (incl reference)
    ref.extensions.all   <- NULL
    curves.all           <- NULL
    curves.scaled.all    <- NULL
    obs.all              <- NULL

    my.ref.lm.all   <- NULL
  my.query.lm.all <- NULL
  lm.all <- NULL

  tab.ucov <- vachette_data$tab.ucov
  stopifnot(!is.null(tab.ucov))
  typ.data <- vachette_data$typ.data
  stopifnot(!is.null(typ.data))
  obs.orig <- vachette_data$obs.orig
  stopifnot(!is.null(obs.orig))
  sim.orig <- vachette_data$sim.orig
  stopifnot(!is.null(sim.orig))

  covariates <- vachette_data$covariates
  ref.region <- vachette_data$ref.region
  ref.dosenr <- vachette_data$ref.dosenr
  #lines 399-760
  # # Loop for all combinations of covariates
  for(i.ucov in c(1:dim(tab.ucov)[1])) {
    # A. ----- Define reference and query typical curves and observations to transform  ----------

    # Ref: may change (by extensions), so define (again) every new combination
    filter_query <-
      paste0(
        "!is.na(region) & ",
        paste0(
          'as.character(',
          names(covariates),
          ')',
          "==",
          "'",
          covariates,
          "'",
          collapse = " & "
        ),
        " & region == ref.region"
      )

    ref <- typ.data %>% filter(eval(parse(text = filter_query)))

    ref.ucov <- unique(ref$ucov)

    if(length(ref.ucov) != 1) stop("Error: length ref.ucov != 1")     # A single ref only possible

    ref.region.type <- tab.ucov$region.type[tab.ucov$ucov==ref.ucov]
    # ref = Typical reference curve
    query.region       <- tab.ucov$region[tab.ucov$ucov==i.ucov]
    query.region.type  <- tab.ucov$region.type[tab.ucov$ucov==i.ucov]

    # Typical query curve
    query <- typ.data %>%
      filter(ucov == i.ucov) %>%
      mutate(ref = tab.ucov$ref[i.ucov])    # Flag for reference

    # 230418 - extrapolate region last-x region by one gridstep size if region.type = 'closed'
    # Currently simple extra x value with same y value ("horizontal" extrapolation - LOCF)
    if(ref.region.type == 'closed')
    {
      message(paste0("Reference region-",ref.region," gap extrapolation"))

      ref.grid.step.size   <- ref$x[dim(ref)[1]] - ref$x[dim(ref)[1]-1]
      ref.curve.next       <- ref[dim(ref)[1],]
      # Increase x:
      ref.curve.next$x     <- ref.curve.next$x + ref.grid.step.size

      # Last 6 datapoints
      last6 <- ref %>% slice((n()-5):n()) %>% select(x,y)

      # Exponential fit
      y  = last6$y
      x  = last6$x
      exp.model <-lm(log(y) ~ x)

      # Extrapolate y
      ref.curve.next$y     <- exp(predict(exp.model,list(x=ref.curve.next$x)))

      ref <- rbind(ref,ref.curve.next)

    }
    if(query.region.type == 'closed')
    {
      message(paste0("Query region-",query.region," gap extrapolation"))

      query.grid.step.size   <- query$x[dim(query)[1]] - query$x[dim(query)[1]-1]
      query.curve.next       <- query[dim(query)[1],]
      # Increase x:
      query.curve.next$x     <- query.curve.next$x + query.grid.step.size

      # Last 6 datapoints
      last6 <- query %>% slice((n()-5):n()) %>% select(x,y)

      # Exponential fit
      y  = last6$y
      x  = last6$x
      exp.model <-lm(log(y) ~ x)

      # Extrapolate y
      query.curve.next$y     <- exp(predict(exp.model,list(x=query.curve.next$x)))

      query <- rbind(query,query.curve.next)
    }

    # We have to run Vachette twice if observation AND simulated replicates have to be generated (for VPC)

    VVPC <- vachette_data$VVPC
    # Observations
    if (!VVPC)
    {
      obs.query   <- obs.orig %>%
        filter(ucov == i.ucov) %>%
        mutate(ref = tab.ucov$ref[i.ucov])  # Flag for reference data
    }
    # Simulated observations for VPC
    if (VVPC)
    {
      obs.query   <- sim.orig %>%
        filter(ucov == i.ucov)  %>%
        mutate(ref = tab.ucov$ref[i.ucov])  # Flag for reference data
    }

    # JL 230606
    # Create PRED for obs.query data points by linear interpolation
    # May be needed for reference extrapolation, see obs.query$PRED
    obs.query$PRED <- approx(x=query$x,y=query$y,xout=obs.query$x)$y

    #  B. -------- Get landmarks approximate positions ----

    # Tolerance to apply for finding maximums, minimums and inflection points
    # Small stepsize -> small tolerance
    # Large stepsize -> large tol
    #tolapply = distance between second and first x obs * tol.noise
    tolapply <- tol.noise*(typ.data$x[2]-typ.data$x[1])
    # get ref landmarks
    my.ref.lm.init    <- get.x.multi.landmarks(ref$x,ref$y,w=window,tol=tolapply) #contiguous inflec cause issue?
    my.ref.lm.init$y  <- approx(ref$x,ref$y, xout=my.ref.lm.init$x)$y #interpolation
    # get query landmarks
    my.query.lm.init    <- get.x.multi.landmarks(query$x,query$y,w=window,tol=tolapply)
    my.query.lm.init$y  <- approx(query$x,query$y, xout=my.query.lm.init$x)$y

    # stop at this point, perhaps, allow user to visually inspect this plot, to assess what increase in stepsize
    # do not error out, return data that is available, and provide to user for plotting, and provide suggestions for updating argument values

    #Validation check
    #If y contiguous y values in series are the same, provide error, increase tol.noise
    #option 1: if multiple inflec in sequence, try decreasing tol.noise

    if(!lm.refine)
    {
      my.ref.lm.refined      <- my.ref.lm.init
      my.query.lm.refined    <- my.query.lm.init
    }
    if(lm.refine)
    {
      my.ref.lm.refined    <- refine.x.multi.landmarks(x=ref$x,y=ref$y,lm=my.ref.lm.init,tol=0.01*tolapply, w1=window.d1.refine,w2=window.d2.refine)
      my.ref.lm.refined$y  <- approx(ref$x, ref$y, xout=my.ref.lm.refined$x)$y

      my.query.lm.refined    <- refine.x.multi.landmarks(query$x,query$y,lm=my.query.lm.init,tol=0.01*tolapply, w1=window.d1.refine,w2=window.d2.refine)
      my.query.lm.refined$y  <- approx(query$x,query$y, xout=my.query.lm.refined$x)$y
    }

    # Check if query and ref have same landmarks in same order
    if(dim(my.ref.lm.refined)[1] != dim(my.query.lm.refined)[1] |
       !sum(my.ref.lm.refined$type == my.query.lm.refined$type)>0)
      stop("No matching landmarks between reference and query")

    # C. ---------- Get last.x for reference and query typical curves ------------------

    scaling <- 'linear' # other options?

    # ------ Landmark last.x determination depends on segment type ---------

    # 1. ref and query both open ends ("default")
    # 2. ref is open end, query is closed --> get.query.open.end() for last x **ref** only
    # 3. ref is closed, query is open end --> get.query.open.end() for last x query only
    # 4. ref and query both closed        --> get.query.open.end() for last x query only

    query.region.type <- unique(query$region.type)
    if(length(query.region.type) != 1) stop("Error: length query.region.type != 1")

    # ---- Reference and query last x -----

    # 1. Find ref last x, Fit query last x
    if(ref.region.type == 'open' & query.region.type == 'open')
    {
      my.ref.lm   <- get.ref.x.open.end(ref$x,ref$y,my.ref.lm.refined,step.x.factor=step.x.factor,tol=tol.end)
      my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.fit,scaling=scaling)
    } #error is in above line
    # 2. Fit ref last x, Fix query last x
    if(ref.region.type == 'open' & query.region.type == 'closed')
    {
      my.ref.lm   <- get.query.x.open.end(query,ref,my.query.lm.refined,my.ref.lm.refined,ngrid=ngrid.fit,scaling=scaling)
      my.query.lm <- my.query.lm.refined
    }
    # 230418 - Update
    # 3+4. Fix ref last x, Fit query last x
    if(ref.region.type == 'closed')
    {
      my.ref.lm   <- my.ref.lm.refined
      my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.fit,scaling=scaling)
    }

    # @James recent developments
    # Overrides above (above needs reprogramming)
    # Always FIX ref and FIT query to find best fitting x ("open"). Then extrapolate if required
    if(ref.region.type == 'closed')
    {
      my.ref.lm   <- my.ref.lm.refined
      my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.fit,scaling=scaling)
    }

    # Recalc landmark y's
    my.query.lm$y  <- approx(query$x,query$y, xout=my.query.lm$x)$y
    my.ref.lm$y    <- approx(ref$x,ref$y, xout=my.ref.lm$x)$y

    # Collect all landmarks
    my.ref.lm.all   <- rbind(my.ref.lm.all,   my.ref.lm   %>% mutate(i.ucov = i.ucov))
    my.query.lm.all <- rbind(my.query.lm.all, my.query.lm %>% mutate(i.ucov = i.ucov))

    lm.all <- rbind(lm.all, my.query.lm %>% mutate(ucov = i.ucov))

    # D. ----- Define and map segments -----------

    nseg     <- dim(my.ref.lm)[1]-1

    # Initial assignment segments to typical curves, collect all query curves
    ref$seg   <- NA
    query$seg <- NA
    for(iseg in c(1:nseg))
    {
      if(iseg==1) # includes first datapoint
      {
        ref <- ref %>%
          mutate(seg = ifelse(x>=my.ref.lm$x[iseg] & x<=my.ref.lm$x[iseg+1],iseg,seg))
        query <- query %>%
          mutate(seg = ifelse(x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1],iseg,seg))
      }
      if(iseg>1)
      {
        ref <- ref %>%
          mutate(seg = ifelse(x>my.ref.lm$x[iseg] & x<=my.ref.lm$x[iseg+1],iseg,seg))
        query <- query %>%
          mutate(seg = ifelse(x>my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1],iseg,seg))
      }
    }

    # E. ---- x-scaled query segments -----------

    # 230418 New query.scaled curve data frame by contracting/expanding to x-range ref
    query.scaled <- NULL
    for(iseg in c(1:nseg))
    {
      # Ref segment length
      ref.seg.length <- my.ref.lm$x[iseg+1] - my.ref.lm$x[iseg]
      # Query segment length
      query.seg.length <- my.query.lm$x[iseg+1] - my.query.lm$x[iseg]
      # Scaling factor
      my.query.add <- query %>%
        filter((nseg==1 & x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1]) |
                 (nseg>1 & x>my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1])) %>%
        mutate(seg=iseg) %>%
        mutate(y.scaling     = 1, y.scaled =1) %>%                # dummies, To remove
        mutate(x.shift.ref   = 0 - my.ref.lm$x[iseg]) %>%
        mutate(x.shift.query = 0 - my.query.lm$x[iseg]) %>%
        mutate(x.trans       = x + x.shift.query) %>%             # May be redefine to my.ref.lm$x[iseg]
        mutate(x.scaling     = ifelse(query.seg.length>0,
                                      ref.seg.length/query.seg.length,
                                      NA)) %>%                    # NA causes troubles?
        mutate(x.scaled      = x.scaling*x.trans - x.shift.ref)   # Back to corresponding ref position

      query.scaled <- rbind(query.scaled,my.query.add)
    }

    # F. ---- Check for need of extended **reference** curve by extrapolation -----------

    # Default no extension for extrapolation
    # tab.extension <- NULL
    EXTENSION <- FALSE

    # Check if there are observation at all
    OBSERV <- ifelse(dim(obs.query)[1]>0,T,F)

    # If all observations are not before query last x, then extrapolate and add to reference region curve
    if (OBSERV)
    {
      # Only needed when there are observation outside the segment last landmark
      if(max(my.query.lm$x) < max(obs.query$x))
      {

        EXTENSION <- TRUE

        message('*** Extension reference curve by simple exponential extrapolation ***')

        # max x observation on query
        max.obs.x <- max(obs.query$x)
        # y observation on query at x max
        max.obs.y <- obs.query$y[obs.query$x==max.obs.x]
        # y typical on query at x max
        max.obs.y.typ <- obs.query$PRED[obs.query$x==max.obs.x]

        # scaled max x observation: (x+x.shift)*scaling
        lastpoint        <- dim(query.scaled)[1]
        max.obs.x.scaled <- (max.obs.x + query.scaled$x.shift.query[lastpoint])*query.scaled$x.scaling[lastpoint]

        # QUERY
        cur.query.seg <- query$seg[dim(query[!is.na(query$seg),])[1]]
        # Adjust my.query.lm (always last point of last segment)
        my.query.lm.extension <- my.query.lm
        my.query.lm.extension$x[nseg+1] <- max.obs.x
        my.query.lm.extension$y[nseg+1] <- max.obs.y.typ

        # REFERENCE
        # Reference curve has to be extrapolated up to max.obs.x.scaled
        # Last segment reference with last two typical datapoints extrapolated to scaled max x for obs

        # lastpoint        <- dim(ref)[1]
        # add.x <- max.obs.x.scaled
        # add.y <- approxExtrap(c(ref$x[lastpoint-1],ref$x[lastpoint]),
        #                       c(ref$y[lastpoint-1],ref$y[lastpoint]),
        #                       xout = add.x)$y
        #
        # ref.add <- ref %>%
        #   slice(1) %>%
        #   mutate(x=add.x,
        #          y=add.y) %>%
        #   mutate(seg=cur.ref.seg)
        # ref <- rbind(ref,ref.add)

        # Same as for interregion gap extrapolation - multiple points
        cur.ref.seg   <- ref$seg[dim(ref[!is.na(ref$seg),])[1]]

        ref.grid.step.size   <- ref$x[dim(ref)[1]] - ref$x[dim(ref)[1]-1]

        # Increase x by ref grid steps up to last point
        # number of steps: ceiling((max.obs.x - max(ref$x))/ref.grid.step.size)
        # JL 230607. Correction, n.steps now defined as total number of steps, isteps as step counts
        i.steps.x          <- c(1:ceiling((max.obs.x - max(ref$x))/ref.grid.step.size))
        n.steps.x          <- length(i.steps.x)
        steps.x            <- max(ref$x) + i.steps.x*ref.grid.step.size

        ref.curve.next     <- ref %>% slice(rep(n(),max(i.steps.x)))
        ref.curve.next$x   <- steps.x
        # Last x exactly matching max query obs x
        ref.curve.next$x[n.steps.x] <- max.obs.x.scaled

        # Last 6 datapoints
        last6 <- ref %>% slice((n()-5):n()) %>% select(x,y)

        # Exponential fit
        y  = last6$y
        x  = last6$x
        exp.model <-lm(log(y) ~ x)

        # Extrapolate y
        ref.curve.next$y     <- exp(predict(exp.model,list(x=ref.curve.next$x)))

        # JL 230607 Assign extrapolation reference part to seg=-99 for plotting purposes
        ref.curve.next$seg <- -99

        ref <- rbind(ref,ref.curve.next)

        # JL 230607
        # Keep extensions with associated query curve ucov
        ref.curve.next$ucov <- i.ucov  # Included reference curve extrapolation
        ref.extensions.all  <- rbind(ref.extensions.all,ref.curve.next)  # Included reference curve extrapolation

        # Check:
        if(cur.ref.seg != cur.query.seg) stop("Error segment number assignment in reference extrapolation block")
      }
    }

    # Collect all typical curves with scaling factors and scaled x,y values
    query.scaled$y.scaled <- approx(ref$x,ref$y,xout=query.scaled$x.scaled)$y

    # JL 230607 Assign part of curve mapped to extrapolated part of reference curve, if required
    # Only extension part....
    if(EXTENSION) query.scaled$seg <- approx(ref$x,ref$seg,xout=query.scaled$x.scaled,
                                             method = 'constant')$y

    curves.all            <- rbind(curves.all,query)              # Included reference curve extrapolation
    curves.scaled.all     <- rbind(curves.scaled.all,query.scaled)

    if(dim(obs.query)[1]>0)
    {
      # Steps
      # a. (Re-)assign segments
      # b. Determine obs new x-position
      # c. Determine obs y-scaling ADDITIVE OR PROPORTIONAL

      # a. (Re-)assign segments, x.shift.ref, x.shift.query and x.scaling to observations and query curve (ref already done)
      obs.query$seg <- NA
      for(iseg in c(1:nseg))
      {
        # extended.query.x <- ifelse(length(tab.extension$x[tab.extension$seg==iseg])>0,
        #                            tab.extension$x[tab.extension$seg==iseg],NA)
        first.query.x    <- my.query.lm$x[iseg]
        # last.query.x     <- ifelse(!is.na(extended.query.x),extended.query.x,my.query.lm$x[iseg+1])
        if(!EXTENSION) last.query.x     <- my.query.lm$x[iseg+1]
        if(EXTENSION)  last.query.x     <- my.query.lm.extension$x[iseg+1]

        obs.query <- obs.query %>% mutate(seg = ifelse((nseg==1 & x>=first.query.x & x<=last.query.x) |
                                                         (nseg>1 & x>first.query.x & x<=last.query.x),iseg, seg))
        query     <- query     %>% mutate(seg = ifelse((nseg==1 & x>=first.query.x & x<=last.query.x) |
                                                         (nseg>1 & x>first.query.x & x<=last.query.x),iseg, seg))
      }

      # obs.query still OK

      # ----- Steps 3: Vachette x transformation -----

      # b. Determine obs new x-position
      obs.query.orig <- obs.query
      obs.query      <- NULL
      for(iseg in c(1:nseg))
      {
        # Original query.lm needed for scaling factor
        # Ref segment length
        ref.seg.length <- my.ref.lm$x[iseg+1] - my.ref.lm$x[iseg]
        # Query segment length
        query.seg.length <- my.query.lm$x[iseg+1] - my.query.lm$x[iseg]
        # Scaling factor
        if(!EXTENSION) obs.query.tmp <- obs.query.orig %>%
          filter((nseg==1 & x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1]) |
                   (nseg>1 & x>my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1]))
        if(EXTENSION) obs.query.tmp <- obs.query.orig %>%
          filter((nseg==1 & x>=my.query.lm.extension$x[iseg] & x<=my.query.lm.extension$x[iseg+1]) |
                   (nseg>1 & x>my.query.lm.extension$x[iseg] & x<=my.query.lm.extension$x[iseg+1]))
        obs.query.add <- obs.query.tmp %>%
          mutate(seg=iseg) %>%
          mutate(y.scaling     = 1, y.scaled =1) %>%               # dummies, To remove
          mutate(x.shift.ref   = 0 - my.ref.lm$x[iseg]) %>%
          mutate(x.shift.query = 0 - my.query.lm$x[iseg]) %>%
          mutate(x.trans       = x + x.shift.query) %>%            # May be redefine to my.ref.lm$x[iseg]
          mutate(x.scaling     = ref.seg.length/query.seg.length) %>%
          mutate(x.scaled      = x.scaling*x.trans - x.shift.ref)  # Back to corresponding ref position

        obs.query <- rbind(obs.query,obs.query.add)
      }

      # ----- Steps 4: Vachette y transformation -----

      # c. Determine obs y-scaling ADDITIVE OR PROPORTIONAL

      # In case of proportional error model
      PROP_TR <- vachette_data$PROP_TR
      if(PROP_TR)
      {
        # PROPORTIONAL - log scale addition. Take full ref and curve curves ....
        ref            <- ref   %>% mutate(ylog = ifelse(y>0,log(y),NA))
        query          <- query %>% mutate(ylog = ifelse(y>0,log(y),NA))
        obs.query      <- obs.query %>% mutate(ylog = ifelse(y>0,log(y),NA))

        # Position on query curve
        obs.query$ylog.query       <- ifelse(!is.na(obs.query$ylog),
                                             approx(x=query$x, y=query$ylog, xout = obs.query$x)$y,
                                             NA)
        # log difference obs.query to query curve
        obs.query$ylog.diff        <- obs.query$ylog - obs.query$ylog.query

        # Position on ref curve (position at x.scaled!!)
        obs.query$ylog.ref         <- ifelse(!is.na(obs.query$ylog),
                                             approx(x=ref$x, y=ref$ylog, xout = obs.query$x.scaled)$y,
                                             NA)
        # New position
        obs.query$ylog.scaled      <- obs.query$ylog.diff + obs.query$ylog.ref
        obs.query$y.scaled         <- exp(obs.query$ylog.scaled)
      }

      # In case of additive error model
      ADD_TR <- vachette_data$ADD_TR
      if(ADD_TR)
      {
        # ADDITIVE - normal scale addition

        # Position on query curve
        obs.query$y.query       <- ifelse(!is.na(obs.query$y),
                                          approx(x=query$x, y=query$y, xout = obs.query$x)$y,
                                          NA)
        # difference obs.query to query curve
        obs.query$y.diff        <- obs.query$y - obs.query$y.query

        # Position on ref curve (position at x.scaled!!)
        obs.query$y.ref         <- ifelse(!is.na(obs.query$y),
                                          approx(x=ref$x, y=ref$y, xout = obs.query$x.scaled)$y,
                                          NA)
        # New position
        obs.query$y.scaled      <- obs.query$y.diff + obs.query$y.ref
      }

      obs.all <- rbind(obs.all,obs.query)

    } # If there are observations to transform

  }

  n.ucov <- vachette_data$n.ucov
  for(i.ucov in c(1:n.ucov))
  {
    typ.curve <- typ.data %>% filter(ucov==i.ucov)
    ref.curve <- typ.data %>% filter(ucov==ref.ucov)
    obs       <- obs.all    %>% filter(ucov==i.ucov)

    #Init columns
    obs.all$dist.add.orig <- NA
    obs.all$dist.add.transformed <- NA
    obs.all$dist.prop.orig <- NA
    obs.all$dist.prop.transformed <- NA
    # Additive distances to original curves
    obs.all$dist.add.orig[obs.all$ucov==i.ucov] <- approx(typ.curve$x,typ.curve$y,xout=obs$x)$y - obs$y
    # Vachette transformed - distances to ref curve
    obs.all$dist.add.transformed[obs.all$ucov==i.ucov] <- approx(ref.curve$x,ref.curve$y,xout=obs$x.scaled)$y - obs$y.scaled

    # Proportional distances to original curves (assume all y > 0)
    # obs.all$dist.prop.orig[obs.all$ucov==i.ucov]        <- log(approx(typ.curve$x,typ.curve$y,xout=obs$x)$y) - log(obs$y)
    obs.all$dist.prop.orig[obs.all$ucov==i.ucov]        <- approx(typ.curve$x,log(typ.curve$y),xout=obs$x)$y - log(obs$y)
    # Vachette transformed - distances to ref curve
    # obs.all$dist.prop.transformed[obs.all$ucov==i.ucov] <- log(approx(ref.curve$x,ref.curve$y,xout=obs$x.scaled)$y) - log(obs$y.scaled)
    obs.all$dist.prop.transformed[obs.all$ucov==i.ucov] <- approx(ref.curve$x,log(ref.curve$y),xout=obs$x.scaled)$y - log(obs$y.scaled)
  }

  obs.all <- obs.all %>%
    arrange(isim, ID, x)

    update(vachette_data, obs.all = obs.all, curves.all = curves.all,
           curves.scaled.all = curves.scaled.all,
           ref.extensions.all = ref.extensions.all,
           lm.all = lm.all, nseg = nseg)


  }

