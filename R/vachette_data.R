
#' Initialize vachette object with required data
#'
#' @param indivsam.obs Observed data of class \code{data.frame}
#' @param output.typ Typ data of class \code{data.frame}
#' @param indivsam.vpc VPC data of class  \code{data.frame}
#' @param vachette.covs Named character vector of covariate names and reference values
#' @param ref.dose Reference dose value in data
#' @param IIV_CORR Logical; if IIV corrections performed
#' @param error Character; options are \code{"proportional"} or \code{"additive"}. Default is \code{"proportional"}
#' @param model.name Character; optional model name for plot output
#' @param mappings Named character vector specifying model variable name to column name in input data
#' @details Required column names in input \code{data.frame} are:
#' \itemize{
#'  \item{"ID"}{ - Subject ID}
#'  \item{"x"}{ - Typically time}
#'  \item{"PRED"}{ - Individual prediction}
#'  \item{"IPRED"}{ - Population prediction}
#'  \item{"OBS"}{ - DV}
#' }
#' If the above column names are different in your input data, use the \code{mappings} argument.
#'
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
vachette_data <-
  function(indivsam.obs,
           output.typ,
           indivsam.vpc = NULL,
           vachette.covs,
           ref.dose,
           IIV_CORR = FALSE,
           error = c("proportional", "additive"),
           model.name = NULL,
           mappings = NULL) {

  vachette_data_env <- environment()
  # Column Validation Check
  .validate_columns(mappings, indivsam.obs, output.typ, IIV_CORR) %>%
    list2env(envir = vachette_data_env)

  # Process/assign covariates
  vachette.covs <- .process_covariates(vachette.covs, indivsam.obs)

  # If ref.dosenr is missing, provide warning then set to 1
  if (missing(ref.dose)) {
    warning("ref.dosenr argument not given, setting ref.dosenr to 1")
    ref.dose <- 1
    #sigmoid model, does not need ref.dosenr
    #if no AMT,
  }

  indivsam.obs <- .calculate_dose_number(indivsam.obs, ref.dosenr = ref.dose, data_type = "obs.data")
  output.typ <- .calculate_dose_number(output.typ, ref.dosenr = ref.dose, data_type = "typ.data")

  # Region is the Vachette terminology for the time between two dose administrations
  ref.region <- ref.dose
  # "Dummy" observations replicate number
  indivsam.obs$isim <- 1
  error <- match.arg(error)
  # assign error flags
  if (error == "proportional") {
    ADD_TR <- FALSE
    PROP_TR <- TRUE
  } else {
    ADD_TR <- TRUE
    PROP_TR <- FALSE
  }

  stopifnot(names(vachette.covs) %in% names(output.typ))

  if (!is.null(indivsam.vpc)) {
    VVPC <- TRUE
    stopifnot(names(vachette.covs) %in% names(indivsam.vpc))
    if (IIV_CORR) {
      warning("indivsam.vpc argument used with IIV_CORR = TRUE, setting IIV_CORR = FALSE")
      IIV_CORR <- FALSE
    }
  } else {
    VVPC <- FALSE
    # "Dummy" vpc simulated observations dataset
    indivsam.vpc <- indivsam.obs
  }


  # Retain required columns
  indivsam.obs  <-
    indivsam.obs %>% dplyr::select(isim, ID, x, PRED, IPRED, OBS, dplyr::all_of(names(vachette.covs)), dosenr)
  indivsam.vpc  <-
    indivsam.vpc %>% dplyr::select(isim, ID, x, PRED, IPRED, OBS, dplyr::all_of(names(vachette.covs)), dosenr)

  # Extract last observed x (look for max x/time in obs and sim data)
  xstop <- max(indivsam.obs$x,indivsam.vpc$x)
  # Define unique covariate combination as new 'COV' column.
  obs.orig <- indivsam.obs %>%
    mutate(COV = paste(!!!syms(names(vachette.covs))))
  sim.orig <- indivsam.vpc %>%
    mutate(COV = paste(!!!syms(names(vachette.covs))))
  output.typ <- output.typ %>%
    mutate(COV = paste(!!!syms(names(vachette.covs))))

  .define_and_enumerate_regions(output.typ, obs.orig, sim.orig, vachette.covs, ref.dose, ref.region) %>%
    list2env(envir = vachette_data_env)

  # Apply IIV Correction
  output.typ <- output.typ %>%
    mutate(y=PRED)
  # to apply
  if(IIV_CORR) obs.orig$y <- obs.orig$OBS  - obs.orig$IPRED + obs.orig$PRED
  if(IIV_CORR) sim.orig$y <- sim.orig$OBS  - sim.orig$IPRED + sim.orig$PRED

  # to omit
  if(!IIV_CORR) obs.orig$y <- obs.orig$OBS
  if(!IIV_CORR) sim.orig$y <- sim.orig$OBS

  list(
    model.name = model.name,
    indivsam.obs = indivsam.obs,
    output.typ = output.typ,
    indivsam.vpc = indivsam.vpc,
    vachette.covs = vachette.covs,
    VVPC = VVPC,
    ADD_TR = ADD_TR,
    PROP_TR = PROP_TR,
    output.typ = output.typ,
    obs.orig = obs.orig,
    sim.orig = sim.orig,
    tab.ucov = tab.ucov,
    ref.region = ref.region,
    ref.dose = ref.dose,
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

.define_and_enumerate_regions <- function(output.typ, obs.orig, sim.orig, vachette.covs, ref.dose, ref.region) {
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
  comb.ucov <- output.typ %>% tidyr::expand(!!!syms(names(vachette.covs)))

  # Collect unique covariate/dose combinations and define region.type
  # Add to info to data frames
  if(length(vachette.covs)==1)
  {
    for(i in c(1:dim(comb.ucov)[1])) #e.g., unique values of covariates nrow(comb.ucov)
    {
      # Get number of doses for each covariate combination:
      z <- output.typ %>% filter(!!sym(names(vachette.covs)) == comb.ucov[names(vachette.covs)][[1]][i])
      ndose <- length(unique(z$dosenr))
      for(idose in c(1:ndose))
      {
        if (idose == ndose) region.type <- 'open'   # Last dose
        if (idose <  ndose) region.type <- 'closed' # Any dose before last dose

        # New unique combination
        n.ucov <- n.ucov + 1
        add <- data.frame(ucov = n.ucov) %>%
          mutate(!!names(vachette.covs) := comb.ucov[names(vachette.covs)][[1]][i],
                 region = idose,
                 region.type = region.type)
        tab.ucov <- rbind(tab.ucov, add)

        # Selection of new unique combination
        sel.typ  <- output.typ[[names(vachette.covs)]] == comb.ucov[names(vachette.covs)][[1]][i] &
          output.typ$dosenr == idose
        sel.obs  <- obs.orig[[names(vachette.covs)]] == comb.ucov[names(vachette.covs)][[1]][i] &
          obs.orig$dosenr == idose
        sel.sim  <- sim.orig[[names(vachette.covs)]] == comb.ucov[names(vachette.covs)][[1]][i] &
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
  if(length(vachette.covs) > 1)
  {
    vars <- names(vachette.covs)

    for(i in c(1:dim(comb.ucov)[1]))
    {
      filter_query <- paste0(vars, "==", "comb.ucov$", vars, "[i]", collapse = " & ")
      # Get number of doses for covariate combination:
      z <- output.typ %>% filter(eval(parse(text = filter_query)))

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
          paste0("output.typ$", vars, "==", "comb.ucov$", vars, "[i]", collapse = " & "),
          " & output.typ$dosenr == idose"
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
        # sel.typ  <- output.typ$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        #   output.typ$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        #   output.typ$dosenr == idose
        # sel.obs  <- obs.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        #   obs.orig$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        #   obs.orig$dosenr == idose
        # sel.sim  <- sim.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        #   sim.orig$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        #   sim.orig$dosenr == idose

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

  tab.ucov$ref <- "No"

  for(i.ucov in c(1:dim(tab.ucov)[1]))
  {
    # if(length(vachette.covs)==1)
    #   if(tab.ucov[i.ucov, vachette.covs] == ref.cov1 && tab.ucov[i.ucov, "region"] == ref.region)
    #     tab.ucov[i.ucov, "ref"] <- "Yes"
    # # where both covs equal ref values in table, set yes
    # if(length(vachette.covs)>1)
      condition <- paste0(
      paste0("as.character(tab.ucov$", names(vachette.covs), "[i.ucov])", "==", "'", vachette.covs, "'", collapse = " && "),
      " & tab.ucov$region[i.ucov] == ref.region"
    )
      if(eval(parse(text = condition)))
        tab.ucov$ref[i.ucov] <- "Yes"
  }

  # Add reference flag to output.typ
  output.typ <- output.typ %>%
    dplyr::left_join(tab.ucov[,c('ucov','ref')],by='ucov')

  return(
    list(
      output.typ = output.typ,
      obs.orig = obs.orig,
      sim.orig = sim.orig,
      tab.ucov = tab.ucov,
      n.ucov = n.ucov
    )
  )
}


#' Apply vachette transformations
#'
#' @param vachette_data Object of class \code{vachette_data}
#' @param LM_REFINE Logical; Set to \code{TRUE} to refine landmark
#' @param tolend Numeric;
#' @param tolnoise Numeric;
#' @param step.x.factor Numeric;
#' @param ngrid.open.end Numeric;
#' @param w.init Numeric; Savitzky Golay smoothing.
#' @param w1.refine Numeric; Savitzky Golay smoothing, first derivative.
#' @param w2.refine Numeric; Savitzky Golay smoothing, second derivative.
#'
#' @name apply_transformations
#' @export
apply_transformations <- function(vachette_data, ...) UseMethod("apply_transformations")

#' @rdname apply_transformations
#' @export
apply_transformations.vachette_data <-
  function(vachette_data,
           LM_REFINE = FALSE,
           tolend = 0.001,
           tolnoise = 1e-8,
           step.x.factor = 1.5,
           ngrid.open.end = 100,
           w.init = 17,     # Savitzky Golay smoothing - initial landmarks
           w1.refine = 7,     # Savitzky Golay smoothing - refine landmarks - first derivative
           w2.refine = 5) {


  stopifnot(inherits(vachette_data, "vachette_data"))
    # Collect all Vachette query curves and original/transformed observations (incl reference)
    curves.all           <- NULL
    curves.scaled.all    <- NULL
    obs.all              <- NULL

    my.ref.lm.all   <- NULL
  my.query.lm.all <- NULL
  lm.all <- NULL

  tab.ucov <- vachette_data$tab.ucov
  stopifnot(!is.null(tab.ucov))
  output.typ <- vachette_data$output.typ
  stopifnot(!is.null(output.typ))
  obs.orig <- vachette_data$obs.orig
  stopifnot(!is.null(obs.orig))
  sim.orig <- vachette_data$sim.orig
  stopifnot(!is.null(sim.orig))

  vachette.covs <- vachette_data$vachette.covs
  ref.region <- vachette_data$ref.region
  ref.dose <- vachette_data$ref.dose
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
          names(vachette.covs),
          ')',
          "==",
          "'",
          vachette.covs,
          "'",
          collapse = " & "
        ),
        " & region == ref.region"
      )

    ref <- output.typ %>% filter(eval(parse(text = filter_query)))

    ref.ucov <- unique(ref$ucov)

    if(length(ref.ucov) != 1) stop("Error: length ref.ucov != 1")     # A single ref only possible

    ref.region.type <- tab.ucov$region.type[tab.ucov$ucov==ref.ucov]
    # ref = Typical reference curve
    query.region       <- tab.ucov$region[tab.ucov$ucov==i.ucov]
    query.region.type  <- tab.ucov$region.type[tab.ucov$ucov==i.ucov]

    # Typical query curve
    query <- output.typ %>%
      filter(ucov == i.ucov) %>%
      mutate(ref = tab.ucov$ref[i.ucov])    # Flag for reference

    # 230418 - extrapolate region last-x region by one gridstep size if region.type = 'closed'
    # Currently simple extra x value with same y value ("horizontal" extrapolation - LOCF)
    if(ref.region.type == 'closed')
    {
      ref.grid.step.size = ref$x[dim(ref)[1]] - ref$x[dim(ref)[1]-1]
      ref.curve.last.point <- ref[dim(ref)[1],]
      # Change x only, assume change of y is negligible
      ref.curve.last.point$x <- ref.curve.last.point$x + ref.grid.step.size
      # Add
      ref <- rbind(ref,ref.curve.last.point)
    }
    if(query.region.type == 'closed')
    {
      query.grid.step.size = query$x[dim(query)[1]] - query$x[dim(query)[1]-1]
      query.curve.last.point <- query[dim(query)[1],]
      # Change x only, assume change of y is negligible
      query.curve.last.point$x <- query.curve.last.point$x + query.grid.step.size
      # Add
      query <- rbind(query,query.curve.last.point)
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

    #  B. -------- Get landmarks approximate positions ----

    # Tolerance to apply for finding maximums, minimums and inflection points
    # Small stepsize -> small tolerance
    # Large stepsize -> large tol
    #tolapply = distance between second and first x obs * tolnoise
    tolapply <- tolnoise*(output.typ$x[2]-output.typ$x[1])
    # get ref landmarks
    my.ref.lm.init    <- get.x.multi.landmarks(ref$x,ref$y,w=w.init,tol=tolapply) #contiguous inflec cause issue?
    my.ref.lm.init$y  <- approx(ref$x,ref$y, xout=my.ref.lm.init$x)$y #interpolation
    # get query landmarks
    my.query.lm.init    <- get.x.multi.landmarks(query$x,query$y,w=w.init,tol=tolapply)
    my.query.lm.init$y  <- approx(query$x,query$y, xout=my.query.lm.init$x)$y

    # stop at this point, perhaps, allow user to visually inspect this plot, to assess what increase in stepsize
    # do not error out, return data that is available, and provide to user for plotting, and provide suggestions for updating argument values

    #Validation check
    #If y contiguous y values in series are the same, provide error, increase tolnoise
    #option 1: if multiple inflec in sequence, try decreasing tolnoise

    if(!LM_REFINE)
    {
      my.ref.lm.refined      <- my.ref.lm.init
      my.query.lm.refined    <- my.query.lm.init
    }
    if(LM_REFINE)
    {
      my.ref.lm.refined    <- refine.x.multi.landmarks(x=ref$x,y=ref$y,lm=my.ref.lm.init,tol=0.01*tolapply, w1=w1.refine,w2=w2.refine)
      my.ref.lm.refined$y  <- approx(ref$x, ref$y, xout=my.ref.lm.refined$x)$y

      my.query.lm.refined    <- refine.x.multi.landmarks(query$x,query$y,lm=my.query.lm.init,tol=0.01*tolapply, w1=w1.refine,w2=w2.refine)
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
      my.ref.lm   <- get.ref.x.open.end(ref$x,ref$y,my.ref.lm.refined,step.x.factor=step.x.factor,tol=tolend)
      my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.open.end,scaling=scaling)
    } #error is in above line
    # 2. Fit ref last x, Fix query last x
    if(ref.region.type == 'open' & query.region.type == 'closed')
    {
      my.ref.lm   <- get.query.x.open.end(query,ref,my.query.lm.refined,my.ref.lm.refined,ngrid=ngrid.open.end,scaling=scaling)
      my.query.lm <- my.query.lm.refined
    }
    # 230418 - Update
    # 3+4. Fix ref last x, Fit query last x
    if(ref.region.type == 'closed')
    {
      my.ref.lm   <- my.ref.lm.refined
      my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.open.end,scaling=scaling)
    }

    # @James recent developments
    # Overrides above (above needs reprogramming)
    # Always FIX ref and FIT query to find best fitting x ("open"). Then extrapolate if required
    if(ref.region.type == 'closed')
    {
      my.ref.lm   <- my.ref.lm.refined
      my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.open.end,scaling=scaling)
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

        # @James: I have to check this part of the code using a suitable example

        EXTENSION <- TRUE

        message('--------------------------------------------------------')
        message('EXTENSION REF CURVE UPDATED >> further checking required')
        message('--------------------------------------------------------')

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
        lastpoint        <- dim(ref)[1]
        add.x <- max.obs.x.scaled
        add.y <- approxExtrap(c(ref$x[lastpoint-1],ref$x[lastpoint]),
                              c(ref$y[lastpoint-1],ref$y[lastpoint]),
                              xout = add.x)$y

        cur.ref.seg   <- ref$seg[dim(ref[!is.na(ref$seg),])[1]]
        ref.add <- ref %>%
          slice(1) %>%
          mutate(x=add.x,
                 y=add.y) %>%
          mutate(seg=cur.ref.seg)
        ref <- rbind(ref,ref.add)

        # Check:
        if(cur.ref.seg != cur.query.seg) stop("Error segment number assignment in reference extrapolation block")

        print("**** END EXTENSION *****")
      }
    }

    # Collect all typical curves with scaling factors and scaled x,y values
    query.scaled$y.scaled <- approx(ref$x,ref$y,xout=query.scaled$x.scaled)$y
    curves.all            <- rbind(curves.all,query)
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
    typ.curve <- output.typ %>% filter(ucov==i.ucov)
    ref.curve <- output.typ %>% filter(ucov==ref.ucov)
    obs       <- obs.all    %>% filter(ucov==i.ucov)

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

  update(vachette_data, obs.all=obs.all, curves.all = curves.all, lm.all = lm.all, nseg = nseg)


  }

`%notin%` <- Negate(`%in%`)

.validate_columns <- function(mappings, indivsam.obs, output.typ, IIV_CORR) {

  obs_req_cols <- c("ID", "PRED", "OBS", "x") #, "dosenr")
  if (isTRUE(IIV_CORR)) {
    obs_req_cols <- c(obs_req_cols, "IPRED")
  }
  sim_req_cols <- c("ID", "PRED", "x") #, "dosenr")

  obs_cols <- colnames(indivsam.obs)
  sim_cols <- colnames(output.typ)

  if (is.null(mappings)) {
    if (any(obs_req_cols %notin% obs_cols)) {
      missing_cols <- setdiff(obs_req_cols, obs_cols)
      stop(
        paste0(missing_cols, collapse = " "),
        " column(s) not found in observed data. Use the 'mappings' argument to manually specify columns mappings."
      , call. = FALSE)
    }

    if (any(sim_req_cols %notin%  sim_cols)) {
      missing_cols <- setdiff(sim_req_cols, sim_cols)
      stop(
        paste0(missing_cols, collapse = " "),
        " column(s) not found in simulated data. Use the 'mappings' argument to manually specify columns mappings."
        , call. = FALSE)
    }
  } else {
    diff_req_mappings <- setdiff(obs_req_cols, names(mappings))
    cols_not_mapped <- setdiff(diff_req_mappings, obs_cols)

    if (length(cols_not_mapped) > 0) {
      stop(
        paste0(cols_not_mapped, collapse = " "),
        " column(s) not found in observed data and not provided in 'mappings'. Use the 'mappings' argument to manually specify columns mappings."
        , call. = FALSE)
      }

    if (any(mappings %notin% obs_cols)) {
      missing_cols <- setdiff(mappings, obs_cols)
      stop(
        paste0(missing_cols, collapse = " "),
        " column(s) provided in mappings argument are not found in observed data."
        , call. = FALSE)
    }

    sim_mappings <- mappings[names(mappings) %notin% c("IPRED", "OBS")]

    if (any(sim_mappings %notin% sim_cols)) {
      missing_cols <- setdiff(sim_mappings, sim_cols)
      stop(
        paste0(missing_cols, collapse = " "),
        " column(s) provided in mappings argument are not found in simulated data."
        , call. = FALSE)
    }

    indivsam.obs <- dplyr::rename(indivsam.obs, dplyr::all_of(mappings))
    output.typ <- dplyr::rename(output.typ, dplyr::all_of(sim_mappings))
  }
  return(
    list(
      indivsam.obs = indivsam.obs,
      output.typ = output.typ
    )
  )
}

.process_covariates <- function(vachette.covs, indivsam.obs) {
  if (is.null(names(vachette.covs))) {
    names(vachette.covs) <- vachette.covs
    cov_names <- vachette.covs
  } else {
    cov_names <- names(vachette.covs)
  }
  for (i in seq_along(cov_names)) {
    #browser()
    cov_name <- cov_names[i]
    cov_ref_val <- vachette.covs[i]
    # If name not supplied, expect value as name, then assign name
    if (cov_name == "") {
      cov_name <- cov_ref_val
      names(vachette.covs)[i] <- cov_name
    }
    stopifnot(cov_name %in% colnames(indivsam.obs))
    # Automatically set ref value to median/mode given cont/cat covariate type, if ref value not given
    # Account for case median is used with even number of values in dataset
    # - compute median for each subject first, then compute median from there

    if (cov_name == cov_ref_val) {
      if (suppressWarnings(all(is.na(as.numeric(unlist(indivsam.obs[, cov_name])))))) {
        cov_ref_val <- .mode(unlist(indivsam.obs[, cov_name]))
      } else {
        if (nrow(indivsam.obs) %% 2 == 0) {
          #indivsam.obs <- dplyr::slice(indivsam.obs, -1) #before midpoint, drop first value
          indivsam.obs <- dplyr::slice(indivsam.obs, -nrow(indivsam.obs)) #after midpoint, drop last value
        }
        cov_ref_val <- median(unlist(indivsam.obs[, cov_name]))
      }
      stopifnot(cov_ref_val %in% as.character(unlist(indivsam.obs[, cov_name])))
      vachette.covs[i] <- cov_ref_val
    } #else add check for provided cov value, ensuring exists in simulated data
  }
  return(vachette.covs)
}


.calculate_dose_number <- function(data, ref.dosenr, data_type = c("obs.data", "typ.data")) {

  data_type <- match.arg(data_type)

  # If dose number provided, and column in data, use that
  if ("dosenr" %in% colnames(data)) {
    message("`dosenr` column found in ", data_type, ", using `dosenr` column in data for corresponding ref.dosenr value")
  } else {
    if ("EVID" %in% colnames(data)) {
      message("`EVID` column found in ", data_type, ", creating `dosenr` column in data for corresponding ref.dosenr value")
      data <- data %>%
        group_by(ID) %>%
        mutate(dosenr = cumsum(EVID == 1)) %>%
        ungroup()
    } else if (all(c("ADDL", "II", "AMT") %in% colnames(data))) {
      message("`ADDL`, `II`, `AMT` columns found in ", data_type, ", creating `dosenr` column in data for corresponding ref.dosenr value")
      dose_data <- data %>%
        select(x, ID, AMT, ADDL, II) %>%
        filter(ADDL > 0)
      dosing <- list()
      for (i in 1:nrow(dose_data)) {
        dose_data_row <- slice(dose_data, i)
        dosing[[i]] <- data.frame(x = seq(from = dose_data_row$x,
                                          by = as.numeric(dose_data_row$II),
                                          length.out = as.numeric(dose_data_row$ADDL) + 1),
                                  ID = dose_data_row$ID,
                                  AMT = dose_data_row$AMT
        )
      }
      dose_data_expanded <- do.call(rbind, dosing)
      data <- filter(data, ADDL == 0) %>%
        select(-ADDL, -II)

      data_new <- bind_rows(data, dose_data_expanded) %>%
        arrange(ID, x)

      data <- data_new %>%
        group_by(ID) %>%
        mutate(dosenr = cumsum(AMT != 0)) %>%
        ungroup()
    } else if ("AMT" %in% colnames(data)) {
      data <- data %>%
        group_by(ID) %>%
        mutate(dosenr = cumsum(AMT != 0)) %>% #if values NA: cumsum(!is.na(amt))
        ungroup()
    } # final control flow should be understood as no dose information in the data, no ref.dosnr, message such e.g., sigmoid model
  }

  # ensure value supplied to ref.dosenr exists inside dosenr column
  if (ref.dosenr %notin% unique(data[["dosenr"]])) {
    stop("ref.dosenr value of ", ref.dosenr, " not found in dosenr column in ", data_type)
  }

  return(data)
}


.mode <- function(x){
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}
