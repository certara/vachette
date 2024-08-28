
# For same looks and feel of plots:
render <- theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        plot.title = element_text(size=14),
        plot.subtitle=element_text(size=12),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0, vjust = 1, size = 10)
  )


#' Plot \code{vachette} Transformed Typical Curves with Landmarks
#'
#' This function generates a ggplot2 visualization of the \code{vachette} transformed
#' typical curves and their corresponding landmarks for a given pharmacometric model.
#' The plot highlights the segments of the reference and query curves, with a special emphasis
#' on the landmark positions.
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing the \code{vachette} transformed typical curves with
#' their landmarks. The plot displays different segments of the curves in different color,
#' and landmarks marked with black crosses.
#'
#' @details The function plots the typical curves along with their
#' segments and marking the landmarks. If the x-axis is
#' logarithmic, it adjusts the axis label accordingly.
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used. The x-axis range is dynamically set
#' based on the minimum and maximum x values, before and after scaling.
#'
#' @export
p.scaled.typical.curves.landmarks <- function(vachette_data) {

  x <- y <- ucov <- ref <- x.scaled <- y.scaled <- seg <- NULL
  stopifnot(inherits(vachette_data, "vachette_data"))
  log.x <- vachette_data$log.x

  curves.all <- vachette_data$curves.all
  lm.all <- vachette_data$lm.all
  model.name <- vachette_data$model.name
  xstart <- min(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)

  stopifnot(!is.null(curves.all))

  gg <- curves.all %>%
    ggplot(aes(x=x,y=y,group=ucov,
               col = factor(seg)))+
    geom_line(lwd=1)+
    geom_line(data=curves.all %>% filter(ref=="Yes"),col='black',lty=2,lwd=1,alpha=0.5)+
    # Add landmark positions
    geom_point(data=lm.all,pch=3,size=4,col='black',stroke = 2)

  gg <- gg + coord_cartesian(xlim=c(xstart,xstop))

  gg <- gg +
    labs(title=paste0(model.name,"; Typical curve segments"),
         subtitle = "Dashed line: Reference typical curve\nGrey: Unused part of typical curve",
         caption = paste0("Reference Covariate: ",
                          paste0(
                            names(vachette_data$covariates),
                            "=",
                            vachette_data$covariates,
                            collapse = ", "
                          )),
         col="Segment")

  if (log.x) {
    gg <- gg +
      labs(x = "ln(x)")
  }

  gg <- gg + render

  return(gg)
}

#' Plot Full \code{vachette} Transformed Typical Curves with Landmarks
#'
#' This function generates a ggplot2 visualization of the full \code{vachette} transformed typical
#' curves and their corresponding landmarks for a given pharmacometric model.
#' The plot highlights the segments of the reference and query curves, with a special emphasis on
#' specific landmark positions.
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing the full \code{vachette} transformed typical curves
#' with their landmarks. The plot displays different segments of the curves in different colors,
#' and landmarks marked with black crosses.
#'
#' @details The function plots the full typical curves along with
#' their segments, applying specific styling to the reference curve and marking
#' the landmarks. If the x-axis is logarithmic, it adjusts the axis label accordingly.
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used. The x-axis range is dynamically set
#' based on the minimum and maximum x values, before and after scaling.
#'
#' @export
p.scaled.typical.full.curves.landmarks <- function(vachette_data) {

  x <- y <- ucov <- ref <- x.scaled <- y.scaled <- seg <- NULL
  stopifnot(inherits(vachette_data, "vachette_data"))
  log.x <- vachette_data$log.x

  curves.all <- vachette_data$curves.all
  lm.all <- vachette_data$lm.all
  model.name <- vachette_data$model.name
  xstart <- min(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)

  gg <- curves.all %>%
    ggplot(aes(x=x,y=y,group=ucov,
               col = factor(seg)))+
    geom_line(lwd=1)+
    geom_line(data=curves.all %>% filter(ref=="Yes"),col='black',lty=2,lwd=1,alpha=0.5)+
    # Add landmark positions
    geom_point(data=lm.all,pch=3,size=4,col='black',stroke = 2)+
    labs(title=paste0(model.name,"; Typical curve segments"),
         subtitle = "Dashed: Reference typical curve\nGrey: Unused part of typical curve",
         caption = paste0("Reference Covariate: ",
                          paste0(
                            names(vachette_data$covariates),
                            "=",
                            vachette_data$covariates,
                            collapse = ", "
                          )),
         col="Segment")

  if (log.x) {
    gg <- gg +
      labs(x = "ln(x)")
  }

  gg <- gg +
    render

  return(gg)
}

#' Plot Scaling Factors for Typical Curves
#'
#' This function generates a ggplot2 visualization of the scaling factors applied
#' to the typical curves within a pharmacometric model. The plot differentiates between
#' reference and query curves and illustrates how scaling factors vary across
#' segments.
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing the scaling factors for the curves.
#' The plot includes facets for each curve type (reference or query),
#' and it shows the scaling factors as a function of the x-axis values.
#'
#' @details The function creates a plot that displays
#' the x-scaling factors for the curves, with different segments color-coded.
#' The plot is faceted by curve type (reference or query).
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used. The x-axis range is dynamically set
#' based on the minimum and maximum x values, before and after scaling.
#'
#' @export
p.scaling.factor <- function(vachette_data) {
  ref <- x <- x.scaling <- seg <- NULL
  curves.scaled.all <- vachette_data$curves.scaled.all
  obs.all           <- vachette_data$obs.all
  model.name        <- vachette_data$model.name
  xstart <- min(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)

  curves.scaled.all %>%
    mutate(mycurve = ifelse(ref=='Yes','Reference','Query')) %>%
    ggplot(aes(x=x,y=x.scaling,col=factor(seg)))+
    geom_line(lwd=1)+
    facet_wrap(~paste(mycurve," covariate =",COV))+
    coord_cartesian(xlim=c(xstart,xstop),
                    ylim=c(0,max(curves.scaled.all$x.scaling[curves.scaled.all$x<=max(obs.all$x,obs.all$x.scaled)])))+

    labs(title=paste0(model.name,"; x-scaling factors"),
         caption = paste0("Reference Covariate: ",
                          paste0(
                            names(vachette_data$covariates),
                            "=",
                            vachette_data$covariates,
                            collapse = ", "
                          )),
         col="Segment") +
    render
}

#' Plot \code{vachette} Transformed Typical Curves
#'
#' This function generates a ggplot2 visualization of the \code{vachette} transformed typical curves
#' for a given pharmacometric model. It distinguishes between reference curves
#' and query curves, showing how the curves transform.
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing the \code{vachette} transformed typical curves. The plot
#' shows both the reference and query curves, and the curves after
#' the vachette transformation.
#'
#' @details The function creates a plot that distinguishes between
#' reference curves (in red) and query curves (in blue). The dashed
#' lines represent the curves after the vachette transformation.
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used. The x-axis range is dynamically set
#' based on the minimum and maximum x values, before and after scaling.
#'
#' @export
p.scaled.typical.curves <- function(vachette_data) {

  x <- y <- ucov <- ref <- x.scaled <- y.scaled <- seg <- NULL
  stopifnot(inherits(vachette_data, "vachette_data"))
  log.x <- vachette_data$log.x

  curves.scaled.all <- vachette_data$curves.scaled.all
  model.name        <- vachette_data$model.name
  xstart <- min(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)

  gg <- curves.scaled.all %>%
    ggplot(aes(x=x,y=y,group=ucov))+
    # geom_line(lwd=1.5,alpha=0.5)+
    geom_line(data=curves.scaled.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Yes'),lwd=1.5,alpha=0.60)+
    geom_line(data=curves.scaled.all %>% filter(ref=='No'),aes(x=x,y=y,col='No'),lwd=1.5)+
    geom_line(aes(x=x.scaled, y=y.scaled), col='black',lwd=0.75,lty=2)
  # scale_color_manual(values=c('red','blue')) +
  # scale_linetype_manual(values = c(3, 1),
  #                       labels = c("Covariate", "Reference")) +

  gg <- gg + coord_cartesian(xlim=c(min(curves.scaled.all$x),xstop))

  gg <- gg +
    scale_color_manual(name='Reference\ncurve',
                       breaks=c('No', 'Yes'),
                       values=c('No'='blue',
                                'Yes'='red'))+
    labs(title=paste0(model.name,"; Covariate typical curves"),
         subtitle = "Dashed: Covariate typical curves after Vachette transformation",
         caption = paste0("Reference Covariate: ",
                          paste0(
                            names(vachette_data$covariates),
                            "=",
                            vachette_data$covariates,
                            collapse = ", "
                          )),
         col="Covariate value\n(Reference)")

  if (log.x) {
    gg <- gg +
      labs(x = "ln(x)")
  }

  gg <- gg +
    render

  return(gg)
}

#' Plot \code{vachette} Transformed Observation Curves
#'
#' This function generates a ggplot2 visualization of transformed observation curves,
#' including their transformations and overlays with reference curves. The plot is
#' useful for comparing the original and the \code{vachette}
#' transformed observation data against typical reference curves.
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing the \code{vachette} transformed observation curves. The plot
#' includes different layers for query observations, their transformed counterparts,
#' and the reference curves, with appropriate color coding and line styles.
#'
#' @details The function constructs a plot that overlays
#' observation curves with the corresponding reference curves. Query observations
#' and their transformed curves are highlighted, and if applicable, an extrapolated
#' segment of the reference curve is shown as a dashed line.
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used. The x-axis range is dynamically set
#' based on the minimum and maximum x values, before and after scaling.
#'
#' @export
p.scaled.observation.curves <- function(vachette_data) {

  x <- y <- ID <- COV <- ref <- x.scaled <- y.scaled <- ucov <- NULL
  stopifnot(inherits(vachette_data, "vachette_data"))
  log.x <- vachette_data$log.x

  # Copy ref curves to each ID
  obs.all <- vachette_data$obs.all
  myids  <- unique(obs.all$ID)
  curves.all <- vachette_data$curves.all
  curves.scaled.all <- vachette_data$curves.scaled.all
  #curves.all.ids <- NULL
  xstart <- min(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)

  # JL 230607
  ref.extensions.all <- vachette_data$ref.extensions.all
  # Plot longest ref.extension.all only. Pick first is multiple occurences
  if(!is.null(ref.extensions.all)) max.x.ucov <- ref.extensions.all$ucov[ref.extensions.all$x == max(ref.extensions.all$x)][1]

  #for(iid in c(1:length(myids))) {curves.all.ids <- rbind(curves.all.ids,curves.all %>% mutate(ID=myids[iid]))}
  curves.all.ids <- purrr::map_dfr(myids, ~curves.all %>% mutate(ID = .x))

  # Overlay observation curves per ID
  gg <- obs.all %>%
    ggplot(aes(x=x,y=y,group=paste(ID,COV)))+

    geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
    geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'),lty=2)+
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'))+

    # Reference in top layer
    geom_line(data=curves.all.ids %>% filter(ref=="Yes"),aes(x=x,y=y,col='Typical reference'),lwd=1)+
    geom_line(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference observations'))

  if(!is.null(ref.extensions.all))
  {
    gg <- gg + geom_line(
      data = ref.extensions.all %>% filter(ucov==max.x.ucov),
      aes(x = x, y = y, col = 'Typical reference'), lty=2,
      lwd = 0.8, col='grey30'
    )
  }

  gg <- gg + coord_cartesian(xlim=c(min(curves.scaled.all$x),xstop))

  gg <- gg +
    geom_line(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference observations'))+
    geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference observations'))+
    # scale_color_manual(values=c('red','blue')) +
    # scale_linetype_manual(values = c(3, 1),
    #                       labels = c("Covariate", "Reference")) +
    scale_color_manual(name='Data type',
                       breaks=c('Query transformed',
                                'Query observations',
                                'Reference observations',
                                'Typical reference'),
                       values=c('Query transformed'='purple',
                                'Query observations'='blue',
                                'Reference observations'='red',
                                'Typical reference'='black'))

  gg <- gg +
    labs(title=paste0(vachette_data$model.name,"; Individual observation curves"),
         subtitle = "Dashed: Extrapolation of typical curve",
         caption = paste0("Reference Covariate: ",
                          paste0(
                            names(vachette_data$covariates),
                            "=",
                            vachette_data$covariates,
                            collapse = ", "
                          )),
         col="Covariate value\n(Reference)")

  if (log.x) {
    gg <- gg +
      labs(x = "ln(x)")
  }

  gg <- gg +
    render

  return(gg)
}

#' Plot \code{vachette} Transformed Observation Curves by Individual
#'
#' This function generates a ggplot2 visualization of transformed observation curves
#' facetted by individual, including their transformations and overlays
#' with reference curves. The plot is useful for comparing the original and the \code{vachette}
#' transformed observation data against typical reference curves.
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing the \code{vachette} transformed observation curves. The plot
#' includes different layers for query observations, their transformed counterparts,
#' and the reference curves, with appropriate color coding and line styles.
#'
#' @details The function constructs a plot that overlays individual
#' observation curves with the corresponding reference curves. Query observations
#' and their transformed curves are highlighted, and if applicable, an extrapolated
#' segment of the reference curve is shown as a dashed line.
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used. The x-axis range is dynamically set
#' based on the minimum and maximum x values, before and after scaling.
#'
#' @export
p.scaled.observation.curves.by.id <- function(vachette_data) {

  x <- y <- ID <- COV <- ref <- x.scaled <- y.scaled <- ucov <- NULL
  stopifnot(inherits(vachette_data, "vachette_data"))
  log.x <- vachette_data$log.x

  # Observation curve for each ID
  obs.all <- vachette_data$obs.all
  myids  <- unique(obs.all$ID)
  curves.all <- vachette_data$curves.all
  curves.all.ids <- NULL
  xstart <- min(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)

  # JL 230607
  ref.extensions.all <- vachette_data$ref.extensions.all
  # Plot longest ref.extension.all only. Pick first is multiple occurences
  if(!is.null(ref.extensions.all)) max.x.ucov <- ref.extensions.all$ucov[ref.extensions.all$x == max(ref.extensions.all$x)][1]

  for(iid in c(1:length(myids))) {curves.all.ids <- rbind(curves.all.ids,curves.all %>% mutate(ID=myids[iid]))}
  gg <- obs.all %>%
    ggplot(aes(x=x,y=y,group=paste(ID,COV)))+

    geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
    geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'),lty=2)+
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'))+

    geom_line(data=curves.all.ids %>% filter(ref=="Yes"),aes(x=x,y=y,col='Typical reference'),lwd=1)+
    geom_line(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference observations'))
  if(!is.null(ref.extensions.all))
  {
    gg <- gg + geom_line(
      data = ref.extensions.all %>% filter(ucov==max.x.ucov),
      aes(x = x, y = y, col = 'Typical reference'), lty=2,
      lwd = 0.8, col='grey30'
    )
  }
  gg <- gg +
    geom_line(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference observations'))+
    geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference observations'))+
    # scale_color_manual(values=c('red','blue')) +
    # scale_linetype_manual(values = c(3, 1),
    #                       labels = c("Covariate", "Reference")) +
    scale_color_manual(name='Data type',
                       breaks=c('Query transformed',
                                'Query observations',
                                'Reference observations',
                                'Typical reference'),
                       values=c('Query transformed'='purple',
                                'Query observations'='blue',
                                'Reference observations'='red',
                                'Typical reference'='black'))+
    facet_wrap(~ID)

  gg <- gg + coord_cartesian(xlim=c(xstart,xstop))

  gg <- gg +
    labs(title=paste0(vachette_data$model.name,"; Individual observation curves by ID"),
         subtitle = "Dashed: Extrapolation of typical curve",
         caption = paste0("Reference Covariate: ",
                          paste0(
                            names(vachette_data$covariates),
                            "=",
                            vachette_data$covariates,
                            collapse = ", "
                          )),
         col="Covariate value\n(Reference)")

  if (log.x) {
    gg <- gg +
      labs(x = "ln(x)")
  }

  gg <- gg +
    render

  return(gg)

}


#' Plot Distances Before and After \code{vachette} Transformation Assuming Additive Error Model
#'
#' This function generates a ggplot2 visualization of the distances
#' before and after a \code{vachette} transformation between observations and the
#' relative typical curves.
#' It is useful for assessing the impact of additive error
#' model assumption on the distances.
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing the plot of distances.
#' The plot compares the original distances to the transformed distances, with
#' segments color-coded.
#'
#' @details The function starts by setting the x-axis range based on the minimum
#' and maximum values of the original and scaled distances. It then plots the
#' distances before and after transformation, with a reference line
#' indicating where the original and transformed distances are equal.
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used.
#'
#' @export
p.add.distances <- function(vachette_data) {
  dist.add.orig <- dist.add.transformed <- seg <- NULL
  xstart <- min(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)

  obs.all <- vachette_data$obs.all

  # Plot additive distances before and after transformation
  gg <- obs.all %>%
    ggplot(aes(x = dist.add.orig, y = dist.add.transformed, col = factor(seg))) +
    geom_abline(slope = 1) +
    geom_point() +
    labs(title = paste0(vachette_data$model.name, "; Normal distances: original and after transformation"),
         subtitle = paste0(if (vachette_data$ADD_TR) "Additive Error Transformation",
                           if (vachette_data$PROP_TR) "Proportional Error Transformation"),
         caption = paste0("Reference Covariate: ",
                          paste0(
                            names(vachette_data$covariates),
                            "=",
                            vachette_data$covariates,
                            collapse = ", "
                          )),
         x = 'Original distance',
         y = 'Distance after transformation',
         col = "Segment") +
    render

  return(gg)
}


#' Plot Distances Before and After \code{vachette} Transformation Assuming Proportional Error Model
#'
#' This function creates a ggplot2 visualization of the distances
#' before and after the \code{vachette} transformation between observations
#' and the relative typical curves. It is useful for assessing the impact of
#' proportional error model assumption on the distances.
#'
#' @param vachette_data An object of class \code{vachette_data}, which contains
#' the necessary data for plotting.
#'
#' @return A \code{ggplot2} object representing the plot of distances.
#' The plot compares the log of the original distances to the log of the distances
#' after transformation, color-coded by segment.
#'
#' @details The plot includes a reference line with a slope of 1, indicating where
#' the original and transformed distances are equal. The title of the plot is
#' dynamically generated based on the model name provided in the \code{vachette_data} object.
#'
#' The subtitle of the plot indicates whether an additive or proportional error
#' transformation was applied. The caption provides the reference covariate(s)
#' used in the model.
#'
#' @export
#'
p.prop.distances <- function(vachette_data) {
  dist.prop.orig <- dist.prop.transformed <- seg <- NULL
  obs.all <- vachette_data$obs.all
  xstart <- min(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)

  # Proportional distances before and after transformation
  obs.all %>%
    ggplot(aes(x=dist.prop.orig,y=dist.prop.transformed,col=factor(seg)))+
    geom_abline(slope=1)+
    geom_point()+
    labs(title=paste0(vachette_data$model.name,"; Proportional distances: original and after transformation"),
         subtitle = paste0(if(vachette_data$ADD_TR) "Additive Error Transformation",if(vachette_data$PROP_TR) "Proportional Error Transformation"),
         caption = paste0("Reference Covariate: ",
                          paste0(
                            names(vachette_data$covariates),
                            "=",
                            vachette_data$covariates,
                            collapse = ", "
                          )),
         x = 'Log of original distances',
         y = 'Log of vachette-transformed distances',
         col="Segm.") +
    render

}

#' Plot Observations and Typical Curves for Query and Reference Data
#'
#' This function generates a ggplot2 visualization of the observations and typical curves
#' for both query and reference data within a pharmacometric model. It is designed to
#' compare the observed data against the typical reference and query curves.
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing the observations and typical curves.
#' The plot displays the query and reference curves and their corresponding observations,
#' color-coded for easy comparison.
#'
#' @details The function plots the query and reference data,
#' differentiating them using distinct colors. Observations are represented
#' by points, while the typical curves are depicted as lines.
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used. The x-axis range is dynamically set
#' based on the minimum and maximum x values, before and after scaling.
#' If the x-axis is logarithmic, the axis label is adjusted accordingly.
#'
#' @export
p.obs.ref.query <- function(vachette_data) {

  x <- y <- ref <- ucov <- NULL
  stopifnot(inherits(vachette_data, "vachette_data"))
  log.x <- vachette_data$log.x

  obs.all <- vachette_data$obs.all
  curves.all <- vachette_data$curves.all
  xstart <- min(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)

  gg <- obs.all %>%
    ggplot(aes(x=x,y=y)) +
    geom_line(data=curves.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query', group=ucov),lwd=1) +
    geom_line(data=curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference', group=ucov),lwd=1) +
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query', group=ucov),pch=19) +
    geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference', group=ucov),pch=19) +
    scale_color_manual(name='Data type',
                       breaks=c('Query',
                                'Reference'),
                       values=c('Query'='blue',
                                'Reference'='red'))

  gg <- gg + coord_cartesian(xlim=c(xstart,xstop))

  gg <- gg +
    labs(title=paste0(vachette_data$model.name,"; Observations + typical curves"),
         caption = paste0("Reference Covariate: ",
                          paste0(
                            names(vachette_data$covariates),
                            "=",
                            vachette_data$covariates,
                            collapse = ", "
                          )))

  if (log.x) {
    gg <- gg +
      labs(x = "ln(x)")
  }

  gg <- gg +
    render

  return(gg)
}


#' Plot Observations and Typical Curves Faceted by Covariate
#'
#' This function generates a ggplot2 visualization of observations and typical
#' curves within a pharmacometric model, faceted by unique covariate values (ucov).
#' The plot compares the observed data against the typical reference and query curves,
#' with each facet representing a different covariate value.
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing the observations and typical curves,
#' faceted by covariate. The plot displays both the reference and query curves
#' and their corresponding observations, color-coded for easy comparison.
#'
#' @details The function begins by ensuring that the input \code{vachette_data}
#' is of the correct class. It then constructs a plot that displays the reference
#' and query data, differentiating them using distinct colors, and faceting the
#' plot by unique covariate values. Observations are represented by points,
#' while the typical curves are depicted as lines.
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used. The x-axis range is dynamically set
#' based on the minimum and maximum x values, before and after scaling.
#' If the x-axis is logarithmic, the axis label is adjusted accordingly.
#'
#' @export
p.obs.cov <- function(vachette_data) {
  x <- y <- ref <- NULL
  stopifnot(inherits(vachette_data, "vachette_data"))
  log.x <- vachette_data$log.x

  obs.all <- vachette_data$obs.all
  curves.all <- vachette_data$curves.all

  xstart <- min(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x, vachette_data$obs.all$x.scaled)

  gg <- obs.all %>%
    ggplot(aes(x=x,y=y)) +
    geom_line(data=curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),lwd=1) +
    geom_line(data=curves.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query'),lwd=1) +
    geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),pch=19) +
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query'),pch=19) +

    scale_color_manual(name='Data type',
                       breaks=c('Query',
                                'Reference'),
                       values=c('Query'='blue',
                                'Reference'='red'))+
    facet_wrap(~paste("ucov",ucov))

  gg <- gg + coord_cartesian(xlim=c(xstart,xstop))

  gg <- gg +
    labs(title=paste0(vachette_data$model.name,"; Observations + typical curves"),
         caption = paste0("Reference Covariate: ",
                          paste0(
                            names(vachette_data$covariates),
                            "=",
                            vachette_data$covariates,
                            collapse = ", "
                          )))

  if (log.x) {
    gg <- gg +
      labs(x = "ln(x)")
  }

  gg <- gg +
    render

  return(gg)
}


#' Plot Excluded Observations with Reasons for Exclusion
#'
#' This function generates a ggplot2 visualization of observations that were
#' excluded from the transformation process within a pharmacometric model.
#' The plot includes both the included and excluded observations, with the
#' excluded ones color-coded by reason for exclusion.
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing the excluded observations and
#' their corresponding reasons for exclusion. The plot also displays the typical
#' reference and query curves for comparison.
#'
#' @details The function  plots the reference and query curves along
#' with the observations. Observations that were excluded from the transformation
#' process are highlighted and labeled with the reason for exclusion, while
#' included observations are shown in light grey.
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used. The x-axis range is dynamically set
#' based on the minimum and maximum x values from the original observations.
#' If the x-axis is logarithmic, the axis label is adjusted accordingly.
#'
#' @export
p.obs.excluded <- function(vachette_data) {

  x <- y <- ref <- reason <- NULL
  stopifnot(inherits(vachette_data, "vachette_data"))
  log.x <- vachette_data$log.x

  obs.all      <- vachette_data$obs.all
  obs.excluded <- vachette_data$obs.excluded
  curves.all   <- vachette_data$curves.all

  xstart <- min(vachette_data$obs.orig$x)
  xstop  <- max(vachette_data$obs.orig$x)

  # # For better legend
  obs.excluded$reason[obs.excluded$reason=="Less or equal to zero"] <- "Less or equal\nto zero"
  obs.excluded$reason[obs.excluded$reason=="Missing corresponding ref segment"] <- "Missing corresp.\nref segment"
  obs.excluded$reason[obs.excluded$reason=="Missing typical curve"] <- "Missing\ntypical curve"

  shapes <- c("Less or equal\nto zero" = 1,
              "Missing corresp.\nref segment" = 3,
              "Missing\ntypical curve" = 4)

  gg <- obs.excluded %>%
    ggplot(aes(x=x,y=y)) +
    geom_line(data=curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),lwd=1) +
    geom_line(data=curves.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query'),lwd=1) +
    scale_color_manual(
      name = 'Data type',
      breaks = c('Query',
                 'Reference'),
      values = c(
        'Query' = 'blue',
        'Reference' = 'red'
      )
    ) +
    geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y),pch=19,col='lightgrey') +
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y),pch=19,col='lightgrey') +
    geom_point(data=obs.excluded %>% filter(ref=='Yes'),aes(x=x,y=y,pch=factor(reason),col='Reference')) +
    geom_point(data=obs.excluded %>% filter(ref=='No'),aes(x=x,y=y,pch=factor(reason),col='Query')) +
    scale_shape_manual(values = shapes)

  gg <- gg + coord_cartesian(xlim=c(xstart,xstop))

  gg <- gg +
    labs(title=paste0(vachette_data$model.name,";\nObservations excluded from transformation (colored)"),
         caption = paste0("Reference Covariate: ",
                          paste0(
                            names(vachette_data$covariates),
                            "=",
                            vachette_data$covariates,
                            collapse = ", "
                          )),
         shape = "Reason\nexclusion:")

  if (log.x) {
    gg <- gg +
      labs(x = "ln(x)")
  }

  gg <- gg +
    render

  return(gg)
}

#' Plot Observations and \code{vachette} Transformed Observations with Connecting Arrows
#'
#' This function generates a ggplot2 visualization of observations and their transformations
#' within a pharmacometric model. The plot superimposes the observed query and reference data
#' from original to \code{vachette} transformed values, showing direction of transformation.
#' Reference and query curves are additionally overlaid for comparison.
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing the observations, their transformations,
#' and the corresponding reference and query curves. The plot also indicates if any
#' extrapolated reference curves are present and how many observations were excluded
#' from the transformation.
#'
#' @details The function plots the original and \code{vachette} transformed observations,
#' with arrows showing the transformation path. Reference and query curves are
#' also plotted, with extrapolated reference curves displayed as dashed lines
#' if available. The subtitle provides information about the error model used
#' (additive or proportional) and the number of excluded observations.
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used. The x-axis range is dynamically set
#' based on the scaled x values. If the x-axis is logarithmic, the axis label
#' is adjusted accordingly.
#'
#' @export
p.vachette.arrow <- function(vachette_data) {

  x <- y <- x.scaled <- y.scaled <- ref <- ucov <- seg <- NULL
  stopifnot(inherits(vachette_data, "vachette_data"))
  log.x <- vachette_data$log.x
  # JL 230607
  ref.extensions.all <- vachette_data$ref.extensions.all
  # Plot longest ref.extension.all only. Pick first is multiple occurences
  if(!is.null(ref.extensions.all)) max.x.ucov <- ref.extensions.all$ucov[ref.extensions.all$x == max(ref.extensions.all$x)][1]
  if(is.null(ref.extensions.all))  extensiontxt <- ""
  if(!is.null(ref.extensions.all)) extensiontxt <- "Dashed: extrapolation reference curve"

  n.excluded <- nrow(vachette_data$obs.excluded)

  # After scaling
  xstart <- min(vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x.scaled)

  gg <- vachette_data$obs.all %>%
    ggplot(aes(x=x,y=y)) +
    # Transformation arrows
    geom_segment(aes(x=x,y=y,xend=x.scaled,yend=y.scaled),
                 arrow = arrow(length = unit(0.2, "cm")),col='grey') +

    geom_line(data=vachette_data$curves.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query',group=ucov),lwd=1) +
    geom_point(data=vachette_data$obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query',group=ucov),pch=19,alpha=0.25) +
    geom_point(data=vachette_data$obs.all %>% filter(ref=='No'),aes(x=x.scaled,y=y.scaled,col='Transformed'),pch=19) +
    # Ref in top layer:
    geom_line(data=vachette_data$curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference',group=ucov),lwd=1)

  if(!is.null(ref.extensions.all))
  {
    gg <- gg + geom_line(
      data = ref.extensions.all %>% filter(ucov==max.x.ucov),
      aes(x = x, y = y, col = 'Reference',group=seg), lty=2,
      lwd = 0.8, col='grey30'
    )
  }

  gg <- gg + coord_cartesian(xlim=c(xstart,xstop))

  gg <- gg +
    geom_point(data=vachette_data$obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference',group=ucov),pch=19,alpha=0.25) +

    scale_color_manual(name='Data type',
                       breaks=c('Query',
                                'Reference',
                                'Transformed'),
                       values=c('Query'='blue',
                                'Reference'='red',
                                'Transformed' = 'purple'))


  gg <- gg + coord_cartesian(xlim=c(xstart,xstop))

  gg <- gg +
    labs(
      title = paste0(vachette_data$model.name, "; Observations + transformations"),
      subtitle = paste0(if (vachette_data$ADD_TR)
        "Additive Error; ", if (vachette_data$PROP_TR)
          "Proportional Error; ",extensiontxt,
        if(n.excluded>0) "\n",
        if(n.excluded>0) n.excluded,
        if(n.excluded>0) " observations not transformed and not displayed"),
      caption = paste0("Reference Covariate: ",
                       paste0(
                         names(vachette_data$covariates),
                         "=",
                         vachette_data$covariates,
                         collapse = ", "
                       ))
    )

  if (log.x) {
    gg <- gg +
      labs(x = "ln(x)")
  }

  gg <- gg +
    render

  return(gg)

}

#' Main Vachette Plot
#'
#' This function generates the main Vachette plot, which visualizes all transformed
#' data points within a pharmacometric model. The plot makes a clear distinction
#' between reference data points (those that were not moved, shown in red) and
#' transformed query data points (those that were moved, shown in purple).
#'
#' @inheritParams p.prop.distances
#'
#' @return A \code{ggplot2} object representing, which
#' includes all observations, with reference data points shown in red, and
#' transformed query data points shown in purple. The plot also includes
#' reference curves and, if applicable, extrapolated reference curves as dashed lines.
#'
#' @details The function plots all the data points, highlighting
#' reference points in red and transformed query points in purple. The plot
#' also overlays reference curves, with extrapolated segments displayed as
#' dashed lines if available. The subtitle provides information about the
#' error model used (additive or proportional) and the number of excluded
#' observations.
#'
#' The plot's title includes the model name, and the caption provides details
#' about the reference covariate(s) used. The x-axis range is dynamically set
#' based on the scaled x values. If the x-axis is logarithmic, the axis label
#' is adjusted accordingly.
#'
#' @export
p.vachette <- function(vachette_data) {
  x <- y <- x.scaled <- y.scaled <- ref <- ucov <- seg <- NULL

  stopifnot(inherits(vachette_data, "vachette_data"))
  log.x <- vachette_data$log.x

  obs.all            <- vachette_data$obs.all
  curves.all         <- vachette_data$curves.all

  n.excluded <- nrow(vachette_data$obs.excluded)

  # After scaling
  xstart <- min(vachette_data$obs.all$x.scaled)
  xstop  <- max(vachette_data$obs.all$x.scaled)

  # extract errors
  # errors <- vachette_data$errors
  # warning(errors)
  # JL 230607
  ref.extensions.all <- vachette_data$ref.extensions.all
  # Plot longest ref.extension.all only. Pick first is multiple occurences
  if(!is.null(ref.extensions.all)) max.x.ucov <- ref.extensions.all$ucov[ref.extensions.all$x == max(ref.extensions.all$x)][1]
  if(is.null(ref.extensions.all))  extensiontxt <- ""
  if(!is.null(ref.extensions.all)) extensiontxt <- "Dashed: extrapolation reference curve"

  gg <- obs.all %>%
    ggplot(aes(x = x, y = y)) +

    geom_point(
      data = obs.all %>% filter(ref == 'Yes'),
      aes(x = x, y = y, col = 'Reference'),
      pch = 19
    ) +
    geom_point(
      data = obs.all %>% filter(ref == 'No'),
      aes(x = x.scaled, y = y.scaled, col = 'Query\nTransformed'),
      pch = 19
    ) +

    # Line on top
    geom_line(
      data = curves.all %>% filter(ref == 'Yes'),
      aes(x = x, y = y, col = 'Reference'),
      lwd = 1
    )
  if(!is.null(ref.extensions.all))
  {
    gg <- gg + geom_line(
      data = ref.extensions.all %>% filter(ucov==max.x.ucov),
      aes(x = x, y = y, col = 'Reference', group=seg), lty=2,
      lwd = 0.8, col='grey30'
    )
  }
  gg <- gg +
    scale_color_manual(
      name = 'Data type',
      breaks = c('Query',
                 'Reference',
                 'Query\nTransformed'),
      values = c(
        'Query' = 'blue',
        'Reference' = 'red',
        'Query\nTransformed' = 'purple'
      )
    )

  gg <- gg + coord_cartesian(xlim=c(xstart,xstop))

  gg <- gg +
    labs(
      title = paste0(vachette_data$model.name, "; Observations + transformations"),
      subtitle = paste0(if (vachette_data$ADD_TR)
        "Additive Error; ", if (vachette_data$PROP_TR)
          "Proportional Error; ",extensiontxt,
        if(n.excluded>0) "\n",
        if(n.excluded>0) n.excluded,
        if(n.excluded>0) " observations not transformed and not displayed"),
      caption = paste0("Reference Covariate: ",
                       paste0(
                         names(vachette_data$covariates),
                         "=",
                         vachette_data$covariates,
                         collapse = ", "
                       ))
    )

  if (log.x) {
    gg <- gg +
      labs(x = "ln(x)")
  }

  gg <- gg +
    render

  return(gg)

}
