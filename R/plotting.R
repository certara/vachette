
#' Plot theme for vachette
#'
#' A ggplot2 theme for vachette plot output
#'
#' @return \code{theme, gg}
#' @export
#'
#' @examples
#' theme_vachette()
#'
theme_vachette <- function() {
  ggplot2::theme_bw()+
  ggplot2::theme(text = ggplot2::element_text(size=16),
        axis.text.x = ggplot2::element_text(size=14),
        axis.text.y = ggplot2::element_text(size=14),
        plot.title = ggplot2::element_text(size=14),
        plot.subtitle=ggplot2::element_text(size=12))

}


# For same looks and feel of plots:
render <- theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        plot.title = element_text(size=14),
        plot.subtitle=element_text(size=12))

# Multiple plots on one page:
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#' p.scaled.typical.curves.landmarks
#'
#' @param vachette_data Object of class\code{vachette_data}
#'
#' @export
p.scaled.typical.curves.landmarks <- function(vachette_data) {

  curves.all <- vachette_data$curves.all
  lm.all <- vachette_data$lm.all
  model.name <- vachette_data$model.name
  xstop <- vachette_data$xstop

  stopifnot(!is.null(curves.all))

  curves.all %>%
  ggplot(aes(x=x,y=y,group=ucov,
             col = factor(seg)))+
  geom_line(lwd=1)+
  geom_line(data=curves.all %>% filter(ref=="Yes"),col='black',lty=2,lwd=1,alpha=0.5)+
  # Add landmark positions
  geom_point(data=lm.all,pch=3,size=4,col='black',stroke = 2)+
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title=paste0(model.name,". Typical curve segments"),
       subtitle = "Dashed line: Reference typical curve",
      # caption=script,
       col="Segment")+
  render
}

#' p.scaled.typical.full.curves.landmarks
#'
#' @param vachette_data Object of class\code{vachette_data}
#'
#' @export
p.scaled.typical.full.curves.landmarks <- function(vachette_data) {
  curves.all <- vachette_data$curves.all
  lm.all <- vachette_data$lm.all
  model.name <- vachette_data$model.name

  curves.all %>%
  ggplot(aes(x=x,y=y,group=ucov,
             col = factor(seg)))+
  geom_line(lwd=1)+
  geom_line(data=curves.all %>% filter(ref=="Yes"),col='black',lty=2,lwd=1,alpha=0.5)+
  # Add landmark positions
  geom_point(data=lm.all,pch=3,size=4,col='black',stroke = 2)+
  labs(title=paste0(model.name,". Typical curve segments"),
       subtitle = "Dashed line: Reference typical curve",
      # caption=script,
       col="Segment")+
  render
}

#' p.scaling.factor
#'
#' @param vachette_data Object of class\code{vachette_data}
#'
#' @export
p.scaling.factor <- function(vachette_data) {
  curves.scaled.all <- vachette_data$curves.scaled.all
  obs.all           <- vachette_data$obs.all
  model.name        <- vachette_data$model.name

 curves.scaled.all %>%
  mutate(mycurve = ifelse(ref=='Yes','Reference','Query')) %>%
  ggplot(aes(x=x,y=x.scaling,col=factor(seg)))+
  geom_line(lwd=1)+
  facet_wrap(~paste(mycurve," covariate =",COV))+
  coord_cartesian(xlim=c(NA,max(obs.all$x,obs.all$x.scaled)),
                  ylim=c(0,max(curves.scaled.all$x.scaling[curves.scaled.all$x<=max(obs.all$x,obs.all$x.scaled)])))+
  labs(title=paste0(model.name,". x-scaling factors"),
       subtitle = paste0(if(vachette_data$ADD_TR) "Additive Error",if(vachette_data$PROP_TR) "Proportional Error"),
       #caption=script,
       col="Segment") +
  render
}

#' p.scaled.typical.curves
#'
#' @param vachette_data Object of class \code{vachette_data}
#'
#' @return Object of class \code{ggplot2}
#' @export
#'
p.scaled.typical.curves <- function(vachette_data) {
  curves.scaled.all <- vachette_data$curves.scaled.all
  model.name        <- vachette_data$model.name

  curves.scaled.all %>%
  ggplot(aes(x=x,y=y,group=ucov))+
  # geom_line(lwd=1.5,alpha=0.5)+
  geom_line(data=curves.scaled.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Yes'),lwd=1.5,alpha=0.60)+
  geom_line(data=curves.scaled.all %>% filter(ref=='No'),aes(x=x,y=y,col='No'),lwd=1.5)+
  geom_line(aes(x=x.scaled, y=y.scaled), col='black',lwd=0.75,lty=2)+
  # scale_color_manual(values=c('red','blue')) +
  # scale_linetype_manual(values = c(3, 1),
  #                       labels = c("Covariate", "Reference")) +
  coord_cartesian(xlim=c(0,vachette_data$xstop)) +
  scale_color_manual(name='Reference\ncurve',
                     breaks=c('No', 'Yes'),
                     values=c('No'='blue',
                              'Yes'='red'))+
  labs(title=paste0(model.name,". Covariate typical curves"),
       subtitle = "Dashed lines: Covariate typical curves after Vachette transformation",
      # caption=script,
       col="Covariate value\n(Reference)")+
  render
}

#' p.scaled.observation.curves
#'
#' @param vachette_data Object of class \code{vachette_data}
#'
#' @return Object of class \code{ggplot2}
#' @export
#'
p.scaled.observation.curves <- function(vachette_data) {
# Copy ref curves to each ID
obs.all <- vachette_data$obs.all
myids  <- unique(obs.all$ID)
curves.all <- vachette_data$curves.all
curves.all.ids <- NULL
for(iid in c(1:length(myids))) {curves.all.ids <- rbind(curves.all.ids,curves.all %>% mutate(ID=myids[iid]))}

# Overlay observation curves per ID
  obs.all %>%
  ggplot(aes(x=x,y=y,group=paste(ID,COV)))+
  geom_line(data=curves.all.ids %>% filter(ref=="Yes"),aes(x=x,y=y,col='Typical reference'),lwd=1)+
  geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
  geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
  geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'),lty=2)+
  geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'))+
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
  coord_cartesian(xlim=c(0,vachette_data$xstop)) +
  labs(title=paste0(vachette_data$model.name,". Individual observation curves"),
       subtitle = "Dashed lines: Individual observation curves after Vachette transformation",
      # caption=script,
       col="Covariate value\n(Reference)")+
  render
}

#' p.scaled.observation.curves.by.id
#'
#' @param vachette_data Object of class \code{vachette_data}
#'
#' @return Object of class \code{ggplot2}
#' @export
#'
p.scaled.observation.curves.by.id <- function(vachette_data) {
# Observation curve for each ID
  obs.all <- vachette_data$obs.all
  myids  <- unique(obs.all$ID)
  curves.all <- vachette_data$curves.all
  curves.all.ids <- NULL
for(iid in c(1:length(myids))) {curves.all.ids <- rbind(curves.all.ids,curves.all %>% mutate(ID=myids[iid]))}
  obs.all %>%
  ggplot(aes(x=x,y=y,group=paste(ID,COV)))+
  geom_line(data=curves.all.ids %>% filter(ref=="Yes"),aes(x=x,y=y,col='Typical reference'),lwd=1)+
  geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
  geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
  geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'),lty=2)+
  geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'))+
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
  facet_wrap(~ID)+
  coord_cartesian(xlim=c(0,vachette_data$xstop)) +
  labs(title=paste0(vachette_data$model.name,". Individual observation curves"),
       subtitle = "Dashed lines: Individual observation curves after Vachette transformation",
       #caption=script,
       col="Covariate value\n(Reference)")+
  render
}


#' p.transformation
#'
#' @param vachette_data Object of class \code{vachette_data}
#'
#' @return Object of class \code{ggplot2}
#' @export
#'
p.transformation <- function(vachette_data) {
# Select first or second segment as example and first region
  nseg <- vachette_data$nseg
  curves.all <- vachette_data$curves.all
  obs.all <- vachette_data$obs.all
  p.seg    <- ifelse(nseg > 1, 2, 1)
  p.region <- 1
  obs.all %>%
    ggplot(aes(
      x = x,
      y = y,
      col = factor(seg),
      pch = factor(paste(!!!syms(names(vachette_data$covariates))))
    )) +
    # geom_point(data = obs.all %>% filter(seg==2),size=2,alpha=0.25) +
    # geom_point(data = obs.all %>% filter(ref=="No",seg==2),aes(x=x.scaled,y=y.scaled),col='purple',size=2) +
    # Transformation arrows
    # geom_segment(aes(x=x,y=y,xend=x.scaled,yend=y.scaled),
    #              arrow = arrow(length = unit(0.2, "cm")),col='grey') +
    geom_line(
      data = curves.all %>% filter(ref != "No"),
      lwd = 1,
      col = 'grey'
    ) +
    geom_line(
      data = curves.all %>% filter(ref == "No"),
      lwd = 1,
      col = 'grey'
    ) +
    geom_line(
      data = curves.all %>% filter(ref != "No", seg == p.seg),
      lwd = 1,
      col = 'red'
    ) +
    geom_line(
      data = curves.all %>% filter(ref == "No", seg == p.seg, region == p.region),
      lwd = 1,
      col = 'blue'
    ) +
    geom_line(
      data = curves.all %>% filter(ref == "No", seg == p.seg, region == p.region),
      aes(x = x.scaled),
      lwd = 1,
      col = 'cyan'
    ) +
    geom_line(
      data = curves.all %>% filter(ref == "No", seg == p.seg, region == p.region),
      aes(x = x.scaled, y = y.scaled),
      lwd = 1,
      col = 'purple',
      lty = 2
    ) +
    coord_cartesian(xlim = c(0, max(obs.all$x, obs.all$x.scaled))) +
    labs(
      title = paste0(
        vachette_data$model.name,
        " - Example Vachette segment-",
        p.seg,
        " transformation"
      ),
      subtitle = paste0(if (vachette_data$ADD_TR)
        "Additive Error", if (vachette_data$PROP_TR)
          "Proportional Error"),
     # caption = script,
      col = "Segm",
      shape = "Covs"
    ) +
    render
}

#' p.add.distances
#'
#' @param vachette_data Object of class \code{vachette_data}
#'
#' @return Object of class \code{ggplot2}
#' @export
#'
p.add.distances <- function(vachette_data) {
  obs.all <- vachette_data$obs.all
# Additive distances before and after transformation
obs.all %>%
  ggplot(aes(x=dist.add.orig,y=dist.add.transformed,col=factor(seg)))+
  geom_abline(slope=1)+
  geom_point()+
  labs(title=paste0(vachette_data$model.name,". Normal distances: original and after transformation"),
       subtitle = paste0(if(vachette_data$ADD_TR) "Additive Error",if(vachette_data$PROP_TR) "Proportional Error"),
       #caption=script,
       x = 'Original distance',
       x = 'Distance after transformation',
       col="Segm.")+
  render
}

#' p.prop.distances
#'
#' @param vachette_data Object of class \code{vachette_data}
#'
#' @return Object of class \code{ggplot2}
#' @export
#'
p.prop.distances <- function(vachette_data) {
  obs.all <- vachette_data$obs.all
# Proportional distances before and after transformation
obs.all %>%
  ggplot(aes(x=dist.prop.orig,y=dist.prop.transformed,col=factor(seg)))+
  geom_abline(slope=1)+
  geom_point()+
  labs(title=paste0(vachette_data$model.name,". Proportional distances: original and after transformation"),
       subtitle = paste0(if(vachette_data$ADD_TR) "Additive Error",if(vachette_data$PROP_TR) "Proportional Error"),
       #caption=script,
       x = 'Original distance on log scale',
       x = 'Distance on log scale after transformation',
       col="Segm.")+
  render
}

#' p.obs.ref.query
#'
#' @param vachette_data Object of class \code{vachette_data}
#'
#' @return Object of class \code{ggplot2}
#' @export
#'
p.obs.ref.query <- function(vachette_data) {
  #stopifnot(length(vachette_data$covariates) == 1)
  obs.all <- vachette_data$obs.all
  curves.all <- vachette_data$curves.all
  obs.all %>%
  ggplot(aes(x=x,y=y)) +
  geom_line(data=curves.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query', group=ucov),lwd=1) +
  geom_line(data=curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference', group=ucov),lwd=1) +
  geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query', group=ucov),pch=19) +
  geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference', group=ucov),pch=19) +

  scale_color_manual(name='Data type',
                     breaks=c('Query',
                              'Reference'),
                     values=c('Query'='blue',
                              'Reference'='red'))+
  # coord_cartesian(xlim=c(0,max(obs.all$x))) +
    coord_cartesian(xlim=c(0,max(vachette_data$lm.all$x))) +
  labs(title=paste0(vachette_data$model.name,". Observations + typical curves"),
       subtitle = paste0(if(vachette_data$ADD_TR) "Additive Error",if(vachette_data$PROP_TR) "Proportional Error")
       #caption=script
       )+
  render
}

#' p.obs.cov
#'
#' @param vachette_data Object of class \code{vachette_data}
#'
#' @return Object of class \code{ggplot2}
#' @export
#'
p.obs.cov <- function(vachette_data) {
  #stopifnot(length(vachette_data$covariates) == 1)
  obs.all <- vachette_data$obs.all
  curves.all <- vachette_data$curves.all

  obs.all %>%
  ggplot(aes(x=x,y=y)) +
  geom_line(data=curves.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query'),lwd=1) +
    geom_line(data=curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),lwd=1) +
  geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query'),pch=19) +
    geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),pch=19) +

  scale_color_manual(name='Data type',
                     breaks=c('Query',
                              'Reference'),
                     values=c('Query'='blue',
                              'Reference'='red'))+
  facet_wrap(~paste("ucov",ucov)) +
  coord_cartesian(xlim=c(0,max(obs.all$x))) +
  labs(title=paste0(vachette_data$model.name,". Observations + typical curves"),
       subtitle = paste0(if(vachette_data$ADD_TR) "Additive Error",if(vachette_data$PROP_TR) "Proportional Error")
       #caption=script
       )+
  render
}

#' p.vachette.arrow
#'
#' @param vachette_data Object of class \code{vachette_data}
#'
#' @return Object of class \code{ggplot2}
#' @export
#'
#'
p.vachette.arrow <- function(vachette_data) {

  vachette_data$obs.all %>%
  ggplot(aes(x=x,y=y)) +
  # Transformation arrows
  geom_segment(aes(x=x,y=y,xend=x.scaled,yend=y.scaled),
               arrow = arrow(length = unit(0.2, "cm")),col='grey') +

  geom_line(data=vachette_data$curves.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query',group=ucov),lwd=1) +
    geom_line(data=vachette_data$curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference',group=ucov),lwd=1) +
  geom_point(data=vachette_data$obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query',group=ucov),pch=19,alpha=0.25) +
  geom_point(data=vachette_data$obs.all %>% filter(ref=='No'),aes(x=x.scaled,y=y.scaled,col='Transformed'),pch=19) +
    geom_point(data=vachette_data$obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference',group=ucov),pch=19,alpha=0.25) +

  scale_color_manual(name='Data type',
                     breaks=c('Query',
                              'Reference',
                              'Transformed'),
                     values=c('Query'='blue',
                              'Reference'='red',
                              'Transformed' = 'purple'))+

  coord_cartesian(xlim=c(0,max(vachette_data$obs.all$x,vachette_data$obs.all$x.scaled))) +
  labs(title=paste0(vachette_data$model.name,". Observations + transformations"),
       subtitle = paste0(if(vachette_data$ADD_TR) "Additive Error",if(vachette_data$PROP_TR) "Proportional Error")
       #caption=script
       )+
  render
}

#' p.vachette
#'
#' @param vachette_data Object of class \code{vachette_data}
#'
#' @return Object of class \code{ggplot2}
#' @export
#'
p.vachette <- function(vachette_data) {
  #stopifnot(length(vachette_data$covariates) == 1)
  obs.all <- vachette_data$obs.all
  curves.all <- vachette_data$curves.all

  obs.all %>%
    ggplot(aes(x = x, y = y)) +
    geom_line(
      data = curves.all %>% filter(ref == 'Yes'),
      aes(x = x, y = y, col = 'Reference'),
      lwd = 1
    ) +
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
    ) +

    coord_cartesian(xlim = c(0, max(obs.all$x, obs.all$x.scaled))) +
    labs(
      title = paste0(vachette_data$model.name, ". Observations + transformations"),
      subtitle = paste0(if (vachette_data$ADD_TR)
        "Additive Error", if (vachette_data$PROP_TR)
          "Proportional Error")
      #caption=script
    ) +
    render
}
