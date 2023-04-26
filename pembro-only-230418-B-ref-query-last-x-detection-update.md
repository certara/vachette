
# vachette

A method for visualization of PMx models.

Vachette builds off previous work done in
[V2ACHER](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12679),
providing a method to visualize PKPD analyses which are impacted by
covariate effects.

## Installation

**Installation instructions will be available once the Github repository
has been made public**

``` r
library(vachette)
```

## Examples

### Pembro

**Note: The Pembro data set contains a categorical covariate e.g.,
`SCHED`**

``` r
obs.data <- read.csv(system.file(package = "vachette", "examples", "pembro-obs.csv"))
typ.data  <- read.csv(system.file(package = "vachette", "examples", "pembro-typ.csv"))

vd <-
  vachette_data(
    obs.data,
    typ.data,
    covariates =  c(SCHED = 'Q3W', ALB = 53.5),
    ref.dosenr = 3,
    model.name = "pembro"
  )

print(vd)
```

    ## Model Name:      pembro 
    ## Covariate Names: SCHED, ALB 
    ## Reference Values:    SCHED=Q3W , ALB=53.5

If no reference value is specified, the covariate central tendency is
used e.g., median for continuous and mode for categorical covariate.

``` r
vd <-
  vachette_data(
    obs.data,
    typ.data,
    covariates =  c("SCHED" = "Q2W", "ALB"= "16"),
    ref.dosenr = 1,
    model.name = "pembro"
  )

print(vd)
```

    ## Model Name:      pembro 
    ## Covariate Names: SCHED, ALB 
    ## Reference Values:    SCHED=Q2W , ALB=16

``` r
vd <- vd |>
  apply_transformations(window = 17,
                        window.d1.refine = 7,
                        window.d2.refine = 5)
```

    ## [1] "**** EXTEND REFERENCE CURVE FOR query i.ucov 7 (Region: 1) *****"
    ## [1] "**** END EXTENSION *****"
    ## [1] "**** EXTEND REFERENCE CURVE FOR query i.ucov 8 (Region: 2) *****"
    ## [1] "**** END EXTENSION *****"
    ## [1] "**** EXTEND REFERENCE CURVE FOR query i.ucov 9 (Region: 1) *****"
    ## [1] "**** END EXTENSION *****"
    ## [1] "**** EXTEND REFERENCE CURVE FOR query i.ucov 10 (Region: 2) *****"
    ## [1] "**** END EXTENSION *****"

``` r
p.vachette(vd)
```

![](pembro-only-230418-B-ref-query-last-x-detection-update_files/figure-gfm/pembro3-1.png)<!-- -->
