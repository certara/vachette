
# vachette

A method for visualization of PMx models.

Vachette builds off previous work done in
[V2ACHER](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12679),
providing a method to visualize PKPD analyses which are impacted by
covariate effects.

## Installation

**Installation instructions will be available once the Github repository
has been made public**

## Examples

### IV

``` r
# Import Example Data
obs.data <- read.csv(system.file(package = "vachette", "examples", "iv-obs.csv")) |> 
  dplyr::rename(WT = "vachette.cov1")
typ.data <- read.csv(system.file(package = "vachette", "examples", "iv-typ.csv")) |> 
  dplyr::rename(WT = "vachette.cov1")

library(vachette)

vd <-
  vachette_data(
    obs.data,
    typ.data,
    covariates = c("WT" = 70),
    ref.dosenr = 1,
    iiv.correction = FALSE,
    model.name = "Intravenous (i.v)"
  )
```

    ## `dosenr` column found in obs.data, using `dosenr` column in data for corresponding ref.dosenr value

    ## `dosenr` column found in typ.data, using `dosenr` column in data for corresponding ref.dosenr value

``` r
vd <- vd |>
  apply_transformations()

p.vachette(vd)
```

![](README_files/figure-gfm/iv-1.png)<!-- -->

### Oral Two-Cov

``` r
obs.data <- read.csv(system.file(package = "vachette", "examples", "oral-two-cov-obs.csv"))
typ.data  <- read.csv(system.file(package = "vachette", "examples", "oral-two-cov-typ.csv"))

vd <-
  vachette_data(
    obs.data,
    typ.data,
    covariates =  c("WT" = 70, "AGE" = 30),
    ref.dosenr = 1,
    iiv.correction = FALSE,
    model.name = "oral-two-cov"
  )
```

    ## `dosenr` column found in obs.data, using `dosenr` column in data for corresponding ref.dosenr value

    ## `dosenr` column found in typ.data, using `dosenr` column in data for corresponding ref.dosenr value

``` r
vd <- vd |>
  apply_transformations(window = 23,
                        window.d1.refine = 7,
                        window.d2.refine = 5)

p.vachette(vd)
```

![](README_files/figure-gfm/oral_two_cov-1.png)<!-- -->

### Sigmoid

``` r
obs.data <- read.csv(system.file(package = "vachette", "examples", "sigmoid-obs.csv")) |>
  dplyr::rename("ID" = id)
typ.data  <- read.csv(system.file(package = "vachette", "examples", "sigmoid-typ.csv"))

vd <-
  vachette_data(
    obs.data,
    typ.data,
    covariates =  c("vachette.cov1" = 70),
    ref.dosenr = 1,
    iiv.correction = FALSE,
    model.name = "sigmoid"
  )
```

    ## `dosenr` column found in obs.data, using `dosenr` column in data for corresponding ref.dosenr value

    ## `dosenr` column found in typ.data, using `dosenr` column in data for corresponding ref.dosenr value

``` r
vd <- vd |>
  apply_transformations(window = 17,
                        window.d1.refine = 7,
                        window.d2.refine = 5)

p.vachette(vd)
```

![](README_files/figure-gfm/sigmoid-1.png)<!-- -->

### Oral-Absorption

``` r
obs.data <- read.csv(system.file(package = "vachette", "examples", "oral-absorption-obs.csv"))
typ.data  <- read.csv(system.file(package = "vachette", "examples", "oral-absorption-typ.csv"))

vd <-
  vachette_data(
    obs.data,
    typ.data,
    covariates =  c("vachette.cov1" = 70),
    ref.dosenr = 1,
    iiv.correction = FALSE,
    model.name = "oral-absorption"
  )
```

    ## `dosenr` column found in obs.data, using `dosenr` column in data for corresponding ref.dosenr value

    ## `dosenr` column found in typ.data, using `dosenr` column in data for corresponding ref.dosenr value

``` r
vd <- vd |>
  apply_transformations(window = 17,
                        window.d1.refine = 7,
                        window.d2.refine = 5)

p.vachette(vd)
```

![](README_files/figure-gfm/oral_absorption-1.png)<!-- -->

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
    iiv.correction = FALSE,
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
    covariates =  c("SCHED", "ALB"),
    ref.dosenr = 3,
    iiv.correction = FALSE,
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

p.vachette(vd)
```

![](README_files/figure-gfm/pembro3-1.png)<!-- -->
