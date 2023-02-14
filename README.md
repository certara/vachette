# vachette
A method for visualization of PMx models

Vachette builds off previous work done in [V2ACHER](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12679), providing a method to visualize PKPD analyses which are impacted by covariate effects.


## Examples

```r
oral_absorption_ex <- system.file(package = "vachette", "examples", "oral_absorption.R")
# Run script
# Note: "plots/vachette-main-oral-absorption-weight<version>.pdf" created in working directory
source(oral_absorption_ex)

# Alternatively, open file in RStudio for interactive execution of example
file.edit(oral_absorption_ex)
```

## Releases

The following releases prior to those archived on GitHub are described below:

* `v0.1.0` Based on vachette-functions-oral-absorption-v2.R
* `v0.2.0` indirect response example with updated functions
* `v0.3.0` same as v2 (align version with vachette-functions-v3.R)
* `v0.4.0` change choice of xlast ref, add i.v. example (no landmark)
* `v0.5.0` New folder structure, Simple Emax model
* `v0.6.0` Using splines for landmark finding
* `v0.7.0` Using Savitzky-Golay for landmark finding
* `v0.8.0` Two-dose
* `v0.9.0` Inspection 2nd order derivatives for inflection points
* `v0.10.0/v.0.11.0` 2-step landmark finding (approx followed by refined)
* `v0.12.0` 
  - Tolerance noise dependent on grid size "tolnoise"
  - Tolerance last data point (how far off asymptote) by user "tolend"
  - Really pass tolerances to functions
* `v0.13.0` Test several models
* `v0.14.0` Repair single dose "reverse", Repair Emax, Repair indirect response
* `v0.15.0` test use of simulation tolerances indirect response
* `v0.16.0`  Pembro model

## Glossary

* Region: Region is the Vachette terminology for the time between two dose administrations
* Region Type:  May be open or closed
* Covariate Combinations Number: Number of unique covariate combinations
