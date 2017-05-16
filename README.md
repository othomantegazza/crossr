## Crossr

### Overview

*Crossr* is set of functions and a workflow that can be useful for quantitative comparison of transcriptomic data from different species.

*Crossr* is under development and comes with **no warranty**.

### Install

in order to install *crossr* from Github you'll need to have the *devtools*
package installed in R.

To install *devtools* from CRAN and load it, type:

```r
install.packages("devtools")
```
And then use it to install *crossr* with the function

```r
devtools::install_github("othomantegazza/crossr", build_vignettes = TRUE)
```
By setting `build_vignettes = TRUE` (the default is `FALSE`), the installation
takes a little longer but in the process you will build the vignettes that
explain how to use *crossr*.

### Use crossr

The workflow and usage of *crossr* is documented in the vignette, that can be
accessed with:

```r
vignette("quick_overview", "crossr")
```
