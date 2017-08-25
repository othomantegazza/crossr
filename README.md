# Crossr

[![Build Status](https://travis-ci.org/othomantegazza/crossr.svg?branch=master)](https://travis-ci.org/othomantegazza/crossr)

### Overview

*Crossr* is set of functions and a workflow that can be useful for quantitative comparison of transcriptomic data from different species.

*Crossr* is under development and comes with **no warranty**.

### Install

In order to install *crossr* from Github you must first install the *devtools*
package in R.

To install *devtools* from CRAN, type:

```r
install.packages("devtools")
```
Some images in *crossr*'s Vignette require the *ggplot2* and *pheatmap* packages. To build the Vignette, 
those packages must be installed in R.
*ggplot2* is also needed for other plotting functionalities of *crossr*.

If you haven't installed *ggplot2* and *pheatmap* already, you can install them by typing:

```r
install.packages("ggplot2", "pheatmap")
```

Then you can use *devtools* to install *crossr* by typing:

```r
devtools::install_github("othomantegazza/crossr", build_vignettes = TRUE)
```

By setting `build_vignettes = TRUE` (the default is `FALSE`), the installation
will take a little longer but in the process you will build the Vignette that
explains how to use *crossr*.



### Use crossr

First load *crossr* with:

```r
library(crossr)
```

The workflow and usage of *crossr* is documented in the vignette, that can be
accessed with:

```r
vignette("quick_overview", "crossr")
```
