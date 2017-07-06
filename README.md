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

And then use it to install *crossr* by typing:

```r
devtools::install_github("othomantegazza/crossr", build_vignettes = TRUE)
```

By setting `build_vignettes = TRUE` (the default is `FALSE`), the installation
will take a little longer but in the process you will build the vignette that
explain how to use *crossr*.

Some images in the Vignette are built with the *pheatmap* package, which should be installed in R. If you haven't already, you can install the pheatmap package by typing:

```r
install.packages("pheatmap")
```

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
