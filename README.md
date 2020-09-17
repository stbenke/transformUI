# transformUI

A simple Shiny gadget that provides a graphical user interface to some of the transformation functions in the [flowCore package](https://github.com/RGLab/flowCore) for flow cytometry data.

The following transformations are currently supported:

- Log10
- logicle
- asinh

Automatic parameter estimation is implemented for logicle using functions from the flowCore package and for asinh from the [flowVS package](http://bioconductor.org/packages/release/bioc/html/flowVS.html).

## Installation

First install flowCore and flowVS from Bioconductor, then transformUI:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("flowCore", "flowVS"))

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("stbenke/transformUI")
```

## Licence
This R package is licensed under the GPL-3 license (see LICENSE file).
