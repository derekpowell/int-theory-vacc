## this will install all required packages
## untested 6/22/21, 11:03 AM

package_list <- c(
  "tidyverse",
  "brms",
  "bnlearn",
  "HydeNet",
  "betareg",
  "kableExtra",
  "DiagrammeR",
  "DiagrammeRsvg",
  "corrplot",
  "cowplot",
  "caret"
)

install.packages(package_list)

devtools::install_github("thomasp85/patchwork")