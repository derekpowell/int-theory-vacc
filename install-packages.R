## this will install all required packages
## untested 6/22/21, 11:03 AM

devtools::install_bioc("graph")

package_list <- c(
  "tidyverse",
  "brms",
  "bnlearn",
  "HydeNet",
  "betareg",
  "kableExtra",
  "DiagrammeR",
  "DiagrammeRsvg",
  "rsvg",
  "corrplot",
  "cowplot",
  "caret"
)

install.packages(package_list)

devtools::install_github("thomasp85/patchwork", "crsh/papaja")

tinytex::install_tinytex()
