# Modeling and leveraging intuitive theories to improve vaccine attitudes

Authors: Derek Powell, Kara Weisman, & Ellen Markman

Manuscript under review. Code and text may change in light of reviewer comments. Please contact if you have any questions about this project or code. For related work, see our [CogSci paper](https://mindmodeling.org/cogsci2018/papers/0183/0183.pdf) and [this repo](https://github.com/derekpowell/vaccbeliefs-cogsci2018).

## Abstract

> Much of the richness of human thought is supported by people’s intuitive theories---mental frameworks capturing the perceived structure of the world. But intuitive theories can also contain and reinforce dangerous misconceptions. In this paper, we take up the case of misconceptions about vaccine safety that discourage vaccination. These misconceptions constitute a major public health risk that predates the coronavirus pandemic but that has become all the more dire in recent years. We argue that addressing such misconceptions requires awareness of the broader conceptual contexts in which they are embedded. To build this understanding, we examined the structure and revision of people’s intuitive theories of vaccination in five large survey studies (total _N_ = 3196). Based on these data, we present a cognitive model of the intuitive theory surrounding people’s decisions about whether to vaccinate young children against diseases like measles, mumps, and rubella (MMR). Using this model, we were able to make accurate predictions about how people’s beliefs would be revised in light of educational interventions, design an effective new intervention encouraging vaccination, and understand how these beliefs were affected by real-world events (the measles outbreaks of 2019). In addition to presenting a promising way forward for promoting the MMR vaccine, this approach has clear implications for encouraging the uptake of COVID-19 vaccines, especially among parents of young children. At the same time, this work provides the foundation for richer understandings of intuitive theories and belief revision more broadly.

## Repository Overview

* `paper/`: Rmarkdown and supporting files for creation of the manuscript. Must run supplemental materials notebook first.
* `supplement/`: Rmarkdown and supporting files for creation of the supplemental materials.
* `data/`: Data for all studies reported in paper.
* `code/`: Code scripts on which all Rmarkdown notebooks depend.
  * `../custom-structure-learning/`: Code implementing custom scoring algorithms and plotting tools used in Bayesian Network structure learning and inference.
* `outputs/`: stores outputs from running `supplement-main.Rmd` which are read to create paper.
* `pilots/`: Anonymized data and analysis notebooks for pilot studies. These are being released for transparency but have not been commented or edited carefully for public consumption.

## Reproducibility notes

To reproduce the manuscript and all analyses, follow the following steps after cloning this repository.

1. Create a `local/` folder in the repository (at the terminal: `mkdir local`)
2. **RECOMMENDED**: Install [Docker](https://www.docker.com/) and utilize `cogdatasci/rstudio` docker container. Recommend you run with the following options:
```bash
docker run -d -p 8787:8787 -v "`pwd`":/home/rstudio/working \
 -e PASSWORD=my_password_here cogdatasci/rstudio
 ```
then navigate to `localhost:8787` in your browser to access the rstudio interface.

3. Install required packages by running `install-packages.R` to make sure you have all packages needed. You will also need to install `jags`. In the Docker container this can be done with `sudo apt-get install jags` at terminal.
4. Open and knit `supplement/supplement-main.Rmd` to generate supplement PDF and save files needed for reproduction of manuscript.
5. Open and knit `paper/paper-main.Rmd` to generate manuscript PDF.
