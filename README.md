# Modeling intuitive theories surrounding vaccination decisions

Authors: Derek Powell, Kara Weisman, & Ellen Markman

Repository related to psychological studies exploring methods for capturing intuitive theories as Bayesian networks using structure learning techniques.

Manuscript in preparation, title may change. Nothing here is final. please contact if you have any questions about this project or code. For related work, see our [CogSci paper](https://mindmodeling.org/cogsci2018/papers/0183/0183.pdf) and [this repo](https://github.com/derekpowell/vaccbeliefs-cogsci2018).

## _Related abstract_

> How can we leverage the cognitive science of lay theories to inform interventions aimed at correcting misconceptions and changing behaviors? Focusing on the problem of vaccine skepticism, we identified a set of 14 beliefs we hypothesized would be relevant to vaccination decisions. We developed reliable scales to measure these beliefs across a large sample of participants (n = 1130) and employed state-of-the-art graphical structure learning algorithms to uncover the relationships among these beliefs. This resulted in a graphical model describing the system of beliefs relevant to childhood vaccinations, with beliefs represented as nodes and their interconnections as directed edges. This model sheds light on how these beliefs relate to one another and can be used to predict how interventions aimed at specific beliefs will play out across the larger system. Moving forward, we hope this modeling approach will help guide the development of effective, theory-based interventions promoting childhood vaccination.

## Repository Overview

* `paper/`: Rmarkdown and supporting files for creation of the manuscript. Must run supplemental materials notebook first.
* `supplement/`: Rmarkdown and supporting files for creation of the supplemental materials.
* `data/`: Data for all studies reported in paper.
* `code/`: Code scripts on which all Rmarkdown notebooks depend.
  * `../custom-structure-learning/`: Code implementing custom scoring algorithms and plotting tools used in Bayesian Network structure learning and inference.
* `pilots/`: Anonymized data and analysis notebooks for pilot studies. These are being released for transparency but have not been commented or edited carefully for public consumption.

## Reproducibility notes

To reproduce the manuscript and all analyses, follow the following steps after cloning this repository.

1. Create a `local/` folder in the repository (at the terminal: `mkdir local`)
2. Install required packages with or install [Docker](https://www.docker.com/) and utilize `cogdatasci/rstudio` docker container. Recommend you run with the following options:
```bash
docker run -d -p 8787:8787 -v "`pwd`":/home/rstudio/working \
 -e PASSWORD=my_password_here cogdatasci/rstudio
 ```
then navigate to `localhost:8787` in your browser to access the rstudio interface. Either way you can run `install-packages.R` to make sure you have all packages needed.
 
3. Open and knit `supplement/supplement-main.Rmd` to generate supplement PDF and save files needed for reproduction of manuscript.
4. (In the near future) Open and knit `paper/paper-main.Rmd` to generate manuscript PDF.
