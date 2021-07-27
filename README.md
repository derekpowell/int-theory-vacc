# Modeling and leveraging intuitive theories to improve vaccine attitudes

Authors: Derek Powell, Kara Weisman, & Ellen Markman

Repository related to psychological studies exploring methods for capturing intuitive theories as Bayesian networks using structure learning techniques.

Manuscript under review. Code and text may change in light of reviewer comments. Please contact if you have any questions about this project or code. For related work, see our [CogSci paper](https://mindmodeling.org/cogsci2018/papers/0183/0183.pdf) and [this repo](https://github.com/derekpowell/vaccbeliefs-cogsci2018).

## Abstract

> Much of the richness of human thought is supported by people’s intuitive theories---explanatory mental frameworks that capture the perceived causal and relational structure of the world. But serious consequences can arise when these theories are mistaken. In this paper, we consider an important and timely example of these risks: when misconceptions about vaccine safety discourage vaccination. We argue that addressing misconceptions requires awareness of the broader conceptual contexts in which they are embedded. To this end, we sought to develop a cognitive model of the intuitive theory surrounding vaccination decisions. From qualitative research, we identified a set of beliefs relevant to this intuitive theory. Then, we built our model by surveying a large sample of U.S. adults and applying Bayesian network structure learning algorithms to their responses to uncover the connections among these beliefs (Study 1). The resulting model supports understanding of belief revision within this larger conceptual system: Our model helped to explain belief revision following an established educational intervention (Study 2), provided insight that supported our formulation of a novel educational intervention that successfully increased vaccination intentions (Studies 3a, 3b), and also helped to explain how people’s beliefs changed following real-world events (Study 4). Taken together, this work illustrates the promise of understanding belief revision in terms of intuitive theories and of discovering and formalizing intuitive theories via computational models. Moreover, it lays the foundation for theory-based educational interventions that could encourage people to vaccinate their children and themselves against dangerous diseases like measles and COVID-19.

## Repository Overview

* `paper/`: Rmarkdown and supporting files for creation of the manuscript. Must run supplemental materials notebook first.
* `supplement/`: Rmarkdown and supporting files for creation of the supplemental materials.
* `data/`: Data for all studies reported in paper.
* `code/`: Code scripts on which all Rmarkdown notebooks depend.
  * `../custom-structure-learning/`: Code implementing custom scoring algorithms and plotting tools used in Bayesian Network structure learning and inference.
* `outputs/`: stores outputs from running `supplement-main.Rmd`
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

3. Install required packages by running `install-packages.R` to make sure you have all packages needed.
4. Open and knit `supplement/supplement-main.Rmd` to generate supplement PDF and save files needed for reproduction of manuscript.
5. Open and knit `paper/paper-main.Rmd` to generate manuscript PDF.
