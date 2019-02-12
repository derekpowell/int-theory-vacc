# Modeling intuitive theories surrounding vaccination decisions

Authors: Derek Powell, Kara Weisman, & Ellen Markman

Repository related to psychological studies exploring methods for capturing intuitive theories as Bayesian networks using structure learning techniques.

Manuscript in preparation, title may change. Nothing here is final. please contact if you have any questions about this project or code. For related work, see our [CogSci paper](https://mindmodeling.org/cogsci2018/papers/0183/0183.pdf) and [this repo](https://github.com/derekpowell/vaccbeliefs-cogsci2018).

_Related abstract_

How can we leverage the cognitive science of lay theories to inform interventions aimed at correcting misconceptions and changing behaviors? Focusing on the problem of vaccine skepticism, we identified a set of 14 beliefs we hypothesized would be relevant to vaccination decisions. We developed reliable scales to measure these beliefs across a large sample of participants (n = 1130) and employed state-of-the-art graphical structure learning algorithms to uncover the relationships among these beliefs. This resulted in a graphical model describing the system of beliefs relevant to childhood vaccinations, with beliefs represented as nodes and their interconnections as directed edges. This model sheds light on how these beliefs relate to one another and can be used to predict how interventions aimed at specific beliefs will play out across the larger system. Moving forward, we hope this modeling approach will help guide the development of effective, theory-based interventions promoting childhood vaccination.

## Project notes

2/12/19, 12:15 PM

I'm moving analysis notebooks into their own `notebooks/` folder. There's too much data/code overlap between study 1 and study 2, it didn't make sense to organize around studies. Accordingly, the `Study1/` and `Study2/` folders should be considered deprecated, but I'm keeping them around as I transition to a new organizational format.

In this same spirit, I need to do some work to refactor the `scripts/` folder and related files.
