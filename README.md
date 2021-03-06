# Mortality attributable to socioeconomic inequality in England

This repository includes code associated with the article:

Lewer D, Jayatunga W, Aldridge RW et al.  Premature mortality attributable to socioeconomic inequality in England between 2003 and 2018: an observational study. Lancet Public Health 2020;5(1):e33–41. https://doi.org/10.1016/S2468-2667(19)30219-1

Detailed results tables are available in the [supplementary appendix](https://www.thelancet.com/cms/10.1016/S2468-2667(19)30219-1/attachment/81055507-b222-435d-bf59-44848a61e28f/mmc1.pdf) and at [UCL](https://doi.org/10.14324/000.ds.10086658); and there is an [interactive map](https://public.tableau.com/profile/rob.aldridge#!/vizhome/MATI_19_11_25/MATI_dashboard). 

The repository includes:

* [Code](https://github.com/danlewer/masi/blob/main/example_masi_calculation.R) providing an example of how MASI is estimated, using publicly available data (that is also included in the repository and read directly by the script)
* [Code](https://github.com/danlewer/masi/blob/main/pie_function.R) showing how the 'exploded pie chart' (figure 3 in the published article) is drawn

# Abstract

### Background
Low socioeconomic position is consistently associated with increased risk of premature death. The aim of this study is to measure the aggregate scale of inequality in premature mortality for the whole population of England.
### Methods
We used mortality records from the UK Office for National Statistics to study all 2 465 285 premature deaths (defined as those before age 75 years) in England between Jan 1, 2003, and Dec 31, 2018. Socioeconomic position was defined using deciles of the Index of Multiple Deprivation: a measure of neighbourhood income, employment, education levels, crime, health, availability of services, and local environment. We calculated the number of expected deaths by applying mortality in the least deprived decile to other deciles, within the strata of age, sex, and time. The mortality attributable to socioeconomic inequality was defined as the difference between the observed and expected deaths. We also used life table modelling to estimate years-of-life lost attributable to socioeconomic inequality.
### Findings
35·6% (95% CI 35·3–35·9) of premature deaths were attributable to socioeconomic inequality, equating to 877 082 deaths, or one every 10 min. The biggest contributors were ischaemic heart disease (152 171 excess deaths), respiratory cancers (111 083) and chronic obstructive pulmonary disease (83 593). The most unequal causes of death were tuberculosis, opioid use, HIV, psychoactive drugs use, viral hepatitis, and obesity, each with more than two-thirds attributable to inequality. Inequality was greater among men and peaked in early childhood and at age 40–49 years. The proportion of deaths attributable to inequality increased during the study period, particularly for women, because mortality rates among the most deprived women (excluding cardiovascular diseases) plateaued, and for some diseases increased. A mean of 14·4 months of life before age 75 years are lost due to socioeconomic inequality.
### Interpretation
One in three premature deaths are attributable to socioeconomic inequality, making this our most important public health challenge. Interventions that address upstream determinants of health should be prioritised.
## Funding
National Institute of Health Research; Wellcome Trust.

# License
This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

Please cite the published article:

Lewer D, Jayatunga W, Aldridge RW et al.  Premature mortality attributable to socioeconomic inequality in England between 2003 and 2018: an observational study. Lancet Public Health 2020;5(1):e33–41. https://doi.org/10.1016/S2468-2667(19)30219-1

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
