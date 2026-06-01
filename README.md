# Beyond the traditional three: A comprehensive model-fit assessment of phenotypic evolution in the fossil record


__Article:__ Unpublished

__Authors:__ Vilde Bruhn Kinneberg<sup>*1†</sup>, Marion Thaureau<sup>*1‡</sup> and Kjetil Lysne Voje<sup>1</sup>

<sup>*</sup>VBK and MT both contributed equally to this work


__Affiliation:__ <sup>1</sup>Evolution and Paleobiology, Natural History Museum, University of Oslo

__Contact:__ <sup>†</sup>v.b.kinneberg@nhm.uio.no, <sup>‡</sup>marion.thaureau@nhm.uio.no

__Journal:__ NA

__Year:__ 2026  

__Abstract:__ Models of phenotypic evolution help us interpret patterns and rates of evolution in empirical data. In the fossil record, three models are traditionally used to study within-lineage evolution: stasis, random walk, and directional evolution. Previous studies have established that random walk and stasis models commonly provide the best fit to empirical time series. Still, the traditional three models represent only a small subset of the models available to describe trait evolution in the fossil record. In this study, we fit nine single-mode and sixteen mode-shift models to a compilation of 594 fossil time series to assess how often the different models best describe phenotypic evolution. We find that stasis and the unbiased random walk models describe a substantial proportion of the time series, but about a quarter are better described by evolutionary models other than the traditional three. When also considering   mode-shift models, the proportion of time series best described by models other than the traditional three increases to more than 50%. Our results indicate a rich diversity of trait dynamics within lineages in the fossil record that the traditional three models cannot fully capture.


__Info:__ This repository contains scripts and data used for analyses in the publication.

__Responsibility:__ VBK and MT are responsible for analyses of data. KLV has contributed with ideas and comments to the analyses. The time series data are downloaded from the [PETS database](https://pets.nhm.uio.no/). 

__Files__ 

_data –_ this folder conatins data loaded in the scripts used for analyses.

_scripts –_ this folder conatins scripts used to performe analyses. All scripts are commented so that it should be possible to follow the instructions in the scripts to produce the results from the article.
<ul>
  <li>single_mode.R _single_mode.R_ runs the analyses for the single-mode models section of the article.</li>
  <li>stats_single_mode.R runs regression analyses and plotting on the single-mode data generated in single_mode.R.</li>
  <li>mode_shift.R _single_mode.R_ runs the analyses pluss summary statistics for the mode-shift models section of the article.</li>
  <li>stats_mode_shift.R runs regression analyses and plotting on the mode-shift data generated in mode_shift.R.</li>
  <li>shift_OU_adequacy.R is a script run on a HPC cluster for the OU part of adequacy testing in the mode-shift part of the article.</li>
  <li>delta_aicc_gap.R runs analyses and plotting for the delta AICc gap part of the article.</li>
  <li>functions.R contains R functions loaded in all the above scripts.</li>
</ul>
