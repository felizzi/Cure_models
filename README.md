Implementing cure-rate models in oncology: a tutorial using R
====================

What are cure models?
----------

Cancer remains one of the most frequent non-communicable diseases and associated with a substantial health and economic burden around the globe [(Fidler _et al._, 2018)](https://www.ncbi.nlm.nih.gov/pubmed/28669281). In recent years, however, new treatments for cancer, so-called immunotherapies, have become available which help some patients to achieve long-term survival. These patients are considered "cured" in a statistical sense - their mortality is the same as that of the general population without cancer (also called "background mortality") [(Othus _et al._, 2018)](https://www.ncbi.nlm.nih.gov/pubmed/28408015).

These benefits of immunotherapies imply that standard statistical methods for evaluating survival need to be adapted, to account for the fraction of patients cured [(Chen, 2013)](https://www.ncbi.nlm.nih.gov/pubmed/24829754). *Cure models* were developed for exactly this purpose and offer a chance to model the cure proportion in a trial or population [(Othus _et al._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22675175).

Cure models are popular with statisticians but may not be familiar to other potential users, e.g. in health technology assessment or health economics. This tutorial, which is split into a paper (under development) and this repository, provides a (relatively) gentle introduction to cure models in order to make them more accessible and usable for non-technical audiences. This repo hosts the code that is discussed and illustrated in more detail in the accompanying paper, which is freely available from _tbc_.

Who can use the tutorial/code?
---------

Anyone who's interested! The code is provided under _tbc_ license.

How does it work?
---------
If you want to run and adapt the code on your machine, all you need is

+ A working [R](https://www.r-project.org/) installation
+ A couple of R packages, including (in alphabetical order)
  * _flexsurv_
  * _MASS_
  * _RCurl_
+ A copy of the code (obviously)
+ An account with the [Human Mortality Database](https://www.mortality.org/). This is required to get high-quality mortality data to approximate background mortality in the model. Registration with the HMD is free but you can, of course, also use alternative mortality data sources (a bit of formatting might then be required to use the code presented here). As a courtesy, we wrote some code to automate the download of HMD tables. Just insert your account details into the code.

Once you're good to go, there are two example scripts (both also discussed in the paper), one where the cure fraction is estimated from the available trial data and one where an external data source is used to estimate the cure fraction. The example scripts pull in all the required data and functions automatically - just make sure you don't change the folder structure.

Who do I contact if I have questions or suggestions?
---------

If you have questions or suggestions on how to improve the code, feel free to get in touch at federico.felizzi@roche.com.
