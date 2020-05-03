Implementing cure models in oncology: a tutorial using R
====================

What are cure models?
----------

Cancer remains one of the most frequent non-communicable diseases and associated with a substantial health and economic burden around the globe ([Fidler _et al._, 2018](https://www.ncbi.nlm.nih.gov/pubmed/28669281)). In recent years, however, new treatments for cancer, including immunotherapies and targeted therapies, have become available which help some patients to achieve long-term survival. These patients are considered "cured" in a statistical sense - their mortality is the same as that of the general population without cancer (also called "background mortality") ([Othus _et al._, 2018](https://www.ncbi.nlm.nih.gov/pubmed/28408015)).

The benefits of these novel therapies imply that standard statistical methods for evaluating survival need to be adapted in order to account for the fraction of patients cured ([Chen, 2013](https://www.ncbi.nlm.nih.gov/pubmed/24829754)). *Cure models* were developed for exactly this purpose and allow to model the cure proportion for a clinical trial or a real-world population ([Othus _et al._, 2012](https://www.ncbi.nlm.nih.gov/pubmed/22675175), [Lambert _et al._, 2007](https://www.ncbi.nlm.nih.gov/pubmed/17021277)).

Cure models are popular with statisticians but may not be familiar to other potential users, e.g. in health technology assessment or health economics. This tutorial, which is split into a paper (under development) and this repository, provides a (relatively) gentle introduction to cure models in order to make them more accessible and usable for non-technical audiences. This repo hosts the code that is discussed and illustrated in more detail in the accompanying paper, which is freely available from this Github repository.

Who can use the tutorial/code?
---------

Anyone who's interested! The code is provided under a CC BY-NC 4.0 license.

How does it work?
---------
If you want to run and adapt the code on your machine, you need

+ A working [R](https://www.r-project.org/) installation
+ A couple of R packages, including (in alphabetical order)
  * _flexsurv_
  * _MASS_
  * _RCurl_
+ An account with the [Human Mortality Database](https://www.mortality.org/). This is required to get high-quality mortality data to approximate background mortality in the model. Registration with the HMD is free but you can, of course, also use alternative mortality data sources (a bit of formatting might then be required to use the code presented here). As a courtesy, we wrote some code to automate the download of HMD tables. Just insert your account details into the code.

Once you're good to go, there are two example scripts (both also discussed in the paper), one where the cure fraction is estimated from the available trial data and one where an external data source is used to estimate the cure fraction. The example scripts pull in all the required data and functions automatically - just make sure you don't change the folder structure.

Folders and files are as follows

+ Analysis scripts
  * `estimate_cure_brim3.R`: Analysis for an 'uninformed' approach, where cure fractions are *estimated* from clinical data (BRIM3 trial)
  * `input_cure_cobrim.R`: Analysis for an 'informed' approach, where cure fractions are provided as an *input* to the estimation based on clinical data (coBRIM trial)
  * `mortality_table_wrap.R`: Prepares mortality data for the uninformed (BRIM3) and informed (coBRIM) analyses
+ Helper functions (in 'functions/')
  * `funs_hazard.R`: Functions to estimate mortality hazard and survival for each and across observation(s) - see the paper for their explanation and derivation
  * `funs_likelihood.R`: Functions to estimate the likelihood functions - see the paper for their explanation and derivation
  * `funs_load_mort_table.R`: Functions to download mortality data from the Human Mortality Database
  * `funs_long_term_survival.R`: Helper functions to estimate hazard and survival rates
+ Data (in 'data/')
  * `app_aut.txt`: Authentication files, change this to your specific authentication requirements
  * `country_list.csv`: List of countries relevant to the BRIM3 and coBRIM trials
  * 'mortality/': Gender-specific mortality data from the HMD, for countries relevant to the analyses
  * 'trials/':
    - `brim3_simulated.csv`: Extract of BRIM3 data, anonymized with random Gaussian noise added to demographic data
    - `cobrim_simulated.csv`: Extract of coBRIM data, anonymized with random Gaussian noise added to demographic data

Whom do I contact if I have questions or suggestions?
---------

If you have questions or suggestions on how to improve the code, feel free to contact Federio Felizzi at firstnameDOTlastnameATrocheDOTcom.
