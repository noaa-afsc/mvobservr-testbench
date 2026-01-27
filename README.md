# Testing the Multivariate Observer Effects Framework
## Collaborators
The primary author(s) responsible for maintaining this repository are:
* Craig H. Faunce (craig.faunce@noaa.gov)
* Debra Duarte (debra.duarte@noaa.gov)
* Geoff Mayhew (geoff.mayhew@noaa.gov)

# Description
Addressing observer effects, wherein the act of observation has influence on the phenomenon being observed, is of great importance to science because its presence indicates that the results from data collected from observation (the sample) are biased and cannot be used to infer the properties of unobserved nature (the population). The presence of an observer effect in partial coverage fisheries means that the data from monitored trips (the sample) are not representative of the entire fleet (the population), and this bias can have broad implications to catch accounting and stock assessments.

Analysts of the Fisheries Monitoring and Analysis Division (AFSC) have developed a model-based method to test for observer effects in the catch of multiple species from the partial coverage fleet (Appendix A; NMFS and AKRO, [2025](https://meetings.npfmc.org/CommentReview/DownloadFile?p=0fbf1a38-b5bf-4029-9dbc-6c009dc7ab38.pdf&fileName=2024%20Observer%20Program%20Annual%20Report.pdf). The method tests different multivariate generalized linear models (MvGLM) with permutation to evaluate the similarity of species abundances from landed catch between trips that were monitored with EM or observers and those that were not. Because the method employs a MvGLM to detect observer effects, it is abbreviated for convenience as MOE (Multivariate Observer Effects), and the model as mvglm_obs.

Despite its development, the relative performance of the method has not been tested. Its relative gains over other methods are based solely on reported advantages from other sources used in its development (Wharton et al. [2012](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00127.x)).  

# Purpose
The work aims to iterate the work of Duerte and Cadrin [(2004)](https://doi.org/10.1016/j.fishres.2024.107000) to include MOE in its evalution of the relative performance of methods used to measure observer effects in fisheries.

The repo [mvobservr](https://github.com/noaa-afsc/mvglm-obs#) associated with conducting MOE will be used to source the most up to date code for this project.  Potential comparisons to this method are slim, because MOE uses a matrix as a response.  PERMANOVA is one for sure however.
 
# Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project content is provided on an "as is" basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
