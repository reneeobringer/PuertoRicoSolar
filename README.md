# PuertoRicoSolar
Data and code for an analysis of the relationship between solar energy and local weather in Puerto Rico. The results from the analysis are currently under review.

A permanent version of this repository has been archived via [Zenodo](https://doi.org/10.5281/zenodo.17781944).


Two categories of data were collected: solar energy generation data and local weather data. The solar energy generation data (folder: `SolarData`) includes potential solar energy generation data for 76 Puerto Rican municipalities obtained through PVGIS. The weather data (folder: `WeatherData`) were obtained from the Daymet gridded dataset. All data were collected 2023 and 2025.

There are two scripts associated with this project. The first, `solarpowercalcs.ipynb`, is a Jupyter Notebook script developed in Google Colab and last ran on 17 November 2025. The second, `PRsolarpower.R`, was developed in R version 4.5.1 and last run on 17 November 2025. The associated Rdata files can be found in the `rdatafiles` folder. In order to run the script in its entirety, the following R packages are required, with the versions used in paraentheses: 

* corrplot (v0.95)
* cowplot (v1.2.0)
* ggplot2 (v4.0.0)
* ggspatial (v1.1.10)
* humidity (v0.1.5)
* lubridate (v1.9.4)
* randomForest (v4.7.1.2)
* readxl (v1.4.5)
* reshape2 (v1.4.4)
* sf (v1.0.22)
* stringr (v1.5.2)
* tidyverse (v2.0.0)
* transport (v0.15.4)

To run both scripts, users need to update the directory paths to the folder downloaded or cloned from this repository. This change can be made on lines 38-42 of `PRsolarpower.R` and within the second code block of `solarpowercalcs.ipynb`. The R script can then be run sequentially or users can choose to run different sections, provided users load the rdata files in the `rdatafiles` folder prior to running each section. 
