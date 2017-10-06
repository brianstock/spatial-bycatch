## Which species distribution model best predicts fisheries bycatch?

This repository contains code to run the analyses in:

> Stock BC, Ward EJ, Eguchi T, Jannot JE, Thorson JT, Feist BE, and Semmens BX. "Random forests outperform other species distribution models for spatiotemporal fisheries bycatch prediction."

Because the fisheries observer datasets we used are confidential ([WCGOP](https://www.nwfsc.noaa.gov/research/divisions/fram/observation/data_collection/manuals/2017%20WCGOP%20Training%20Manual%20Final%20website%20copy.pdf), [HILL](http://www.nmfs.noaa.gov/pr/interactions/fkwtrt/meeting1/handouts/observer_manual.pdf)), here we perform the same analyses using the publicly available [West Coast Groundfish Trawl Survey](https://www.nwfsc.noaa.gov/research/divisions/fram/groundfish/bottom_trawl.cfm).

### Two variations
  1. Simple (*if you only care to see how the models are implemented*): starts with a prepared.RData file, no data processing, no cross validation, only 1 species (darkblotched rockfish), only 5 years.
  2. Full: includes data processing, cross validation, 3 species, all years.

### 1) Simple

##### 1_simple

Get straight to running the models, starting with a prepared `.RData` file:
  - skips data processing and cross validation
  - 1 species (darkblotched rockfish)
  - 5 years

### 2) Full
  - data processing
  - cross validation
  - 3 species
  - all years
  
##### [2a_process_survey](https://rawgit.com/brianstock/spatial-bycatch/master/2a_process_survey.html)

Details how we downloaded and processed the data, including preparing the environmental covariates. We save the output of this script as `wcann_processed.RData`, so if you are not interested in this step, you can skip ahead to run the spatial models.

##### [2b_run_models](https://rawgit.com/brianstock/spatial-bycatch/master/2b_run_models.html)

Sets up and runs 10-fold cross validation for all spatial models. *Note that in our manuscript we used 5-fold cross validation repeated 10x. We use 1x 10-fold CV here to save model run time.* We save the output of this script as `wcann_models_finished.RData`, so if you are not interested in running this step (takes ~24 hours on my computer), you can see the model output.

##### [2c_output](https://rawgit.com/brianstock/spatial-bycatch/master/2c_output.html)

Creates output from model runs, replicating the figures from our analysis:
  * Fig. 1: Maps of effort and catch (raw data)
  * Fig. 2: Compare model performance with boxplots of AUC and RMSE
  * Fig. 3: Compare the reduction in bycatch-to-target ratio
  * Fig. 4: Maps of predicted density (mean) and variablity (log CV)
  * Fig. 5: Visualize covariate effects for GMRF and RF
  * Fig. 6: Map the GMRF spatial random field by year (one for each year, only 1 species)
  * Fig. S3: Map the GMRF spatial random field (one across all years, for each of 3 species)
