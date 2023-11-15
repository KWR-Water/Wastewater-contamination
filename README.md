## Wastewater-contamination

Pathogen intrusion in drinking water systems can pose severe health risks to human populations. Here we present a benchmark model to assess the health impacts of gross wastewater contamination in drinking water networks and investigate the mitigation potential of chlorination. The model considers organic matter degradation, chlorine decay mechanisms, diverse pathogen inactivation kinetics, as well as stochastic water demands. 
We study wastewater infiltration events that can occur anywhere within a drinking water distribution network, focusing on three pathogens: Enterovirus, Campylobacter, and Cryptosporidium. Their distinct chlorine resistance and infectivity patterns are applied to both chlorinated and non-chlorinated settings, and combined with pathogen-specific dose-response models within the Quantitative Microbial Risk Assessment framework. Synthetic household-level water demand time series were used to isolate the tap water end-use and calculate the infection risk (exposure via ingestion). Several consumption events were distributed throughout the day per individual. Model outcomes indicate that while disinfection aids mitigation, gross contaminations can still lead to infections due to chlorine depletion at the point of contamination. 
In our scenarios, chlorine-susceptible pathogens infected 0.78–26.6 % of the downstream population, while chlorine-resistant ones impacted the entire downstream population. 
Contamination location plays an important role since it affects the dilution of the pathogen concentration and thus exposure. Hydraulic uncertainty had a limited influence on infection risk. Furthermore, the model outcome is sensitive to the initial pathogen concentration, duration of the contamination, and the inactivation rate in chlorinated systems. The model further indicates that the time window for effective mitigation of the size of a waterborne outbreak is short (within hours).

## Requirements
Please install the following software before use (see links for installation instructions):
* [Matlab](http://www.mathworks.com/)
* [EPANET-Matlab-Toolkit] (https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit)

## Instructions
* Run [Run_benchmark.m] in MATLAB to perform simulations
* Run [Display_results.m] to see the results

## Contributors
* [Sotirios Paraskevopoulos](https://github.com/Sotireas), [KWR Water Research Institute](https://www.kwrwater.nl/en/)
* [Stelios Vrachimis](https://github.com/SteliosVr), [KIOS Research and Innovation Center of Excellence, University of Cyprus](http://www.kios.ucy.ac.cy/)
* [Marios Kyriakou](https://github.com/Mariosmsk), [KIOS Research and Innovation Center of Excellence, University of Cyprus](http://www.kios.ucy.ac.cy/)
