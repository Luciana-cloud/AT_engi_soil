AT_engi_soil (atrazine degradation model for engineered systems and extension for soil simulations)

This respository includes a collection of models to simulate the degradation of atrazine and formation of the metabolite hydroxyatrazine in engineered systems (chemostat/retentostat) and soil systems.

Engineered systems
The models for engineered systems can be found in the folder "chemo_retent_models".
The first two model approaches use a simple bacterial growth based on hydroxyatrazine. One is described by Monod-kinetics (variant_1.m), the other by thermodynamically constrained growth (variant_2.m) (1).
The third model is a complex degradation model including bacterial differentiation based on substrate concentration (variant_complex.m).

Soil systems
The soil system models can be found in the folder "soil_models".
The model extension for soil simulations is based on variants 1 and 2 from the engineered systems, and on adding a leaching flux.
Additionally, we included one model version without any biological degradation of atrazine.

The experimental data used for model calibration of the engineered systems can be found in the folder "experimental data". In the same folder we included soil sorption data, as well as the regression to determine the sorption parameters.
For further information about the data and experimental setup, please see references (2), (3), (4) and (5).

We provided several MATLAB live script files:
The folders "calibration" and "sensitivity_anaylysis" contain different examples for calibration and sensitivity analysis respectivly. We presented one example per model variant, but more information for additional examples is given.
The figures of the paper can be found in the folder "figures_paper".

Data files over 1 GB in size are not included, but are available upon request.

Contact:

In case of questions, please contact Luciana Chavez (lucianagoku@hotmail.com)

References:

1. Desmond-Le Quéméner, E., & Bouchez, T. (2014). A thermodynamic theory of microbial growth. The ISME Journal, 8(8), 1747–1751. https://doi.org/10.1038/ismej.2014.7

2. Ehrl, B. N., Kundu, K., Gharasoo, M., Marozava, S., & Elsner, M. (2019). Rate-Limiting Mass Transfer in Micropollutant Degradation Revealed by Isotope Fractionation in Chemostat. Environmental Science & Technology, 53(3), 1197–1205. https://doi.org/10.1021/acs.est.8b05175

3. Gharasoo, M., Ehrl, B. N., Cirpka, O. A., & Elsner, M. (2019). Modeling of Contaminant Biodegradation and Compound-Specific Isotope Fractionation in Chemostats at Low Dilution Rates. Environmental Science & Technology, 53(3), 1186–1196. https://doi.org/10.1021/acs.est.8b02498

4. Kundu, K., Marozava, S., Ehrl, B., Merl-Pham, J., Griebler, C., & Elsner, M. (2019). Defining lower limits of biodegradation: atrazine degradation regulated by mass transfer and maintenance demand in Arthrobacter aurescens TC1. The ISME Journal, 13(9), 2236–2251. https://doi.org/10.1038/s41396-019-0430-z

5. Kundu, K., Weber, N., Griebler, C., & Elsner, M. (2020). Phenotypic heterogeneity as key factor for growth and survival under oligotrophic conditions. Environmental Microbiology, n/a(n/a). https://doi.org/10.1111/1462-2920.15106
