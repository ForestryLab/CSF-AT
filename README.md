# CSF-AT Shiny web-app
# Background
<div>
<img src="https://www3.unimol.it/assets/images/unimol/images/header/unimol_on.svg", alt="logo_UNIMOL" width="160" align="right" style="margin: 0 0 50px 50px;" />
<img src="https://forwards-project.eu/wp-content/themes/forwards_theme/images/FWD-LOGO-RGB_MAIN_COLOR.svg", alt="logo_CSF" width="150" align="right" style="margin: 0 0 50px 50px;" />
The Climate Smart Forestry Assessment Tool (CSF-AT) is a user-friendly Shiny web-app designed to assess and quantify the CSF management of a forest stand using forest inventory data. The app allows  users to upload standardized datasets, calculate the composite climate-smart index and its three components, and explore the results through interactive visualizations.
  
The project is funded by **FORWARDS (ForestWard Observatory) EU-funded project**.
</div>

# Launching the web-app 
CSF-AT  app is available with the following link: https://ricercaforestale-unimol.shinyapps.io/CSF-AT/

## Sections of web-app
The web-app consists of three main panel: 
1.	App Description — This panel outlines the CSF concept, the objectives of the application, the methodological framework implemented, and the operational workflow, providing users with an overview of the tool’s functionality and expected outputs.
2.	Smartness Calculation — This is the operational core of the tool. Users can upload one or multiple standardized datasets to compute the Composite Climate-Smart Index (ICSF) and its three pillars at both plot and stand levels. In this panel, the methodological framework developed for ICSF assessment is fully digitalized, including data normalization, indicator weighting, and hierarchical aggregation. All resulting outputs can be exported in .csv format. 
3.	Graphic Comparison — This panel provides an interactive ternary plot displaying the proportional contributions of the three pillars for each observation, enabling rapid visual comparison across stands, years, or management contexts.

# Methodological framework behind the CSF assessment
The CSF assessment is based on a weighted hierarchical approach designed to evaluate forest smartness by integrating ecological, social, and management-related information into a single composite index. The full list of indicators included in the CSF assessment, along with possible calculation methods, is available in the repository folder "tables -> indicators_description".

<img src="/www/CSF_framework.png" width="600" align="center"/>

The process starts from a set of predefined indicators describing key structural and functional attributes of forest stands. Indicator values are first normalized using min–max scaling, so that all variables are expressed on a common scale ranging from 0 (least desirable condition) to 1 (most desirable condition). This step ensures comparability among indicators with different units and ranges.
Normalized indicators are then weighted according to their relative importance, as defined through expert judgment using the Analytic Hierarchy Process (AHP). The weighting reflects the contribution of each indicator to the three main CSF pillars: Mitigation, Adaptation, and Socio-economic dimension.
Finally, indicators are aggregated hierarchically. First, weighted indicators are combined within each pillar to obtain pillar-specific smartness scores. These pillar scores are then aggregated to compute the **Composite Climate-Smart Forestry Index**, which represents the overall smartness of the forest stand (Alfieri et al., 2024).

# User Interface and Tutorial
1.	After launching the app, the user defines the number of datasets to analyze.
2.	The user uploads one .csv file per forest stand, each containing at least two sample plots. The uploaded .csv files must follow a defined structure with pre-set column headers. An example file is provided in the repository folder  “tables -> examples_datasets” for reference.
3.	The user clicks “Calculate Smartness” to compute the CSF indices.
d.	Results are displayed at plot level (original indicators + csf, csf_mit, csf_adp, csf_sd).
e.	Results are summarized at stand level (mean smartness across plots).
f.	If yearly data are provided, the user can select the year to view corresponding results.
4.	The user can click “Download all” to export a .zip archive containing plot-level and stand-level CSV files.
5.	The user can click “Reset App” to clear all inputs and start over.

<img src="/www/smartness_calculation.png" width="600" align="center"/>

7.	To explore smartness visually, the user switches to the “Graphic Comparison” window:
d.	View the interactive ternary plot showing contributions of mitigation, adaptation, and social dimension.
e.	Hover over points to see dataset name, year, number of plots, and absolute smartness scores.
f.	Zoom, pan, and export the graph as a .png file.

<img src="/www/graphic_comparison.png" width="600" align="center"/>

# References
Alfieri, D., Tognetti, R., & Santopuoli, G. (2024). Exploring climate-smart forestry in Mediterranean forests through an innovative composite climate-smart index. Journal of Environmental Management, 368, 122002.
