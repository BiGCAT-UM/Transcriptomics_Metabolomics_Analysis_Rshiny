## A workflow for functional analysis of transcriptomic and metabolomic data in R shiny 

This workflow, developed in R Shiny, includes differential gene expression analysis, statistical analysis of metabolomics data, as well as identifier mapping and pathway enrichment analysis for both transcriptomics and metabolomics data followed by integration of this data through network analysis to identify disease-related processes and visualization of multi-omics data. A publicly available (https://ibdmdb.org/) gut-transcriptomic and stool-metabolomic dataset of the gut microbial ecosystem in inflammatory bowel diseases was used to test the proposed workflow.

#### Transcriptomics analysis  
[1-Data preprocessing](/1-data_filtering)<br /> 
[2-Differential gene expression analysis](/2-data_normalization)<br />
[3-Identifier mapping](/3-differential_gene_expression_analysis)<br />
[4-Identifier mapping](/4-identifier_mapping)<br />
[5-Gene Pathway analysis (ORA)](/5-pathway_analysis/)<br />
[6-Heatmap creation](/6-create_heatmap/)<br />
[7-Network analysis](/7-network_analysis)<br />

#### Metabolomics analysis  
[7-Data preprocessing](/7-metabolite_data_preprocessing)<br />
[8-Significantly changed metabolites analysis](/8-significantly_changed_metabolites_analysis)<br />
[9-Identifier mapping](9-metabolite_identifier_mapping)<br />
[10-Metabolite Pathway Analysis (ORA)](/10-metabolite_pathway_analysis)<br />

#### Multi-omics visualization
[11-Pathway selection](11-pathway_selection)<br />
[12-Multi-omics visualization](12-multiomics_visualization)<br />

The workflow is an example of how to bring together different software tools and methods to analyze transcriptomics and metabolomics data to reveal the underlying mechanism behind IBD disease as a use case study.
![workflow pptx](https://user-images.githubusercontent.com/65600609/210248763-ae312fec-4df9-43f0-9cc6-c995629dd2c2.jpg)

## Setup and Preparation
The setup of this project has been tested with:
- OS Windows 10, R-studio 2021.09.02, R 4.1.3.
- OS Linux (Debian), R-studio 2022.02.2, R 4.2.

### 1. Download required software tools
* [R](https://cran.r-project.org/bin/windows/base/)
* [RStudio](https://www.rstudio.com/products/rstudio/download/#download)
* [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
* [Cytoscape](https://cytoscape.org/)

### 2. Install required apps in Cytoscape
1. Open Cytoscape.
2. Go to the app manager: Apps -> App Manager
3. Install the following apps:
    * WikiPathways
    * stringApp
    * clusterMaker2

Detailed instructions can be found [here](https://bigcat-um.github.io/Transcriptomics_Metabolomics_tutorials/pages/prep).

### 3. Run the app
To run the app in RStudio, click on **"Run App"** in the top right corner when having either the `ui.R`, `server.R`, or `global.R` file open in the RStudio window.

If this is not possible, please try running the following commands to start the app:
```r
# Install the shiny package
install.packages("shiny")

# Load the shiny package
library(shiny)

# Run the shiny app
runApp(..PATH..)
```

## Acknowledgment
This research was undertaken by Maastricht University (UM, Netherlands), a beneficiary in FNS-Cloud, which has received funding from the European Union’s Horizon 2020 Research and Innovation programme (H2020-EU.3.2.2.3. – A sustainable and competitive agri-food industry) under Grant Agreement No. 863059.
