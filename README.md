# PathFX

## Update March 3, 2021 - Version 2.0 has officially been merged with master. 
PathFX version2.0 uses the same core algorithm, but the input data have been updated using data pulled from PheGeni, DisGeNet, OMIN, and Clinvar as of March of 2020. PathFX version2.0 uses DrugBank data from May of 2020. To use PathFXv1.0 as publised in Wilson et al, PLoS Comp Bio, 2018, please checkout the branch version1.0. 

## High level description
PathFX is an interaction-network tool to search for the most relevant protein-protein interactions around a drug's target(s), and then analyzes for which phenotypes the network is enriched relative to the entire interaction network. The algorithm provides tabular results of associated phenotypes, a network file, and sources for phenotypic associations. Running the algorithm with phenotype clustering will produce a further summary file of associated phenotypes, grouped by semantic similarity using UMLS, and a clustergram figure showing associations between these clusters.

We describe the method further and testing and application of the method in Wilson et al, PLoS Comp Bio, 2018.

## Dependencies
PathFX was created in Python and is written in python version 3.6.0 and the necessary dependencies and versions are included as a .txt file in the rscs dir.

## Running PathFX
We have included an example script 'run_PathFX.py' to demonstrate usage of algorithm parameters.

## Phenotype Clustering with UMLS, umls-interface.pl, and umls-similarity.pl
Enabling the phenotype clustering feature requires a local installation of the UMLS metathesaurus and two perl packages. Please register for UMLS access through the NLM and follow their guidance for installation (We used version 2017AA). We have also provided the configuration file for calling the umls-similarity.pl script and a protocol for installing UMLS and the perl packages. Installation of UMLS, and perl packages is not required if you do not wish to use the phenotype clustering feature.
Update 3/3/21 - We have created dockerized versions of PathFX + Phenotype Clustering from the original publication (e.g. PathFX version1.0 + UMLS 2017AA) and an updated version (e.g. PathFX version2.0 + UMLS 2020AA). Both containers are available through PathFXweb (https://www.pathfxweb.net/).


## PathFXweb
We have also created PathFXweb, a full-service web application for using the PathFX algorithm. Please visit https://www.pathfxweb.net/ for more information. You can run analyses via the webserver or download a dockerized version of PathFX + UMLS. Authentication through NLM is required for this download.

## Contact/Support
If you have questions, comments, or feedback about PathFX or PathFXweb, please contact pathfx_support@googlegroups.com.
