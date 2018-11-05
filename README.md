# PathFX

## High level description
PathFX is an interaction-network tool to search for the most relevant protein-protein interactions around a drug's target(s), and then analyzes for which phenotypes the network is enriched relative to the entire interaction network. The algorithm provides tabular results of associated phenotypes, a network file, and sources for phenotypic associations. Running the algorithm with phenotype clustering will produce a further summary file of associated phenotypes, grouped by semantic similarity using UMLS, and a clustergram figure showing associations between these clusters.

We describe the method further and testing and application of the method in Wilson et al, PLoS Comp Bio, (in press).

## Dependencies
PathFX was created in Python and is written in python version 3.6.0 and the necessary dependencies and versions are included as a .txt file in the rscs dir.

## Phenotype Clustering with UMLS
Enabling the phenotype clustering feature requires a local installation of the UMLS metathesaurus. Please register for UMLS access through the NLM and follow their guidance for installation. We used version 2017AA. Installation of UMLS is not required if you do not wish to use the phenotype clustering feature.
*** Ask Mike to share his UMLS installation document
