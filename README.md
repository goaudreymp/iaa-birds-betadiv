# Birds and Barriers: present and past seas are dominant correlates of avian turnover in the Indo-Australia Archipelago
This folder includes data and code to run the analysis for the manuscript in progress Prasetya _et al_., _"Birds and Barriers: present and past seas are dominant correlates of avian turnover in the Indo-Australia Archipelago"_.

The associated code files includes:
- Rproject for easy management of files and code (`iaa-birds-betadiv.Rproj`)
- R script to run the code analysis (`iaa-birds-betadiv.R`)

Also included is a 'master' file with the collated list of all species and presence/absence data:
- `Prasetya_birdslistmaster.xlsx`

Below is a description for each data file included in this folder that is required to run the above code:

## Beta-diversity 
- Presence/absence matrix and species list matching Jetz (2009) phylogeny name
  - Species-level (`spcommatrix.csv`)
  - Genus-level (`gencommatrix.csv`)
  - Family-level (`famcommatrix.csv`)

## Connectivity Analyses
- calculated connectivity values 
  - Species-level (`edges.sp.csv`)
  - Genus-level (`edges.g.csv`)
  - Family-level (`edges.f.csv`) 

## Correlates of Turnover
- latitude longitude of each area (`region.latlong.csv`)
- environmental data for each area obtained from WorldClim (`biomeans.csv`)
- matrix of past historical connectivity to continental shelfs (`shelfmatrix.csv`)
- matrix of minimum geographic distance (`region.shp.distance.csv`)
- matrix of total land area between area pairs (`area.mx.csv`)
- combined results of MMRR analysis for plotting (`all.dist.matrix.melted.csv`)