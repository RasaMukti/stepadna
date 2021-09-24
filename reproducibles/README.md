
The purpose of the files and data in this folder is to be able to reproduce results described in the paper "Modelling the spatiotemporal spread of beneficial alleles using ancient genomes"

doi: https://doi.org/10.1101/2021.07.21.453231

## Data
The file "tyr_data.csv" contains genotype information for rs1042602(A) allele. The ancient samples are obtained from Allen Ancient DNA Resource data [1].

and modern samples are from Human Genome Diversity Panel [2].

The two files in the LCT_data folder named "ancient_data.xlsx" and "Lactase_all_GP.txt" contain genotype information for rs4988235(T) allele for ancient and modern day samples, respecitvely. The data set has been compiled in Segurel et al. (2020) [3].

The geographical boundaries used in the study are between 30°N to 75°N, and between 10°W and 80°E. All samples outside of this region have been filtered out.

## Code
In order to obtain parameter estimates for rs1042602(A) allele run stepadna.R using command:
```
Rscript stepadna.R -f tyr_data.csv -a 26361 -o tyr_output.csv -i 50 -c 50
```
Parameter estimates for rs4988235(T) allele can be obtained running the command:
```
Rscript stepadna_lct.R -a 7441 -o lct_output.csv -i 50 -c 50
```
The "-a" flag specifies the age of the allele in years. The commands will be run initiating 50 points in the simulated annealing algorithm parallelised across 50 cores.

The results described in the afore-mentioned paper have been produced using R verision 3.6.

The scripts are dependant on the following libraries:
deSolve, ReacTran, RColorBrewer, geosphere, fields, stats4, GEOmap, geomapdata, plsgenomics, paramtest, bbmle, parallel, maps, mapdata, mapplots, rworldmap, scales, tidyverse, FME, optparse

These will be installed automatically upon running the script, unless they have already been installed.

## References
[1] AADR: https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data

[2] Bergstrom A, McCarthy SA, Hui R, Almarri MA, Ayub Q, Danecek P, Chen Y, Felkel S, Hallast P, Kamm J, et al. (2020). Insights into human genetic variation and population history from 929 diverse genomes. Science, 367(6484)

[3] Segurel L, Guarino-Vignon P, Marchi N, Lafosse S, Laurent R, Bon C, Fabre A, Hegay T, & Heyer E (2020). Why and when was lactase persistence selected for? insights from central asian herders and ancient dna. PLoS Biology, 18(6):e3000742
