# stepadna
Spatio-temporal modelling of beneficial allele spread using low-coverage Ancient DNA

Input file should be a .csv file, where each row corresponds to a sample and contains the following columns:
- age of the sample in years BP named "age"
- degrees latitude named "latitude"
- degrees longitude named "longitude"
- pseudohaploid genotype named "genotype", where 0 indicates ancestral allele state and 2 stands for derived allele.

Example of an input file:

   | age | latitude | longitude | genotype
   | --- | --- | --- | ---
s1 | 6590 | 34.3 | 44.5 | 0


The input also asks for the age of the allele.
Optional parameters are
- number of cores (default is 1)
- number of initial points for the simulated annealing algorithm (defaults to 2)

The output is a .csv file containing the inferred:
- selection coefficient
- longitudinal diffusion coefficient
- latitudinal diffusion coefficient
- longitudinal advection coefficient
- latitudinal advection coefficient
- longitude of the origin of the allele
- latitude of the origin of the allele

 Additionally, a .gif file with the inferred allele dynamics in produced.

