# stepadna
Spatio-temporal modelling of beneficial allele spread using low-coverage Ancient DNA

You can read about the method in more detail in: doi: https://doi.org/10.1101/2021.07.21.453231

### Input
Input file should be a .csv file, where each row corresponds to a sample and contains the following columns:
- age of the sample in years BP named "age"
- degrees latitude named "latitude"
- degrees longitude named "longitude"
- pseudohaploid genotype named "genotype", where 0 indicates ancestral allele state and 2 stands for derived allele.

Example of an input file:

 sampleID  | age | latitude | longitude | genotype
 --- | --- | --- | --- | ---
s1 | 6590 | 35.3 | 44.5 | 0
s2 | 3200 | 54 | 51.8 | 0
s3 | 0 | 52.2 | 18.6 | 2


Input parameters:
```
-f --file – input file containing the data
-a --age – age of the allele in years 
-o --out - output file name (defaults to out.csv)
-i --init - number of points to initialize for the simulated annealing algorithm (defaults to 2)
-c --cores - number of cores to use (defaults to 1)
```
Example to run the programme from the command line for an allele that is 29000 years old
with 20 initial points in the simulated annealing algorithm distributed across 20 cores:
```
Rscript stepadna.R -f data.csv -a 29000 -o result.csv -i 20 -c 20
```
### Output
The output is a .csv file containing the inferred:
- selection coefficient
- longitudinal diffusion coefficient
- latitudinal diffusion coefficient
- longitudinal advection coefficient
- latitudinal advection coefficient
- longitude of the origin of the allele
- latitude of the origin of the allele

for time periods before and after 5000 years. 95% Confidence intervals are produced as well.

It is possible to produce a GIF with the inferred allele dynamics.
This is done using the file stepplots.R and it takes as an input the output file produced by stepadna.R.
An example of a command for doing this:
```
Rscript stepplots.R -f result.csv -n 10
```

The input parameters are:
```
-f --file – file produced by stepadna.R
-o --out - output file name (defaults to outplot.gif)
-n --nimage - number of images to display in the GIF (default is 20)
```

Example of an output GIF showigng tempo-spatial allele frequency dynamics:
![](https://github.com/RasaMukti/stepadna/blob/main/outplot.gif)
