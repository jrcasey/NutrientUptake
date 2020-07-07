# NutrientUptake

Accompanying Matlab code to Casey, J.R., Follows, M.J., 2020. A steady-state model of microbial acclimation to substrate limitation. *PLOS Computational Biology* (in review)

![Figure 3](https://github.com/jrcasey/NutrientUptake/blob/master/assets/Figure_3.jpg)

An application of our steady-state acclimation model of nutrient transport: *Escherichia coli* K12 BW25113 growth in glucose-limited chemostats and glucose-replete batch cultures ([Schmidt *et al*., 2016](https://www.nature.com/articles/nbt.3418)) using the latest genome-scale metabolic model of *E. coli* K12 MG1566 (*i*ML1515; [Monk *et al*., 2017](https://www.nature.com/articles/nbt.3956)). The code will compute parameters in this figure (*e.g*., n<sup>\*</sup>, n<sup>max</sup>, S<sup>\*</sup>, S<sup>l\*</sup><sub>lb</sub>, S<sup>\*</sup><sub>ub</sub>, v, v<sub>max</sub><sup>G</sup>, v<sup>D</sup>, &mu;, &mu;<sub>max</sub>, *etc*).

No bugs found so far with Matlab 2019b. Parameters determined by FBA and molecular modeling are not included here (too many dependencies to wrangle) but I'd be happy to guide you through our approach to that if you're interested in applying this model to another microbe or another transport process. 

Please email me (jrcasey at mit dot edu) with feedback or questions!

## Instructions
1. Fork or clone this repo `https://github.com/jrcasey/NutrientUptake.git` to a local directory
2. Add the new directory to Matlab's path (you could use `pathtool` or `addpath(path/to/repo))`)
3. Navigate to `src/NutrientUptake.m` and hit run!

You should get a plot that looks like this:
![Figure 5](https://github.com/jrcasey/NutrientUptake/blob/master/assets/Figure_5.jpg)
