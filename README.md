# NutrientUptake
![license badge](https://img.shields.io/github/license/jrcasey/NutrientUptake) ![languange badge](https://img.shields.io/github/languages/top/jrcasey/NutrientUptake) ![size badge](https://img.shields.io/github/repo-size/jrcasey/NutrientUptake)

Accompanying Matlab code to Casey, J.R., Follows, M.J., 2020. A steady-state model of microbial acclimation to substrate limitation. *PLOS Computational Biology* (accepted; link to ms soon)

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


### Captions to figures
Top figure:
*Schematic of steady-state nutrient acclimation. A) – Optimal transporter abundances $n^{\*}$ (solid green line) lie at the porter limitation (blue shaded region) boundary. All points above and below this boundary fail to maximize growth rate; for all points below, the uptake rate is sub-maximal; for all points above, some transporters are unoccupied. Over an interval of bulk substrate concentrations $S^{\infty} \leq S^{\*}$, this boundary is met with diffusion limitation (orange shaded region), where the catalytic rate exceeds the encounter rate; as concentrations exceed $S^{\*}$, the porter limitation boundary transitions to an internal growth rate limit (yellow shaded region), where the catalytic rate exceeds the rate of some downstream reaction. B) – For smaller cells, or for transporters with slower kcat, membrane surface area limitation (a special case of porter limitation) may be encountered. In the interval bounded by $S_{SA}^{lb}$ and $S_{SA}^{ub}$, $n^{\*}$ is constrained to $n_{max}$ (black dashed line). C) – Growth rates of optimally acclimated cells (solid green line) follow the intersection of two limits: the diffusive limit (purple dashed line) and the internal growth limit $\mu_{max}^{G}$ (black dashed line). These rate limits are, again, bisected by $S^{\*}$ with a sharp transition. D) – If surface area limitation is encountered, growth rates in the interval $S_{SA}^{lb} \leq S^{\infty} \leq S_{SA}^{ub}$ follow a more gradual, hyperbolic transition from diffusion limitation to growth limitation.*

Bottom figure:
*Model predictions and observations of \textit{Escherichia coli} acclimation to growth on glucose. Top panel -- Modeled and experimental abundances of the glucose transporter PtsG across steady-state concentrations of glucose (log scale). Contours indicate the corresponding uptake rates. The vertical green dashed line represents the critical substrate concentration $S^{\*}$. $n_{max}$ is above the plotted range. \textit{Bottom panel} -- Modeled and experimental uptake rates over an interval of glucose concentrations. The uptake rate profile for a glucose-replete batch acclimated culture is shown in cyan. Uptake rates for all acclimated phenotypes over the concentration range are shown in red, with model predictions shown as a continuous line and experimentally derived values from published data (Schmidt et al., 2016) shown as markers. For guidance, the maximum diffusive flux is shown as a dashed blue line, $S^{\*}$ is shown as a vertical dashed green line, and $v_{max}^{G}$ is indicated by the horizontal dashed orange line.*