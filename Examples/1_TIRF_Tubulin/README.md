# Tubulin - fair-SIM

[1] M. Müller, V. Mönkemöller, S. Hennig, W. Hübner, and T. Huser, “Open-source image reconstruction of super-resolution structured illumination microscopy data in ImageJ,” Nat Commun, vol. 7, no. 1, Art. no. 1, Mar. 2016, doi: 10.1038/ncomms10980
[2] P. Kner, B. B. Chhun, E. R. Griffis, L. Winoto, and M. G. L. Gustafsson, “Super-resolution video microscopy of live cells by structured illumination,” Nat Methods, vol. 6, no. 5, Art. no. 5, May 2009, doi: 10.1038/nmeth.1324.

### Overview

This dataset was presented first by Gustafsson in [2] and used as an example in the fairSIM original publication [1]. It was acquired with TIRF illumination, presenting a challenging pattern estimation. It can be downloaded from [here](https://github.com/fairSIM/test-datasets/releases/download/TIRF-SIM-Georgia/TIRF_Tubulin_525nm.tif).

### Imaging Parameters

<img src="https://latex.codecogs.com/gif.latex?\lambda: " /> 525

Pixel size: 63 nm

_NA_: 1.49

### Patter Parameters Estimation

We reconstructed this data set with _flexSIM_, _fairSIM_ and _HiFi-SIM_, estimating with the pattern estimation parameters described below.

| Method  | K_x  [px] | K_y [px] | Phase [rad] | OTF Att Strength (FWHM) | OTF Comp | 
| ------------- |:-------------:| ------------- | ------------- | | |
| flexSIM (or 1) | 1.200 | 89.800 | 0.495 | 0.999 (0.3) | 0.3 |
| flexSIM (or 2) | 79.105 | 43.725 | 0.965 | | |
| flexSIM (or 3) | 78.010 | -45.822 | 2.435 | | |
| fairSIM (or 1) | 1.167 | 89.767 | 1.218 | - | 0.15 |
| fairSIM (or 2) | 79.100 | 43.700 | 2.088 |  | |
| fairSIM (or 3) | 77.989 | -45.922 | -1.0125 | |
| HiFi-SIM (or 1) | 1.216 | 89.722 | 1.198 | - | 0.05 |
| HiFi-SIM (or 2) | 79.093 | 43.660 | 2.263 |  | |
| HiFi-SIM (or 3) | 77.969 | -45.907 | -1.084 | |