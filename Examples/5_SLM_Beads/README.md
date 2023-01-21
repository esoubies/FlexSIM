# Beads - ML-SIM

[1] C. N. Christensen, E. N. Ward, P. Lio, and C. F. Kaminski, “ML-SIM: A deep neural network for reconstruction of structured illumination microscopy images,” arXiv, arXiv:2003.11064, Mar. 2020. doi: 10.48550/arXiv.2003.11064.

### Overview

Presented originally in the ML-SIM publication [1], this dataset shows a standard 9-image SIM acquisition of bead structures imaged with SLM illumination. The dataset can be downloaded from [here](https://github.com/charlesnchr/ML-SIM/raw/master/Test_data/SLM-SIM%20beads.tif)

### Imaging Parameters

<img src="https://latex.codecogs.com/gif.latex?\lambda: " /> 550

Pixel size: 64 nm

_NA_: 1.2

### Patter Parameters Estimation

We reconstructed this data set with _flexSIM_, _fairSIM_ and _HiFi-SIM_, with the pattern estimation parameters described below.

| Method  | K_x  [px] | K_y [px] | Phase [rad] | OTF Att Strength (FWHM) | OTF Comp | 
| ------------- |:-------------:| ------------- | ------------- | | |
| flexSIM (or 1) | 0.151 | -146.157 | 1.708 | - (-) | 0.3 |
| flexSIM (or 2) | -126.940 | 72.617 | 0.662 | | |
| flexSIM (or 3) | 126.832 | 73.561 | 1.969 | | |
| fairSIM (or 1) | 0.144 | -146.211 | -2.621 | - | 0.15 |
| fairSIM (or 2) | -126.967 | 72.633 | 1.339 |  | |
| fairSIM (or 3) | 126.811 | 73.544 | -2.220 | |
| HiFi (or 1) | 0.154 | -146.216 | -2.626 | 0.97 (0.48) | 0.9| 
| HiFi (or 2) | 126.969 | -72.661 | 1.180 |  | |
| HiFi (or 3) | 126.846 | 73.525 | -2.211 |  | | |