# Endoplasmic - ML-SIM

[1] C. N. Christensen, E. N. Ward, P. Lio, and C. F. Kaminski, “ML-SIM: A deep neural network for reconstruction of structured illumination microscopy images,” arXiv, arXiv:2003.11064, Mar. 2020. doi: 10.48550/arXiv.2003.11064.

### Overview

Presented originally in the ML-SIM original publication, this dataset shows a standard 9-image SIM acquisition of endoplasmic structures imaged with SLM illumination. The dataset can be downloaded from [here](https://github.com/charlesnchr/ML-SIM/raw/master/Test_data/SLM-SIM%20syrlyso640-5.tif)

### Imaging Parameters

<img src="https://latex.codecogs.com/gif.latex?\lambda: " /> 510

Pixel size: 64 nm

_NA_: 1.2

### Patter Parameters Estimation

We reconstructed this data set with _flexSIM_, _fairSIM_ and _HiFi-SIM_, with the pattern estimation parameters described below.

| Method  | K_x  [px] | K_y [px] | Phase [rad] | OTF Att Strength (FWHM) | OTF Comp | 
| ------------- |:-------------:| ------------- | ------------- | | |
| flexSIM (or 1) | 114.146 | 0.268 | 1.192 | 0.99 (0.5) | 0.3 |
| flexSIM (or 2) | -56.656 | -99.832 | 3.0348 | | |
| flexSIM (or 3) | -57.626 | 99.479 | 0.399 | | |
| fairSIM (or 1) | 114.144 | 0.256 | 2.549 | - | 0.15 |
| fairSIM (or 2) | -56.678 | -99.722 | 2.831 |  | |
| fairSIM (or 3) | -57.589 | 99.367 | -1.971 | |
| HiFi (or 1) | 114.154 | 0.278 | 2.410 | 0.998 (0.37) | - | 
| HiFi (or 2) | 56.661 | 99.722 | 2.785 |  | |
| HiFi (or 3) | 57.599 | -99.401 | -2.101 |  | | |