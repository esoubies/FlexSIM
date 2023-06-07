# Microtubules - ML-SIM

[1] E. N. Ward et al., “MAI-SIM: interferometric multicolor structured illumination microscopy for everybody.” arXiv, Jun. 01, 2022. doi: 10.48550/arXiv.2206.00532.

### Overview

Presented originally in the MAI-SIM publication, this dataset shows COS7 Microtubules imaged with standard 9-image SIM and MAI illumination. The dataset can be downloaded [here](https://github.com/charlesnchr/ML-SIM/raw/master/Test_data/AtheiSIM%2013-43_488.tif)

### Imaging Parameters

<img src="https://latex.codecogs.com/gif.latex?\lambda: " /> 530

Pixel size: 64 nm

_NA_: 1.2

### Patter Parameters Estimation

We reconstructed this data set with _flexSIM_, _fairSIM_ and _HiFi-SIM_, estimating with the pattern estimation parameters described below.

| Method  | K_x  [px] | K_y [px] | Phase [rad] | OTF Att Strength (FWHM) | OTF Comp | 
| ------------- |:-------------:| ------------- | ------------- | | |
| flexSIM (or 1) | 135.869  | 2.789 | 0.226 | 0.99 (0.1) | 0.3 |
| flexSIM (or 2) | 70.174 | 124.696 | 2.412 | | |
| flexSIM (or 3) | -77.260 | 114.963 | 1.224 | | |
| fairSIM (or 1) | 136.078 | 2.944 | -0.995 | - | 0.15 |
| fairSIM (or 2) | 70.433 | 124.766 | -2.606 |  | |
| fairSIM (or 3) | -77.278 | 114.833 | -2.33 | |
| HiFi (or 1) | 136.031 | 2.969 | -0.9202 | 0.97 (0.48) | 0.9| 
| HiFi (or 2) | 70.401 | 124.784 | 2.5428 |  | |
| HiFi (or 3) | 77.278 | -114.846 | -2.3761 |  | | |