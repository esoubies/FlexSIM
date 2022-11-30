# Membrane - ML-SIM

[1] E. N. Ward et al., “MAI-SIM: interferometric multicolor structured illumination microscopy for everybody.” arXiv, Jun. 01, 2022. doi: 10.48550/arXiv.2206.00532.

### Overview

Presented originally in the MAI-SIM publication, this dataset shows a cell membrane imaged with MAI illumination and standard 9-image SIM. It can be downloaded from [here](https://github.com/charlesnchr/ML-SIM/raw/master/Test_data/AtheiSIM%2014-52_488.tif). 

### Imaging Parameters

<img src="https://latex.codecogs.com/gif.latex?\lambda: " /> 530

Pixel size: 64 nm

_NA_: 1.2

### Patter Parameters Estimation

We reconstructed this data set with _flexSIM_, _fairSIM_ and _HiFi-SIM_, estimating with the pattern estimation parameters described below.

| Method  | K_x  [px] | K_y [px] | Phase [rad] | OTF Att Strength (FWHM) | OTF Comp | 
| ------------- |:-------------:| ------------- | ------------- | | |
| flexSIM (or 1) | 114.607  | 76.973 | 0.555 | 0.99 (0.1) | 0.3 |
| flexSIM (or 2) | 125.030 | -70.161 | 1.200 | | |
| flexSIM (or 3) | -1.559 | -136.092 | 2.163 | | |
| fairSIM (or 1) | 114.5 | 76.811 | -0.461 | - | 0.15 |
| fairSIM (or 2) | 125.078 | -70.166 | 1.295 |  | |
| fairSIM (or 3) | 1.656 | -136.078 | 0.085 | |
| HiFi (or 1) | 114.525 | 59.062 | -0.4373 | 0.97 (0.48) | 0.9| 
| HiFi (or 2) | 125.092 | -70.154 | 1.1858 |  | |
| HiFi (or 3) | 1.660 |  -1.361 | 0.1106 |  | | |