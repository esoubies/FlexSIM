# COS7 Microtubules - HiFi-SIM

[1] _HiFi-SIM_, G. Wen et al., “High-fidelity structured illumination microscopy by point-spread-function engineering,” Light Sci Appl, vol. 10, no. 1, Art. no. 1, Apr. 2021, doi: 10.1038/s41377-021-00513-w.

### Overview

This dataset was presented first by Wen et al in [1] to test the _HiFi-SIM_. It was acquired with SLM illumination and a standard 9-image SIM. It consists of Microtubule structures, and can be downloaded [here](https://static-content.springer.com/esm/art%3A10.1038%2Fs41377-021-00513-w/MediaObjects/41377_2021_513_MOESM2_ESM.rar).

### Imaging Parameters

<img src="https://latex.codecogs.com/gif.latex?\lambda: " /> 525

Pixel size: 65 nm

_NA_: 1.49

### Patter Parameters Estimation

We reconstructed this data set with _flexSIM_, _fairSIM_ and _HiFi-SIM_, estimating with the pattern estimation parameters described below.

| Method  | K_x  [px] | K_y [px] | Phase [rad] | OTF Att Strength (FWHM) | OTF Comp | 
| ------------- |:-------------:| ------------- | ------------- | | |
| flexSIM (or 1) | 49.009 | 125.022 | 0.926 | 0.999 (0.3) | 0.2 |
| flexSIM (or 2) | 132.637 | 20.149 | 2.036 | | |
| flexSIM (or 3) | 83.571 | -104.859 | 3.072  | | |
| fairSIM (or 1) | 150.011 | -144.833 | -2.140 | 0.97 (1.5) | 0.15 |
| fairSIM (or 2) | -49.344 | -202.256 | 2.203 |  | |
| fairSIM (or 3) | -200.211 | -57.788 | 0.423 | |
| HiFi-SIM (or 1) | 150.031 | -144.846 | -2.110 | 0.97 (1.5) | |
| HiFi-SIM (or 2) | 49.340 | 202.216 | 2.091 |  | |
| HiFi-SIM (or 3) | 200.216 | 57.784 | 0.401 | |