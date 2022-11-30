# Microtubules - OpenSIM

[1] A. Lal, C. Shan, and P. Xi, “Structured Illumination Microscopy Image Reconstruction Algorithm,” IEEE Journal of Selected Topics in Quantum Electronics, vol. 22, no. 4, pp. 50–63, Jul. 2016, doi: 10.1109/JSTQE.2016.2521542.


### Overview

Presented in the OpenSIM original publication, this dataset shows Microtubules structures imaged with standard 9-image SIM. The dataset can be downloaded [here](https://github.com/LanMai/OpenSIM/raw/master/SIMexpt/sim01z4.tif)

### Imaging Parameters

<img src="https://latex.codecogs.com/gif.latex?\lambda: " /> 515

Pixel size: 60 nm

_NA_: 1.49

### Patter Parameters Estimation

We reconstructed this data set with _flexSIM_, _fairSIM_ and _HiFi-SIM_, estimating with the pattern estimation parameters described below.

| Method  | K_x  [px] | K_y [px] | Phase [rad] | OTF Att Strength (FWHM) | OTF Comp | 
| ------------- |:-------------:| ------------- | ------------- | | |
| flexSIM (or 1) | 120.647 | 33.586 | 1.124 | 0.99 (0.1) | 0.3 |
| flexSIM (or 2) | 89.255 | -87.630 | 1.955 | | |
| flexSIM (or 3) | 31.421 | 121.177 | 2.485 | | |
| fairSIM (or 1) | 120.744 | 33.522 | 2.287 | - | 0.15 |
| fairSIM (or 2) | 89.3 | -87.7 | -2.241 |  | |
| fairSIM (or 3) | 31.411 | 121.233 | -1.601 | |
| HiFi (or 1) | 120.722 | 33.524 | 2.332 | 0 (-) | - | 
| HiFi (or 2) | 89.278 | -87.722 | -2.091 |  | |
| HiFi (or 3) | 31.401 | 121.216 | -1.5144 |  | | |