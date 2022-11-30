# Tetraspeck - fair-SIM

[1] M. Müller, V. Mönkemöller, S. Hennig, W. Hübner, and T. Huser, “Open-source image reconstruction of super-resolution structured illumination microscopy data in ImageJ,” Nat Commun, vol. 7, no. 1, Art. no. 1, Mar. 2016, doi: 10.1038/ncomms10980

### Overview

Presented in the fairSIM original publication, this dataset shows tetraspeck structures imaged with SLM-SIM. For a more homogeneus resolution improvement, 12 images (4 orientations with 3 phases each) were acquired. The dataset can be dowloaded from [here](https://github.com/fairSIM/test-datasets/releases/download/SLM-SIM-Bielefeld/SLM-SIM_Tetraspeck200_680nm.tif).

### Imaging Parameters

<img src="https://latex.codecogs.com/gif.latex?\lambda: " /> 525

Pixel size: 97.6 nm

_NA_: 1.49

### Patter Parameters Estimation

We reconstructed this data set with _flexSIM_ and _fairSIM_, estimating with the pattern estimation parameters described below.

| Method  | K_x  [px] | K_y [px] | Phase [rad] | OTF Att Strength (FWHM) | OTF Comp | 
| ------------- |:-------------:| ------------- | ------------- | | |
| flexSIM (or 1) | 2.123 | -91.227 | 1.632 | 0.99 (0.1) | 0.3 |
| flexSIM (or 2) | 92.066 | 1.836 | 0.556 | | |
| flexSIM (or 3) | 61.982 | 64.372 | 0.720 | | |
| flexSIM (or 4) | 65.046| -62.165 | 0.456 | | |
| fairSIM (or 1) | 2.1 | -91.211 | 0.310 | - | 0.15 |
| fairSIM (or 2) | 92.078 | 1.833 | 1.423 |  | |
| fairSIM (or 3) | 62.033 | 64.389 | 0.060 | |
| fairSIM (or 4) | 65.122 | -62.100 | 0.954 | |