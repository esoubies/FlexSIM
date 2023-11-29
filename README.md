# FlexSIM

FlexSIM, for *flexible* SIM reconstruction, aims at providing reliable SIM reconstructions for a variety of SIM data, going from ``ideal'' ones acquired under standardized protocols and configurations, to ones obtained under more challenging settings and subject to reconstruction artifacts. 
More details can be found in the following paper:

[1] <a href=" " target="_blank">Handling Challenging Structured Illumination Microscopy Data with FlexSIM.</a>, 
Preprint (2023), <br />
E. Soubies, A. Nogueron, F. Pelletier, T. Mangeat, C. Leterrier, M. Unser, and D. Sage.

## Getting started 

Download or clone this repository and run the script **InstallFlexSIM.m** which will download the required <a href="https://biomedical-imaging-group.github.io/GlobalBioIm/">GlobalBioIm library</a> v1.2 (or more recent releases) and make neceesary changes in your Matlab path.

## Repository content

The repository is organized as follows.
* File **InstallFlexSIM.m**: FlexSIM installation script
* File **FlexSIM.m**: Main function of FlexSIM, the one that needs to be run with parameters as input (see provided examples)
* Folder **src**: Matlab source files of FlexSIM
* Folder **Simulation**: Matlab scripts for patterns estimation experiments in [1]
* Folder **Examples**:  Scripts to download and reconstruct the 20 open datasets described in Table S1 of [1] (see [Examples](#examples) below).

## FlexSIM parameters

### General parameters

| Parameter | Description |
|------|------|
| DataPath | *Path to the SIM raw stack.* |
| pathToFlexSIM | *Path to FlexSIM root folder.* |
| displ | *From 0 to 2 with increasing number of display.* |
| verbose | *From 0 to 2 with increasing text displays.* |
| sav | *Boolean on whether to save the reconstructed image and estimated patterns.* |
| parallelProcess | *Boolean on whether to use parallel computing. Requires the parallel computing toolbox.*  |




## Examples

This folder contains scripts to download and reconstruct the 20 open 2D-SIM datasets described in Table S1 of [1]. Each script is made of the following two steps
* download raw data and set FlexSIM parameters
* run FlexSIM.m function

The 20 open 2D-SIM datasets are sourced from 7 publications including **FairSIM** [2], **OpenSIM** [3], **HiFi-SIM** [4], **ML-SIM** [5], **JSFR-SIM** [6], **Direct-SIM** [7], and **PCA-SIM** [8]. This corresponds to a collection of SIM data  acquired with a diversity of SIM systems and configurations.
Each subfolder of the Example folder corresponds to one dataset with the following naming convention
* _Reference_SIM-type_Bio-structure_

[2] _<a href="https://www.nature.com/articles/ncomms10980#citeas" target="_blank">FairSIM</a>_, M. Müller, V. Mönkemöller, S. Hennig, W. Hübner, and T. Huser, Open-source image reconstruction of super-resolution structured illumination microscopy data in ImageJ, Nat. Commun., vol. 7, no. 1, no. 1, Mar. 2016.

[3] _<a href="https://ieeexplore.ieee.org/document/7400963" target="_blank">OpenSIM</a>_, A. Lal, C. Shan, and P. Xi, Structured Illumination Microscopy Image Reconstruction Algorithm, IEEE Journal of Selected Topics in Quantum Electronics, vol. 22, no. 4, Jul. 2016.

[4] _<a href="https://www.nature.com/articles/s41377-021-00513-w" target="_blank">HiFi-SIM</a>_, G. Wen et al., High-fidelity structured illumination microscopy by point-spread-function engineering, Light Sci Appl, vol. 10, no. 1, no. 1, Apr. 2021.

[5] _<a href="https://opg.optica.org/boe/fulltext.cfm?uri=boe-12-5-2720&id=450173" target="_blank">ML-SIM</a>_, C. N. Christensen, E. N. Ward, P. Lio, and C. F. Kaminski, ML-SIM: A deep neural network for reconstruction of structured illumination microscopy images, Biomed. Opt. Express, vol. 12, no. 5, May 2021. 

[6] _<a href="https://www.spiedigitallibrary.org/journals/advanced-photonics/volume-4/issue-02/026003/High-speed-image-reconstruction-for-optically-sectioned-super-resolution-structured/10.1117/1.AP.4.2.026003.full?SSO=1" target="_blank">JSFR-SIM</a>_, Z. Wang et al., High-speed image reconstruction for optically sectioned, super-resolution structured illumination microscopy, Advanced Photonics, vol. 4, no. 2, Mar. 2022.

[7] _<a href="https://photonix.springeropen.com/articles/10.1186/s43074-023-00092-6" target="_blank">Direct-SIM</a>_, G. Wen et al., Spectrum-optimized direct image reconstruction of super-resolution structured illumination microscopy, PhotoniX, vol. 4, no. 1, June 2023.

[8] _<a href="https://elight.springeropen.com/articles/10.1186/s43593-022-00035-x" target="_blank">Direct-SIM</a>_, J. Qian, et al., Structured illumination microscopy based on principal component analysis, eLight, vol. 3, no. 1, Feb. 2023.


## Conditions of use

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/> .

Whenever you present or publish results that are based on this library, please cite [1].

