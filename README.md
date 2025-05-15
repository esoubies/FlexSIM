# FlexSIM

FlexSIM, for *flexible* SIM reconstruction, aims at providing reliable SIM reconstructions for a variety of SIM data, going from ``ideal'' ones acquired under standardized protocols and configurations, to ones obtained under more challenging settings and subject to reconstruction artifacts.  
More details can be found in the following paper:

<a href="https://onlinelibrary.wiley.com/doi/10.1111/jmi.13344" target="_blank">Surpassing Light Inhomogeneities in Structured-Illumination Microscopy with FlexSIM.</a>,  <br />
Journal of Microscopy (2024), <br />
E. Soubies, A. Nogueron, F. Pelletier, T. Mangeat, C. Leterrier, M. Unser, and D. Sage.

## News

* May 2025: New saving options and integration of the newly proposed _<a href="https://www.nature.com/articles/s41592-025-02667-6" target="_blank">Dark-Sectioning</a>_ method.
* Feb 2025: New examples have been added (see below) 
* Dec 2024: The handling of temporal stacks has been improved. A new parameter *eqOrr* has been introduced to impose a soft constraint of equally spaced orientations during patterns estimation. (see parameters below).
* Nov 2024: FlexSIM is now compatible with GPU computation (with Matlab Parrallel Toolbox). Set params.GPU = 1 and it's done!
* March 2024: New parameters to ease the handling of temporal stacks.
* Feb 2024: Integration of a pure Matlab version of the VMLMB optimization method in GlobalBioIm which is used in FlexSIM. As such, no need anymore to compile mex files.

## Getting started 

Download or clone this repository and run the script **InstallFlexSIM.m** which will download the required <a href="https://biomedical-imaging-group.github.io/GlobalBioIm/">GlobalBioIm library</a> v1.2 (or more recent releases) and make neceesary changes in your Matlab path.

## Repository content

The repository is organized as follows.
* File **InstallFlexSIM.m**: FlexSIM installation script
* File **FlexSIM.m**: Main function of FlexSIM, the one that needs to be run with parameters as input (see provided examples)
* Folder **src**: Matlab source files of FlexSIM
* Folder **Examples**:  Scripts to download and reconstruct the 20 open datasets described in Table S1 of [1] (see [Examples](#examples) below).

## FlexSIM parameters

### General parameters

| Parameter | Description |
|------|------|
| DataPath | *Path to the SIM raw stack.* |
| pathToFlexSIM | *Path to FlexSIM root folder.* |
| displ | *From 0 to 2 with increasing number of display.* |
| verbose | *From 0 to 2 with increasing text displays.* |
| sav | *Saving options (0: no saving / 1: save reconstruction / 2: save patterns and reconstruction).* |
| dataType | *Data type for saving: 'uint8', 'uint16', 'uint32', 'single' (default)* |
| GPU | *Boolean on whether to use GPU or not* |
| parallelProcess | *Boolean on whether to use parallel computing. Requires the parallel computing toolbox.*  |

### Physical parameters and pre-processing

| Parameter | Description |
|------|------|
|  | ----  **Reconstruct on ROI / specific frames** | 
| SzRoi | *Size (px) of the ROI of the data considered for reconstruction* |
| posRoi | *Position of the top-left corner of the ROI.* |
| frameRange | *To treat only a subset of temporal frames of the stack (e.g., [1;5]).* |
|  | ---- **Patterns**|
| StackOrder | *Stack order: ap, pa, apw, paw, wap, wpa, axp, pxa (with a=angles, p=phases, w=widefield). The last two correspond to 9x9 montage with angles (resp phases) in row axp (resp pxa).* |
| nbOr | *Number of orientations.* |
| nbPh | *Number of phases.* |
|  | ---- **OTF Approximation** | 
| lamb | *Emission wavelength (nm).* |
| res | *Pixel size (nm).* |
| Na | *Objective numerical aperture.* |
| damp | *damping parameter in [0,1] (1= no damping) to attenuate middle freq in the approx OTF.* |
|  | ----  **Background estimation** | 
| SzRoiBack | *Size (px) of the ROI for background estimation (position automatically detected to minimize the intensity within the ROI).* |

### Patterns estimation

| Parameter | Description |
|------|------|
|  SzRoiPatt  |  *Size (px) of the ROI for pattern estimation.*  |
|  posRoiPatt | *Position of the top-left corner of the ROI for pattern estimation. (if empty automatically detected to maximize the intensity within the ROI)* |
|  maskWF |    *Radius (as a factor of the cutoff freq.) of the disk used to mask central Fourier frequencies.* |
|  ringRegionSearch |  *Lower and upper limits of Fourier ring region to search peaks (given as factor of the cutoff freq.).*  |
|  eqPh |  *Boolean, if true equally-spaced phases are assumed.*  |
|  eqOrr | *Boolean, if true equally-spaced orientation are assumed (constraint imposed in a soft way).* |
|  estiPattLowFreq | *Boolean, if true, estimate the low-freq. component of the patterns.*   |
|  doRefinement | *If false, do not performs the refinement step* |
|  pattAmp | *Amplitude a of the pattern (to be adjusted manually).*   | 
|  framePattEsti | *To use only a subset of frames for estimating a common pattern to all frames (empty to use all frames). Only used when cstTimePatt =1.* |
|  cstTimePatt | *If true, a single set of patterns is estimated and used for all frames.* |
|  rollMed | *Size of the rolling median interval used to adjust patterns parameters accros frames. Can only be used with cstTimePatt=0.* |    

### Image reconstruction

| Parameter | Description |
|------|------|
|   | ---- **OTF Attenuation**  |
| OTFAttStr  |  *Strength of the OTF attenuation (in [0,1]). If 0 no OTF attenuation.* |
| OTFAttwdth  | *Width of the OTF attenuation (>0). If 0 no OTF attenuation.* |
|   | ---- **Dark-Sectioning** (from this _<a href="https://www.nature.com/articles/s41592-025-02667-6" target="_blank">paper</a>_)  |
| DarkSec  |  *If >0 activate the Dark-sectioning (1 for middle and 2 for strong out-of-focus).* |
| DarkSecThres  | *Threshold used to define a mask on the image that contains the background.* |
|   | ---- **Cost function** |
|  apodize | *Boolean on whether to use apodization on boundaries.* |
|  sepOrr  | *Boolean on whether to treat each orientation separately.* |
|  padSz | *Padding size (px) used in the forward operator.* |
|  mu | *Regularization parameter.* |
|  regType | *Regularizer: 1 - 1st-order Tikhonov, 2 - Total Variation, 3 - Good roughness.*  |
|   | ---- **Optimization** |
| maxIt  | *Maximum number of iterations (stopping criteria).* |
| stepTol  | *elative error tolerance between two iterates (stopping criteria).* |




## Examples

This folder contains scripts to download and reconstruct the 25 open 2D-SIM datasets described in Table S1 of [1]. Each script is made of the following two steps
* download raw data and set FlexSIM parameters
* run FlexSIM.m function

The 20 open 2D-SIM datasets are sourced from 8 publications including **FairSIM** [1], **OpenSIM** [2], **HiFi-SIM** [3], **ML-SIM** [4], **JSFR-SIM** [5], **Direct-SIM** [6], **PCA-SIM** [7], **B-SIM** [8], as well as the dark-sectioning paper [9]. This corresponds to a collection of SIM data  acquired with a diversity of SIM systems and configurations.
Each subfolder of the Example folder corresponds to one dataset with the following naming convention
* _Reference_SIM-type_Bio-structure_

[1] _<a href="https://www.nature.com/articles/ncomms10980#citeas" target="_blank">FairSIM</a>_, M. Müller, V. Mönkemöller, S. Hennig, W. Hübner, and T. Huser, Open-source image reconstruction of super-resolution structured illumination microscopy data in ImageJ, Nat. Commun., vol. 7, no. 1, no. 1, Mar. 2016.

[2] _<a href="https://ieeexplore.ieee.org/document/7400963" target="_blank">OpenSIM</a>_, A. Lal, C. Shan, and P. Xi, Structured Illumination Microscopy Image Reconstruction Algorithm, IEEE Journal of Selected Topics in Quantum Electronics, vol. 22, no. 4, Jul. 2016.

[3] _<a href="https://www.nature.com/articles/s41377-021-00513-w" target="_blank">HiFi-SIM</a>_, G. Wen et al., High-fidelity structured illumination microscopy by point-spread-function engineering, Light Sci Appl, vol. 10, no. 1, no. 1, Apr. 2021.

[4] _<a href="https://opg.optica.org/boe/fulltext.cfm?uri=boe-12-5-2720&id=450173" target="_blank">ML-SIM</a>_, C. N. Christensen, E. N. Ward, P. Lio, and C. F. Kaminski, ML-SIM: A deep neural network for reconstruction of structured illumination microscopy images, Biomed. Opt. Express, vol. 12, no. 5, May 2021. 

[5] _<a href="https://www.spiedigitallibrary.org/journals/advanced-photonics/volume-4/issue-02/026003/High-speed-image-reconstruction-for-optically-sectioned-super-resolution-structured/10.1117/1.AP.4.2.026003.full?SSO=1" target="_blank">JSFR-SIM</a>_, Z. Wang et al., High-speed image reconstruction for optically sectioned, super-resolution structured illumination microscopy, Advanced Photonics, vol. 4, no. 2, Mar. 2022.

[6] _<a href="https://photonix.springeropen.com/articles/10.1186/s43074-023-00092-6" target="_blank">Direct-SIM</a>_, G. Wen et al., Spectrum-optimized direct image reconstruction of super-resolution structured illumination microscopy, PhotoniX, vol. 4, no. 1, June 2023.

[7] _<a href="https://elight.springeropen.com/articles/10.1186/s43593-022-00035-x" target="_blank">PCA-SIM</a>_, J. Qian, et al., Structured illumination microscopy based on principal component analysis, eLight, vol. 3, no. 1, Feb. 2023.

[8] _<a href="https://www.nature.com/articles/s44303-024-00066-8" target="_blank">B-SIM</a>_,  A. Saurabh, et al., Approaching maximum resolution in structured illumination microscopy via accurate noise modeling, npj Imaging, vol. 3, no. 5, Jan. 2025.

[9] _<a href="https://www.nature.com/articles/s41592-025-02667-6" target="_blank">Dark-SIM</a>_,  R. Cao et al., Dark-based optical sectioning assists background removal in fluorescence microscopy, Nature Methods, 2025.

## Conditions of use

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/> .

Whenever you present or publish results that are based on this repository, please cite:

E. Soubies, A. Nogueron, F. Pelletier, T. Mangeat, C. Leterrier, M. Unser, and D. Sage. Surpassing Light Inhomogeneities in Structured-Illumination Microscopy with FlexSIM. Journal of Microscopy, 2024.

