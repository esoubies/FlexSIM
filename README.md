# FlexSIM

FlexSIM, for *flexible* SIM reconstruction, aims at providing reliable SIM reconstructions for a variety of SIM data, going from ``ideal'' ones acquired under standardized protocols and configurations, to ones obtained under more challenging settings and subject to reconstruction artifacts. 
More details can be found in the following paper:

[1] <a href=" ">Handling Challenging Structured Illumination Microscopy Data with FlexSIM.</a>, 
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
* Folder **Examples**:  Scripts to download and reconstruc the 20 open datasets described in [Table 1, 1]. See Section [Examples](#examples) below.

## FlexSIM parameters


## Examples


## Conditions of use

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/> .

Whenever you present or publish results that are based on this library, please cite [1].

