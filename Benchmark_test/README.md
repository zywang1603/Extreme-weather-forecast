# Large Sample Nowcasting Evaluation
This repository contains the Python scripts that were used for the large sample nowcasting evaluation in the publication of Imhoff et al. (2020). The repository contains three main folders: pysteps, rainymotion and HPCrunScript. The first two folders contain the model configurations of pysteps (Pulkkinen et al., 2019) and rainymotion (Ayzel et al., 2019) as used in this study. See also https://github.com/pySTEPS/pysteps and https://github.com/hydrogo/rainymotion for the latest versions of these nowcasting algorithms. Other than applied in the original distributions of these algorithms, the folders also contain a 'catchment_slice.py' script which saves the nowcasting outputs clipped to the given catchment extents, as based on provided shapefiles. This step is implemented to lower the data storage, since the focus of this study was on the use of nowcasting for hydrological basins and polder areas. The scripts can be found in ./pysteps/pysteps/utils/ and ./rainymotion/rainymotion/. 

The folder HPCrunScript contains the Python scripts that were used on the cluster to run the models. Scripts are available for the runs with: the pysteps deterministic setup, the pysteps probabilistic setup, the pysteps setup with only advection and all rainymotion models. Note that the mentioned events in the scripts for the pysteps runs are just a sample of the 1533 events.    

# Usage
For the use of rainymotion and pysteps, please see the aforementioned github pages. Both have example run scripts. 

# Terms of use
pysteps (Copyright (c) 2019, PySteps developers) is provided under the BSD 3-Clause. rainymotion (Copyright (c) 2018, Georgy Ayzel) is provided under an MIT license. Both algorithms are free to use. The pysteps and rainymotion folders in this repository follow the mentioned terms of use for pysteps and rainymotion. This repository is in no way a new means of distributing the source code of both algorithms. It is meant for the reproducibility of the study conducted by Imhoff et al. (2020). For the latest distribution and examples of both algorithms, please see the aforementioned github pages. 

The scripts in the 'HPCrunScripts' folder are free to use: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your opinion) any later version. These scripts are redistributed in the hope that they will be useful, but without any warranty and without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licences/. The producer of this tool is by no means to be held responsible for the use of and the results of this tool.

# References
Ayzel, G., Heistermann, M., & Winterrath, T. (2019). Optical flow models as an open benchmark for radar-based precipitation nowcasting (rainymotion v0.1). Geoscientic Model Development, 12 (4), 1387-1402. doi: 10.5194/gmd-12-1387-2019.

Imhoff, R. O., Brauer, C. C., Overeem, A., Weerts, A. H., & Uijlenhoet, R. (2020). Spatial and temporal evaluation of radar rainfall nowcasting techniques on 1,533 events. Water Resources Research, 56, e2019WR026723. doi: 10.1029/2019WR026723.

Pulkkinen, S., Nerini, D., Pérez Hortal, A. A., Velasco-Forero, C., Seed, A., Germann, U., & Foresti, L. (2019). Pysteps: an open-source Python library for probabilistic precipitation nowcasting (v1.0). Geoscientific Model Development, 12 (10), 4185-4219. doi: 10.5194/gmd-12-4185-2019. 
