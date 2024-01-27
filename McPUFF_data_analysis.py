#-*- coding utf-8 -*-
#
# This file is part of McPUFF.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""Script for running data analysis with functions in 
``McPUFF_Perturbed_Data`` using information in McPUFF `pickle` files.

The function ``McPUFF_Perturbed_Data.Correlation_division_Eexc_light_heavy_Single_Parameters()``
function contains an example of how to load data from a `pickle` file 
with McPUFF data from a ``Single_Parameters`` simulation. Data from a 
``TMC`` simulation can be loaded in the same manner, but the data will be
stored in other lists of the ``Reaction`` object. An example for the 
loading of data from a ``TMC`` simulation is shown in the function
``McPUFF_Perturbed_Data.Correlation_division_Eexc_light_heavy_TMC()``.

See the 
'McPUFF_object_structure' in the `See Also` section for more information.

See Also
--------
``McPUFF_Perturbed_Data``
McPUFF_object_structure: url: <https://github.com/UPTEC-F-23065/McPUFF/blob/ac8c29ce13729ff85d598f8165555891a9e56ea4/ObjectStructureMcPUFF.png>
````
"""
from package_McPUFF.McPUFF_Perturbed_Data import McPUFF_Perturbed_Data as MPD
pth_PKL = "/__local_path_to_pickle_file(s)__/"                                                  # Add local path to pickle file(s).
MPD_object = MPD(pth_PKL) 
#===============================================================================================================================#
#-------------------- EXAMPLE FOR LOADING PICKLE DATA FROM A SINGLE_PARMETERS MODE SIMULATION IN McPUFF ------------------------#
# MPD.Correlation_division_Eexc_light_heavy_Single_Parameters(MPD_object)
#-------------------- EXAMPLE FOR LOADING PICKLE DATA FROM A TMC MODE SIMULATION IN McPUFF -------------------------------------#
# MPD.Correlation_division_Eexc_light_heavy_TMC(MPD_object)
#===============================================================================================================================#

#------------------------------------------------------- END OF MODULE ---------------------------------------------------------#