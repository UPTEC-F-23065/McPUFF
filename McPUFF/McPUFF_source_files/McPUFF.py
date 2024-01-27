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
"""Main script for McPUFF

This script runs the perturbed simulations using McPUFF.
All inputs to GEF and TALYS and all necessary local paths are declared
here. There are other parameters that control the simulations that the
user can adjust. These are described under the heading 
"Additional McPUFF simulation parameters" in the `Notes` section.

Notes
-----
There are a number of internal variables in McPUFF that controls various
aspects of the simulations, such as e.g. scaling numbers for parameters 
and ranges for distributions. It was decided against a GUI for setting
parameters and variables as McPUFF is intended to be run on computer 
clusters where it may not be possible to install the necessary modules. 
The variables that the user can adjust are specified in the list below. 
On the left are listed the classes and functions in "McPUFF_program" in
which the user can adjust the variables. On the right is a short 
desription of the variable. More information can be found in the 
respective classes and functions in ``McPUFF_program.py``.

The word "special" in the variable names on the right indicate GEF
parameters whose default values is zero. E.g ``st_dev_special_param``.

#------------------------------------------ ADDITIONAL McPUFF SIMULATION PARAMETERS --------------------------------------------#
Reaction.print_TALYS_ff_files()                           # Set E-value in printed in ".ff"-file name. (TALYS interpolates between files with different energies in the file name).
Modified_Parameter.__init__.                              # Set values for "max-min" distribution in "Single_Parameters" mode.
Modified_Parameter.__init__.                              # Set "scaling_number_uniform" for random numbers from uniform distribution for ordinary parameters in "Single_Parameters" mode.
Modified_Parameter.create_perturbed_parameter_value()     # Set "scaling_special_case_parameters" for random number from uniform distribution for special parameters in "Single_Parameters" and "TMC" modes.
Modified_Parameter.create_perturbed_parameter_value()     # Set "st_dev_special_param" for random numbers from normal distribution for special parameters in "Single_Parameters" and "TMC" mode.
TMC_Mod_Param_object.__init__.                            # Set "scaling_number_uniform" for random number from uniform distribution for ordinary parameters in "TMC" mode.
Reaction.__init__.                                        # Set "number_of_workers" for CPU's to use in "Single_Parameters" mode.
Reaction.__init__.                                        # Set "max_multithreads_TMC" for CPU's to use in "TMC" mode.
Reaction.delete_GEF_result_folder()                       # Turn off the deletion of GEF runtime data.  Be warned that the data amount can be quite large (GB).
Reaction.delete_TALYS_result_files()                      # Turn off the deletion of TALYS runtime data. Be warned that the data amount can be quite large (GB). 
Reaction.delete_TALYS_ff_files()                          # Turn off the deletion of ".ff" files in the GEF library during runtime. Be warned that the number of files equals twice the number of "TMC" simulations.
Reaction.perturbed_calculations_single_parameter()        # Set "simultaneous_threads_per_param" to assign number of CPU's to use for multi-threading. (Divide "number_of_workers" between multi-threads).
#-------------------------------------------------------------------------------------------------------------------------------#

See Also
--------
P. Karlsson, "Total Monte Carlo of the fission model in
GEF and its influence on the nuclear evaporation in TALYS," 
Technical report UPTEC F 23065, Uppsala University, 2023, 
url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>, 
accessed 1 Janauary 2024.

``McPUFF_program.Reaction``
``McPUFF_program.Modified_Parameter``
``McPUFF_program.Random_Parameter_value``
``McPUFF_program.TMC_Object``
``McPUFF_program.TMC_Mod_Param_object``

Examples
--------
The script instatiates a ``Reaction`` object in the script and produces a
`pickle` file at the specified location. The `pickle` file can be loaded
later with the class ``Perturbed_Parameter_TMC_Data``.

>>> reac = Reaction(92,236,2.53e-8,1e6,500,"normal","pth_GEF_program",
                                "pth_TALYS_program","pth_main",True,"TMC")
[reac (Reaction object), "Z92_A236_n_E2.53e-08MeV.pkl" ]
"""
from package_McPUFF import * 
from package_McPUFF.McPUFF_program import Reaction
#================================================= MODULE VARIABLES ============================================================#

#------------------------------------------------- PATH DECLARATIONS -----------------------------------------------------------#
pth_GEF_program = "/__local_path_to_modified_GEF"                               # No ´/´ at end of paths. 
"""Local path to modified GEF software."""               
pth_TALYS_program = "/__local_path_to_modified_TALYS"                           # No ´/´ at end of paths.
"""Local path to modified TALYS software."""
pth_main = "/__local_path_to_McPUFF_source_files"                               # No ´/´ at end of paths.
"""Local path to McPUFF main script "McPUFF.py" and "package_McPUFF.py"."""
#------------------------------------------------- SIMULATION INFORMATION ------------------------------------------------------#
Z_target = 92
"""Atomic number of target nucleus (`int`)."""
A_compound = 236
"""Atomic mass number of compound nucleus (`int`)."""
E_reaction = 2.53e-08 		
"""Kinetic energy of incident fission neutron (MeV) (`float`)."""
MC_runs = 1e5
"""Number of Monte Carlo simulation in modified GEF (`int`)."""
number_of_randoms = 3
"""Number of Total Monte Carlo (TMC) simulations in McPUFF (`int`).

One simulation is performed for every random number used to perturb
the GEF parameters. Hence, the number of random numbers equals the
number of "TMC" simulations."""
#--------------------------------------------- McPUFF SIMULATION MODE ----------------------------------------------------------#
#program_flag = 'Single_Parameters'             # Run McPUFF in "Single_Parameters" mode.
program_flag = 'TMC'                            # Run McPUFF in "TMC" mode.
"""Determines which McPUFF mode is used for "TMC" simulations."""
#---------------------------------- RUN McPUFF SIMULATIONS WITH OR WITHOUT TALYS SIMULATIONS -----------------------------------#
#TMC_with_TALYS = True                          # Run McPUFF with both GEF and TALYS simulations.
TMC_with_TALYS = False                          # Run McPUFF with GEF simulations only.
"""Determines whether TALYS simulations should be performed or not in McPUFF.

If variable is set to "False", McPUFF only performs simulation with GEF."""
#---------------------------------- SPECIFY WHICH DISTRIBUTION TO DRAW RANDOM NUMBERS FROM -------------------------------------#
#distribution_flag = 'uniform'                  # Draws random numbers from a uniform distribution.
distribution_flag = 'normal'                    # Draws random numbers from a normal (Gaussian) distribution.
#distribution_flag = 'max-min'                  # Performs McPUFF simulations using scaling provided by user.
"""Specifies from which distribution to draw random numbers used in simulation.

See the "See Also" section.
"""
#======================================================= RUN McPUFF SIMULATION =================================================#
reac = Reaction(Z_target,A_compound,E_reaction,MC_runs,number_of_randoms,distribution_flag,pth_GEF_program,pth_TALYS_program,pth_main,TMC_with_TALYS,program_flag)
#-------------------------------------------------------- END OF SCRIPT --------------------------------------------------------#