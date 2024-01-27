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
"""Create GEF output folder containing GEF input file.

Create uniquely named GEF output folder and uniquely named GEF input
file. User chooses which GEF input file options to use and their values
by entering them in the file. More options are available in GEF and can
be added to the list of 'GEF INPUT OPTIONS' (See the info section of the
function ``create_GEF_input_file_and_folder`()`` below).
"""
import os
#------------------------------------------------------- GEF INPUT OPTIONS ------------------------------------------------------#
global_opt      = 'global'
"""Option for using user defined parameter values in GEF.

'global' must be used. 'local' option uses hard-coded 'optimal'
parameter values that overrides the parameter values of the GEF function
'MyParameters'.
"""
MyParam_opt     = 'MyParameters'
"""Activates the GEF function 'MyParameters'."""
LMD_opt         = 'lmd'
"""Produces 'list-mode-output'-file. File contains event-by-event 
fission fragment data produced by GEF simulation.
"""
LMD_plus_opt         = 'lmd+'
"""Produces 'list-mode-output'-file. File contains event-by-event 
fission fragment data produced by GEF simulation plus neutron and
gamma spectra.
"""
projectile_opt  = 'EN'
"""For simulating neutron-induced fission. GEF can simulate 'spontaneous'
or 'neutron-induced' fission.
"""
#------------------------------------------------------- MODULE VARIABLES ------------------------------------------------------#
elements = {'92':'U','94':'Pu'}
"""Dictionary containing name and mass number of chemical elements used as
targets in simulation of fission.

This dictionary can be extended for investigations of other elements.
"""
#-------------------------------------------------------- MODULE METHODS -------------------------------------------------------#
def create_GEF_input_file_and_folder(pth_GEF_cwd,Z_target,A_compound,E_reaction,MC_runs):
        """Create folder and GEF input file for individual GEF simulation
        when multi-threading.
        
        To avoid ambiguity when multi-threading, each GEF output folder
        and GEF input file has a unique name (See the `See Also` section).

        Parameters
        ----------
        pth_GEF_cwd : `str`
                Local path to individual GEF output folder used during 
                multi-threading. (See the `See Also` section).
        Z_target : `int`
                Atomic number of target chemical element. 
                
                Needed to create file name in prescribed TALYS format of 
                fission fragment files and for GEF simulations.
        A_compound : `int`
                Mass number of compound target nuclei. I.e A + 1.
                
                Needed to create file name in prescribed TALYS format of  
                fission fragment files and for GEF simulations.
        E_reaction : `float`
                Kinetic energy in MeV of incident neutron in fission event.
                
                Needed to create file name in prescribed TALYS format of
                fission fragment files and for GEF simulations.
        MC_runs : `int`
                Given in the GEF input file as a multiplication factor
                of the default number of Monte Carlo simulations. 
                
                The GEF default value is 10^5 MC simulations. 

        Returns
        -------
        Function has no return value.

        Notes
        -----
        The choices for the GEF options in this file are the ones 
        necessary for simulations using McPUFF. More options are 
        available. (See GEF 'ReadMe' file in `See Also` section).

        See Also
        --------
        GEF 2023-V1.1 'ReadMe' file: <https://www.khschmidts-nuclear-web.eu/GEF_code/GEF-2023-1-1/Standalone/Readme.txt>
        ``McPUFF_program.create_GEF_workingdir_and_inputfile()``
        ``TALYS_Input.py``

        Examples
        --------
        Example generates no output. Creates uniquely named file 
        containing uniquely named GEF input file at path provided.

        >>> GEF_input.create_GEF_input_file_and_folder("/local_path/", 92,
                                                         236, 2.53e-8, 10)
        []
        """
        #-----------------------------------------------------------------------------------------------------------------------#
        if str(E_reaction).endswith('.0'):                                                      # GEF cannot handle zero decimals. It converts that float to int. Non-zero decimals work fine.
                E_reaction = int(E_reaction)
        with open(os.path.join(pth_GEF_cwd,'MyParameters.dat'),'w') as myparam_dat:			
                pass                                                                            # If no MyParameters.dat file, create. If exist, replace with empty
        GEF_in_file_name = f'{elements[str(Z_target)]}{A_compound-1}{projectile_opt}.in'        # Input file named after target element (I.e A-1).
        with open(os.path.join(pth_GEF_cwd,'file.in'),'w') as in_file:                          # If no file, create one. If exist, replace contents
                in_file.write(f'"in/{GEF_in_file_name}"\nEND')                                  # End file with EOF statement "END".
        GEF_in_folder = os.path.join(pth_GEF_cwd,'in',"")                                       # Must end with empty string "".
        if not os.path.exists(GEF_in_folder):                                                   # If directory doesn't exist, create it
                os.mkdir(GEF_in_folder)
        for file_name in os.scandir(GEF_in_folder):                                             # Remove previous files.
                if file_name.is_file():
                        os.remove(os.path.join(GEF_in_folder,file_name.name))
        iterations = int(MC_runs/1e5)                                                           # Converts number of GEF MC-simulations given by user to input file multiplication factor.
        with open(os.path.join(GEF_in_folder,GEF_in_file_name),'w') as GEF_input_file:          # Write GEF input data to file. LMD_plus_opt adds n and g energies. (Takes about 4 times longer).
                        GEF_input_file.write(f'{str(iterations)}\n' \
                                             f'{float(E_reaction)}\n'    \
                                             #f'Options({global_opt},{MyParam_opt},{LMD_opt})\n' \
                                             f'Options({global_opt},{MyParam_opt},{LMD_plus_opt})\n' \
                                             f'{Z_target}, {A_compound}, "{projectile_opt}"\n'  \
                                             f'END')                                                   # End file with EOF statement "END".
        del GEF_in_file_name,in_file,GEF_in_folder,GEF_input_file 
#------------------------------------------------------- END OF METHOD ---------------------------------------------------------#
                        
#------------------------------------------------------- END OF MODULE ---------------------------------------------------------#