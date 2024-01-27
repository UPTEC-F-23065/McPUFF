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
"""Create TALYS input file, assign values to chosen keywords and
place file in TALYS output folder. 

User can add or remove TALYS keywords to be used in the TALYS input file.
The values for the keyword are passed as arguments to the function 
``write_TALYS_input_file``, or simply entered as values in the file. 
"""
import os
def create_TALYS_input_file(pth_TALYS_input,unique_thread,TALYS_input_file_name,Z_target,A_compound,E_reaction):
        """Create TALYS input file.
        
        Creates uniquely named folder for individual TALYS run at location
        chosen by user. The ``write_TALYS_input_file`` function is called 
        to create a TALYS input file in prescribed format which is placed
        in the newly created TALYS output folder. 

        Parameters
        ----------
        pth_TALYS_input : `str`
                Path to location where the TALYS output folder is to be
                created. 
        unique_thread : `str`
                Unique ID for individual TALYS simulation during multi-
                thread simulations. 
                
                Given as value for the especially created TALYS keyword 
                `geffissionfileid`. See `Notes for more details.`
        TALYS_input_file_name : `str`
                Unique name for TALYS input file. 
                
                Created by function in `McPUFF_progam.py`. See `Notes` 
                for more details.
        Z_target : `int`
                Atomic number of target chemical element. 
                
                Needed to create file name in prescribed TALYS format of 
                fission fragment files and for TALYS simulations.
        A_compound : `int`
                Mass number of compound target nuclei. I.e A + 1. 
                
                Needed to create file name in prescribed TALYS format of 
                fission fragment files and for TALYS simulations.
        E_reaction : `float`
                Kinetic energy in MeV of incident neutron in fission event. 
                
                Needed to create file name in prescribed TALYS format of
                fission fragment files and for TALYS simulations.

        Returns
        -------
        Function has no return value.
        
        See Also
        --------
        ``McPUFF_program.py.create_perturbed_FY()``
        ``McPUFF_program.py.create_perturbed_TMC_FY()``

        Notes
        -----
        The ``McPUFF_program.py`` script performs multi-thread 
        simulations. To avoid ambiguity at run time, each TALYS simulation 
        has its own output folder containing its own input file. These 
        have unique names that are dependent on which mode is used. 
        For ``Single_Parameters`` mode the unique name is made up of a
        combination of parameter name and the enumeration of the 
        simulation. For ``TMC`` mode, the unique name is made up of string 
        ``TMC`` plus enumeration number of simulation. These name are 
        created in ``McPUFF_program.py`` (See the `See Also` section). The 
        unique name is passed to TALYS as the value for the keyword 
        ``geffissionfileid`` [1]_.

        References
        ----------
        .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
        GEF and its influence on the nuclear evaporation in TALYS." 
        Technical report UPTEC F 23065, Uppsala University, 2023, accessed
        1 Janauary 2024, url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>.
        
        Examples
        --------
        Examples generate no output. Creates uniquely named file 
        containing uniquely named TALYS input file at path provided.

        Example for `TMC` simulation (enumeration 5):
        >>> TALYS_Input.create_TALYS_input_file("/.../local_path/.../", 
                             "TMC_5", "TMC_5_input.in" ,92 , 236, 2.53e-8)
        []

        Example for ``Single_Parameters`` simulation (enumeration 5):
        >>> TALYS_Input.create_TALYS_input_file("/.../local_path/.../", 
         "_P_Z_Curve_S2_5", "_P_Z_Curve_S2_5_input.in" ,92 , 236, 2.53e-8)
        [] 
        """
        #-----------------------------------------------------------------------------------------------------------------------#
        TALYS_folder = os.path.join(pth_TALYS_input+f'{unique_thread}',"") # Must end with empty string "".
        if not os.path.exists(TALYS_folder):    # If directory doesn't exist, create it
                os.mkdir(TALYS_folder)
        write_TALYS_input_file(TALYS_folder,TALYS_input_file_name,unique_thread,Z_target,A_compound,E_reaction)         
        del TALYS_folder,TALYS_input_file_name 
        #------------------------------------------------------- END OF METHOD ---------------------------------------------------------#

def write_TALYS_input_file(folder, file_name,unique_thread_id,Z_target,A_compound,E_reaction):
        """Assign values to TALYS keywords and create TALYS input file.

        To exclude a keyword from the simulation, simply comment it away. 
        E.g:  # f'element {element}\\n'
        
        Parameters
        ----------
        folder : `str`
                Path to TALYS output folder. 
        file_name : `str`
                Unique TALYS input file name. 
                
                See ``TALYS_Input.create_TALYS_input_file()``. 
                .
        unique_thread_id : `str`
                Value for keyword `geffissionfileid`. See `Notes`.
        Z_target : `int`
                Atomic number of target chemical element. 
                
                Needed to create file name in prescribed TALYS format of 
                fission fragment files and for TALYS simulations.
        A_compound : `int`
                Mass number of compound target nuclei. I.e A + 1. 
                
                Needed to create file name in prescribed TALYS format of 
                fission fragment files and for TALYS simulations.
        E_reaction : `float`
                Kinetic energy in MeV of incident neutron in fission event. 
                
                Needed to create file name in prescribed TALYS format of 
                fission fragment files and for TALYS simulations.

        Returns
        -------
        Function has no return value.

        See Also
        --------
        `Modification-of-TALYS-for-TMC-simulations.git 
        <https://github.com/UPTEC-F-23065/Modification-of-TALYS-for-TMC-simulations.git>`_.

        "TALYS-1.96/2.0. Simulation of nuclear reactions". 
        url:<https://www-nds.iaea.org/talys/tutorials/talys_v1.96.pdf> 

        Notes
        -----        
        The `unique_thread_id` is passed to TALYS as the value of an 
        especially created keyword named ``geffissionfileid`` [1]_ (See 
        the `See Also` section). It is inserted before file suffix `.ff` 
        of perturbed fission fragment files in order to distinguish the 
        file from unperturbed fission fragment files in TALYS library.

        There are in the vicinity of 400 keywords to choose from in the 
        TALYS software. More information about the keywords in TALYS can 
        be found in the TALYS documentation in the "See Also" section.

        References
        ----------
        .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
        GEF and its influence on the nuclear evaporation in TALYS." 
        Technical report UPTEC F 23065, Uppsala University, 2023, accessed
        1 Janauary 2024, url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>.

        Examples
        --------
        Examples generate no output. Creates uniquely named file 
        containing uniquely named TALYS input file at path provided.

        'TMC' simulation (enumeration 5):
        >>> TALYS_Input.write_TALYS_input_file("/local_path/",
                             'TMC_5_input.in', 'TMC_5', 92 , 236, 2.53e-8)
        []

         'Single_Parameters' simulation (enumeration 5):
        >>> TALYS_Input.write_TALYS_input_file("/local_path/",
          '_P_Z_Curve_S2_5_input.in', '_P_Z_Curve_S2_5', 92, 236, 2.53e-8)
        []
        """            
        #--------------------------------------------------- METHOD VARIABLES ------------------------------------------------------#
        elements = {'92':'U','94':'Pu'}                                                 # This list can be extended for the investigation of other elements.
        """Dictionary containing mass number and symbol of chemical elements
        used as targets in simulation of fission.

        This dictionary can be extended for investigations of other elements.
        """
        #------------------------------------------------ TALYS INPUT FILE KEYWORDS --------------------------------------------#
        projectile  = 'n'
        element     = elements[str(Z_target)]
        mass        = str(A_compound-1)                                                 # Minus 1 because GEF uses compound mass and TALYS uses target mass.
        energy      = str(float(E_reaction))
        fission     = 'y'
        ejectiles   = 'g n'
        massdis     = 'y'
        fymodel     = '4'
        ffmodel     = '1' 
        elow        = '0.000001'
        Rfiseps     = '0.000000001'
        outspectra  = 'y'
        bins        = '100'
        channels    = 'n'
        maxchannel  = '8'
        Rspincutff  = '4'
        Rspincut    = '0.4'
        geffissionfileid   = None       # New keyword for passing unique file name to TALYS. See the `See Also` section.
        #-----------------------------------------------------------------------------------------------------------------------#
        with open(os.path.join(folder,file_name),'w') as TALYS_input_file:            
                        TALYS_input_file.write( 
                                        f'projectile {projectile}\n'            \
                                        f'element {element}\n'                  \
                                        f'mass {mass}\n'                        \
                                        f'energy {energy}\n'                    \
                                        f'fission {fission}\n'                  \
                                        f'ejectiles {ejectiles}\n'              \
                                        f'massdis {massdis}\n'                  \
                                        f'fymodel {fymodel}\n'                  \
                                        f'ffmodel {ffmodel}\n'                  \
                                        f'elow {elow}\n'                        \
                                        f'Rfiseps {Rfiseps}\n'                  \
                                        f'outspectra {outspectra}\n'            \
                                        f'bins {bins}\n'                        \
                                        f'channels {channels}\n'                \
                                        f'maxchannel {maxchannel}\n'            \
                                        f'Rspincutff {Rspincutff}\n'            \
                                        f'Rspincut {Rspincut}\n'                \
                                        f'geffissionfileid {unique_thread_id}'  \
                                        )
#------------------------------------------------------- END OF METHOD ---------------------------------------------------------#
                        
#------------------------------------------------------- END OF MODULE ---------------------------------------------------------#