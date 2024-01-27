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
"""Program for simulating fission in GEF using perturbed parameters.

Evaporation of perturbed fission fragment data simulated by both
GEF and TALYS. The script takes user information about: 
- The fission reaction that is to be simulated.
- Which program mode to use.
- If TALYS simulations should be included. 
- GEF and TALYS input data.
- Output path for results. 
The script stores the results in a `pickle` file.
The script deletes all GEF and TALYS simulation files and folders
after the results are stored in the ``Reaction`` object. 

See Also
--------
``package_McPUFF.GEF_input``
``package_McPUFF.TALYS_Input``
``package_McPUFF.Gaussian_GEF_Param``
``package_McPUFF.Parameters_to_vary``
"""
from package_McPUFF import Gaussian_GEF_Param
from package_McPUFF import GEF_input
from package_McPUFF import Parameters_to_vary
from package_McPUFF import TALYS_Input
from subprocess import PIPE, DEVNULL
from threading import Thread
import concurrent.futures
import numpy as np
import subprocess
import shutil
import pickle
import math
import sys
import os
class Reaction:
    """Main object that holds all information about reaction being 
    investigated as well as all simulation results.
    
    This object holds all unperturbed simulation results (for comparison).
    The simulation objects are stored in different lists depending on 
    which program mode that is used (See the 'McPUFF_object_structure' in 
    the `See Also` section).

    Parameters
    ----------
    Z_target : `int`
        Atomic number of target chemical element. 
                
        Needed to create file name in prescribed TALYS format of 
        fission fragment files and for GEF simulations.
    A_compound : `int`
        Mass number of compound target nuclei. I.e A + 1.
                
        Needed to create file name in prescribed TALYS format of  
        fission fragment files and for GEF simulations.
    E_reaction : `float`
        Kinetic energy (MeV) of incident neutron in fission event.
                
        Needed to create file name in prescribed TALYS format of
        fission fragment files and for GEF simulations.
    runs_MC : `int`
        Given in the GEF input file as a multiplication factor
        of the default number of Monte Carlo simulations. 
                
        The GEF default value is 10^5 MC simulations.  
    number_of_randoms : `int`
        The number of simulations with randomly perturbed parameters
        that are to be performed with McPUFF.
    dist_flag : `str`
        Specifies if random numbers are to be drawn from "uniform", 
        "normal" or "max-min" distributions.
    pth_GEF_program : `str`
        Local path to modified GEF software [1]_.
    pth_TALYS_program : `str`
        Local path to modified TALYS software [2]_. See the `See Also` section.
    pth_main : `str`
        Local path to main Python script. Used to create McPUFF outout 
        folder.
    with_TALYS : `boolean`
        Flag to perform evaporation simulations using TALYS or not. 
    prog_flag : `str`
        Flag to set program mode (``TMC`` or ``Single_Parameters``).
      
    See Also
    --------
    ``McPUFF_program.Modified_Parameter``
    ``McPUFF_program.TMC_Mod_Param_object``
    ``package_McPUFF.Parameters_to_vary``
    ``McPUFF_program.read_and_clear_GEF_results()``
    ``McPUFF_program.read_and_clear_TALYS_results()``
    McPUFF_object_structure url: <https://github.com/UPTEC-F-23065/McPUFF/blob/ac8c29ce13729ff85d598f8165555891a9e56ea4/ObjectStructureMcPUFF.png>

    Notes
    -----
    For each simulation loop in McPUFF, following the procedure 
    developed in [3]_, fission events that only occur once are removed. 
    For more information, see [1]_. Hence, the actual simulation results 
    are fewer than specified by user. The actual number of simulations in 
    a loop can be determined using ``unperturbed_ignored_events``. The 
    fission yield is presented as a fraction. Using the 
    ``unperturbed_ignored_events``, it is possible to recreate the number
    of times a certain fission event has occured.

    References
    ----------
    .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
    GEF and its influence on the nuclear evaporation in TALYS," 
    Technical report UPTEC F 23065, Uppsala University, 2023, 
    url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>, 
    accessed 1 Janauary 2024.

    .. [2] P. Karlsson, "Modification-of-TALYS-for-TMC-simulations 
    (Version 1.0.0)" [Computer software], GitHub, 2023, url:
    <https://github.com/UPTEC-F-23065/Modification-of-TALYS-for-TMC-simulations>,
    accessed 4 January 2024.

    .. [3] F. Nordström, "Benchmark of the fission channels in TALYS.", 
    Technical Report UPTEC ES 21016, Uppsala university, 2021.

    .. [4] Python Software Foundation. "Built-in Functions---Python 3.12.1 
    documentation", 2024, 
    url: <https://docs.python.org/3/library/functions.html#vars>, 
    accessed 4 January 2024.

    Examples
    --------    
    Example generates no terminal output. Produces `pickle` file at
    output location provided by user after simulation is complete.

    >>> reac = Reaction(92,236,2.53e-8,1e6,500,"normal",
                                 "/local/path/modified/GEF/",
                                    "/local/path/modified/TALYS/",
                                        "/local/path/main/Python/script",
                                                            True,"TMC")
    ["Z92_A236_n_E2.53e-08MeV.pkl"]
    """
    #--------------------------------------------------- CLASS ATTRIBUTES ------------------------------------------------------#
    distribution_flag = None
    """Determines simulation mode (`str`)"""
    With_TALYS_flag = None 
    """`True` for including simulation of evaporation process with TALYS
    (`boolean`)."""
    program_flag = None 
    """Flag for determining which program mode to use in simulations 
    (`str`)."""
    dictionary_unpert_param_name_val = None 
    """Dictionary of parameter name(s) and default parameter value(s) 
    (`dict` [`str`, `float`]).
    
    Parameters to perturb are chosen by user in 
    ``package_McPUFF.Parameters_to_vary``."""
    list_of_Mod_Param_objects = None
    """Holds objects created in the ``Single_Parameters`` mode 
    (`list` [``McPUFF_program.Modified_Parameter``])."""
    list_of_TMC_Objects = None 
    """Holds objects created in the ``TMC`` mode 
    (`list` [``McPUFF_program.TMC_Object``])."""
    unperturbed_FY = None 
    """Holds unperturbed fission fragment yield GEF
    results, dtype = `float32`. 

    For GEF 'lmd'-option:  (`numpy.ndarray`, (300,11)). 
    Holds: Z, A, Yff, avg(TXE), avg(TKE), E* and std(E*) for light and 
    heavy fission fragment per fission event.
    
    For GEF 'lmd+'-option: (`numpy.ndarray`, (300,19)).
    Holds: Z, A, Yff, avg(TXE), avg(TKE), E* and std(E*), E(n), std(E(n)), 
    E(g), std(E(g)) for light and heavy fission fragment per fission event.
    
    For more information, see [2]_ and [3]_."""
    unperturbed_ignored_events = None 
    """Fission events in GEF 'lmd' file that are removed from results 
    (`int`).

    See `Notes` for more information."""
    unperturbed_TALYS_results = None 
    """The keys of this dictionary are subject to change, use the Python
    command `vars(object)` [4]_ to see the contents (`dict`). 

    See also ``McPUFF_program.read_and_clear_TALYS_results()``."""
    unperturbed_GEF_results = None  
    """The keys of this dictionary are subject to change, use the Python
    command `vars(object)` [4]_ to see the contents (`dict`). 

    See also ``McPUFF_program.read_and_clear_GEF_results()``."""
    reaction_info = None
    """Holds information about:``Mc_runs``,``Z_target``,``A_compound``
    ,``E_reaction``,``num_of_rand`` 
    (`dict` [`int`, `int`, `int`, `float`, `int`])."""
    #------------------------------------------------------ CONSTRUCTOR --------------------------------------------------------#
    def __init__(self,Z_target,A_compound,E_reaction,runs_MC,number_of_randoms,dist_flag,pth_GEF_program,pth_TALYS_program,pth_main,with_TALYS,prog_flag):   
        self.distribution_flag = dist_flag                                              # String: Determines which distribution random numbers are drawn from.
        self.With_TALYS_flag = with_TALYS                                               # Boolean:Determines if run TALYS too
        self.program_flag = prog_flag                                                   # String: Determines if single_parameter or TMC-run
        self.dictionary_unpert_param_name_val = {}                                      # Dictionary: To store parmeter names and unperturbed parameter values.
        self.list_of_Mod_Param_objects = []                                             # List: To store modified parameter objects using append.
        self.list_of_TMC_Objects = []                                                   # List: To store TMC objects using append.
        self.unperturbed_FY = None                                                      # np.array: To store unperturbed fission yields (Declared later)
        self.unperturbed_ignored_events = None                                          # Int: Multichance fission events are removed. Number for each simulations used for processing data.
        self.unperturbed_TALYS_results = None                                           # Dictionary: To store unperturbed TALYS results (Declared later)
        self.unperturbed_GEF_results = None                                             # Dictionary: To store unperturbed GEF results (Declared later)
        self.reaction_info = {"MC_runs":int(runs_MC),                                   # Dictionary: To store reaction information used by functions.
                              "Z_target":int(Z_target),
                              "A_compound":int(A_compound),
                              "E_reaction":float(E_reaction),
                              "num_of_rand":int(number_of_randoms)}
        #------------------------------------------------- CREATE PATH DICTIONARY ----------------------------------------------#
        path_dict = Reaction.create_paths_and_directories(Z_target,A_compound,pth_GEF_program,pth_TALYS_program,pth_main)  
        """Create dictionary of local paths to output folder, modified GEF
        and modified TALYS."""     
        #------------------------------------------- DICTIONARY OF PARAMETER KEY-VALUES ----------------------------------------#
        param_to_vary = Parameters_to_vary.parameters_to_vary()  
        """Reads which of all 94 parameters to use.
        
        GEF parameters in "MyParameters_sample.dat" file: 
        [5,8,10,11,12,14,16,20,24,26,27,28,30,38,39,40,41,48,49,50,51,52,
        53,54,56,60,65,81,82,83,84,85,86,89,96]"""
        dict_unpert_param_name_val = {}                                         
        with open(path_dict['GEF_program_path'] +'Parameters.bas', 'r') as f:           # Open parameter source file in GEF. 
            list_of_parameters = f.readlines()                                          # Create list of lines in file.
            for i in param_to_vary:  
                if i!=55 and i!=98:                                                     # SKip lines without parameters in source file.
                    rad=list_of_parameters[i].split()                                   # Split line into separate strings.
                    dict_unpert_param_name_val[rad[0]]=rad[2]                           # Store name and unperturbed value.
        number_of_parameters = len(param_to_vary)                                       # Number of parameters used. Needed to determine how many CPU's to use (see below).
        self.dictionary_unpert_param_name_val = dict_unpert_param_name_val              # Store in Reaction object.
        del f,i, list_of_parameters, param_to_vary                            
        #------------------------------------ CREATE UNPERTURBED FISSION YIELDS (FY) FOR COMPARISON ----------------------------#
        try:
            unperturbed_thread = Custom_Thread(target=Reaction.create_unperturbed_FY, args=[Z_target,A_compound,E_reaction,runs_MC,with_TALYS,path_dict]) 
            """Create Custom thread for unperturbed FY. 

            Must be custom to be able to retrieve results."""      
            unperturbed_thread.start()                                                                                                                         
        except Exception as e:
            sys.exit(e)                                                                                                                                         
        #------------------------------------ CREATE PERTURBED FY IN 'SINGLE_PARAMETERS' MODE ----------------------------------#
        if self.program_flag == 'Single_Parameters':                                    # Program slow if all processors used. Computer internal processing can use available CPU's for multi-thread processing. 
            number_of_workers = math.floor((os.cpu_count()-2)/3)                        # Determines number of parameters simulated in parallel.
            #----------- PERFORM MULTI-THREAD SIMULATIONS ASYNCHRONOUSLY USING PYTHONS 'CONCURRENT FUTURES' MODULE -------------#
            try:
                with concurrent.futures.ThreadPoolExecutor(max_workers=number_of_workers) as parameter_executor:   
                    future_mod_param_obj = {parameter_executor.submit(Reaction.perturbed_calculations_single_parameter,key,dict_unpert_param_name_val,
                                                                    Z_target,A_compound,E_reaction,runs_MC,number_of_randoms,number_of_workers,with_TALYS,dist_flag,path_dict): key for key in dict_unpert_param_name_val}  
                future_done = concurrent.futures.as_completed(future_mod_param_obj)     # Completion time of simulations vary. Collect results asynchronously.
                for mod_param_obj in future_done:                                       # Simulation results for each parameter collected in its own 'Modified_Parameter' object.
                    self.list_of_Mod_Param_objects.append(mod_param_obj.result())       # 'Modified_Parameter' objects stored in list in main 'Reaction' object.    
            except Exception as e:
                print('Encountered a problem in Threadpool in main.\n')
                sys.exit(e)
            del number_of_parameters,number_of_workers,parameter_executor,future_mod_param_obj,future_done
        #----------------------------------------------- CREATE PERTURBED FY IN 'TMC' MODE -------------------------------------#
        elif self.program_flag == 'TMC':                                                # Program slow if all processors used. Computer internal processing can use available CPU's for multithread processing.
            max_multithreads_TMC = math.floor(os.cpu_count()*(2/3))                     # One GEF-run per random number. All param at once. Use 2/3 of available CPU's.
            #----------- PERFORM MULTI-THREAD SIMULATIONS ASYNCHRONOUSLY USING PYTHONS 'CONCURRENT.FUTURES' MODULE -------------#
            try:                                                
                with concurrent.futures.ThreadPoolExecutor(max_workers=max_multithreads_TMC) as perturbed_TMC_executor:   # Max_workers is really max-2 because main and unpert thread.
                    future_tmc_obj = {perturbed_TMC_executor.submit(Modified_Parameter.create_perturbed_TMC_FY,dict_unpert_param_name_val,n,Z_target,A_compound,
                                        E_reaction,runs_MC,path_dict,with_TALYS,dist_flag):n for n in range(int(number_of_randoms))}
                future_done = concurrent.futures.as_completed(future_tmc_obj)           # Completion time of simulations vary. Collect results asynchronously.
                for tmc_obj in future_done:                                             # Simulation results for all parameters collected in 'TMC_Object' object.
                    self.list_of_TMC_Objects.append(tmc_obj.result())                   # 'TMC_Object' objects stored in list in main 'Reaction' object.
            except Exception as e:
                sys.exit(e) 
            del max_multithreads_TMC,perturbed_TMC_executor,future_tmc_obj,future_done   
        else:
            print('Incorrect program flag- Exiting program')
            sys.exit()
        #------------------------------------------ COLLECT DATA FROM UNPERTURBED THREAD ---------------------------------------#
        unperturbed_thread.join()                                                       # Main thread waits for unperturbed thread to complete.
        self.unperturbed_FY = unperturbed_thread.unperturbed_FY
        self.unperturbed_GEF_results = unperturbed_thread.unperturbed_GEF_results
        self.unperturbed_TALYS_results = unperturbed_thread.unperturbed_TALYS_results
        self.unperturbed_ignored_events = unperturbed_thread.unperturbed_ignored_events
        del unperturbed_thread
        #--------------------------------------------------- PICKLE RESULTS ----------------------------------------------------#
        with open(os.path.join(path_dict['PKL_path'],f'Z{Z_target}_A{A_compound}_n_E{E_reaction}MeV.pkl'), 'wb') as Reac_object:
            pickle.dump(self, Reac_object)   
        del Reac_object
    #------------------------------------------------------ CLASS METHODS ------------------------------------------------------#
    def clear_MyParameters_dat(GEF_cwd_path):
        """Writes empty string into file, replacing whatever was there to
        make room for new perturbed parameter values.
        
        Parameters
        ----------
        GEF_cwd_path : `str`
            Path to modified GEF current working directory.
        
        Returns
        -------
        Function has no return value.

        Examples
        --------
        Function replaces content of 'MyParameters.dat' file with 
        empty string.
        >>> Reaction.clear_MyParameters_dat("/local/path/modified/GEF")
        []
        """
        #-----------------------------------------------------------------------------------------------------------------------#
        path_MyParameters_file = os.path.join(GEF_cwd_path,'MyParameters.dat')
        if os.path.isfile(path_MyParameters_file):
            with open(path_MyParameters_file, 'w') as MyParam_dat_file:
                MyParam_dat_file.write('')
            del path_MyParameters_file,MyParam_dat_file
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#
            
    def create_paths_and_directories(Z_target,A_compound,pth_GEF_program,pth_TALYS_program,pth_main):
        """Function creates directories and path strings needed for 
        simulations.

        Parameters
        ----------
        Z_target : `int`
                Atomic number of target chemical element. 
                
                Needed to create file name in prescribed TALYS format of 
                fission fragment files in GEF ''.ff'-file library.
        A_compound : `int`
                Mass number of compound target nuclei. I.e A + 1. 
                
                Needed to create file name in prescribed TALYS format of 
                fission fragment files in GEF ''.ff'-file library.
        pth_GEF_program : `str`
            Local path to modified GEF software.
        pth_TALYS_program : `str`
            Local path to modified TALYS software.
        pth_main : `str`
            Local path to McPUFF main script.

        Returns
        -------
        path_dict : `dict` 
            path_dict dictionary with keys:

            ``"GEF_program_path"``
                Path to modified GEF executable (`str`).
            ``"TALYS_program_path"``
                Path to modified TALYS executable (`str`).
            ``"Output_McPUFF_path"``
                Path to McPUFF output folder (`str`).
            ``"GEF_working_dir_path"``
                Path to modified GEF output folder in McPUFF output folder
                (`str`).
            ``"TALYS_working_dir_path"``
                Path to modified TALYS output folder in McPUFF output 
                folder (`str`).
            ``"TALYS_input_path"``
                Path to where to place modified TALYS input file (`str`).
            ``"PKL_path"``
                Path to where to place `pickle` file with McPUFF results
                (`str`).
            ``"TALYS_ff_file_path"``
                Path to modified GEF library of ".ff" files (`str`).

            Notes
            -----
            The ``TALYS_ff_file_path`` below must be in the format of the 
            default ".ff" files in the GEF library. The name includes
            the symbol of the chemical element. This `if` statement
            can be extended if simulations are to be performed for
            other elements.
            
            Examples
            --------
            >>> path_dict = Reaction.create_paths_and_directories(92,236,
                               "/path/mod/GEF","/path/mod/TALYS","/path/McPUFF")
            [path_dict = {"GEF_program_path":"/local_path/","TALYS_program_path":
                        "/local_path/","Output_McPUFF_path":"/local_path/",
                        "GEF_working_dir_path":"/local_path/","TALYS_working_dir_path":
                        "/local_path/","TALYS_input_path":"/local_path/",
                        "PKL_path":"/local_path/","TALYS_ff_file_path":"/local_path/"}]
        """
        #-----------------------------------------------------------------------------------------------------------------------#
        GEF_program_path = os.path.join(pth_GEF_program,"")                                     # Path GEF executable. The empty string at the end adds a '/'
        TALYS_program_path = os.path.join(pth_TALYS_program,"")		                            # Path TALYS executable
        if not os.path.isdir(os.path.join(os.path.dirname(pth_main),'Output_McPUFF')):          # If not exist, create folder for storing all program output.     
            prog_out_dir_name = os.path.join(os.path.dirname(pth_main),'Output_McPUFF')         # Location of output folder can be chosen by user.    
            os.mkdir(prog_out_dir_name)
            Output_McPUFF_path = os.path.join(prog_out_dir_name,"")
            del prog_out_dir_name
        else:
            Output_McPUFF_path = os.path.join(os.path.dirname(pth_main),'Output_McPUFF',"")    
        if not os.path.isdir(os.path.join(Output_McPUFF_path,'GEF_working_directory')):			# If not exist, create GEF-working_directory and create path
            dir_path = os.path.join(Output_McPUFF_path,'GEF_working_directory')			        # Creates GEF output folder within McPUFF output folder.
            os.mkdir(dir_path)
            GEF_working_dir_path = os.path.join(dir_path,"")
            del dir_path
        else:
           GEF_working_dir_path = os.path.join(Output_McPUFF_path,'GEF_working_directory',"")   
        if not os.path.isdir(os.path.join(Output_McPUFF_path,'TALYS_working_directory')):       # If not exist, create TALYS-working_directory and create path
            dir_path = os.path.join(Output_McPUFF_path,'TALYS_working_directory')               # Creates TALYS output folder within McPUFF output folder.
            os.mkdir(dir_path)
            TALYS_working_dir_path = os.path.join(dir_path,"")
            del dir_path
        else:
            TALYS_working_dir_path = os.path.join(Output_McPUFF_path,'TALYS_working_directory',"") 		
        TALYS_input_path = os.path.join(TALYS_working_dir_path,'TALYS_folder_')                 # Create path to where to place TALYS input file.          
        if not os.path.isdir(os.path.join(Output_McPUFF_path,'pickle_results')):                # If not exist, create pickle file folder and create path.
            dir_path = os.path.join(Output_McPUFF_path,'pickle_results')                        # Creates pickle file folder within McPUFF output folder.
            os.mkdir(dir_path)
            PKL_path = os.path.join(dir_path,"")
            del dir_path
        else:
            PKL_path = os.path.join(Output_McPUFF_path,'pickle_results')
        PKL_path = os.path.join(Output_McPUFF_path,'pickle_results',"")
        # This if-statement can be extended to include other elements. See 'Notes'.
        if Z_target == 92:                                                                       
            TALYS_ff_file_path = os.path.join(os.path.dirname(pth_TALYS_program),'structure','fission','ff','gef',f'U{str(A_compound)}',"")
        if Z_target == 94:
            TALYS_ff_file_path = os.path.join(os.path.dirname(pth_TALYS_program),'structure','fission','ff','gef',f'Pu{str(A_compound)}',"")
        #----------------------------------------- CREATE DICTIONARY TO RETURN ----------------------------------------------------#
        path_dict = {'GEF_program_path':GEF_program_path,'TALYS_program_path':TALYS_program_path,'Output_McPUFF_path':Output_McPUFF_path,
                     'GEF_working_dir_path':GEF_working_dir_path,'TALYS_working_dir_path':TALYS_working_dir_path,'TALYS_input_path':TALYS_input_path,
                     'PKL_path':PKL_path,'TALYS_ff_file_path':TALYS_ff_file_path}
        return path_dict
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#

    def create_unperturbed_FY(Z_target,A_compound,E_reaction,MC_runs,With_TALYS_flag,path_dict):
        """Perform simulation with unperturbed GEF parameters.
        
        Parameters
        ----------
        Z_target : `int`
            Atomic number of target chemical element. 
                
            Needed to create file name in prescribed TALYS format of 
            fission fragment files and for GEF simulations.
        A_compound : `int`
            Mass number of compound target nuclei. I.e A + 1.
                    
            Needed to create file name in prescribed TALYS format of  
            fission fragment files and for GEF simulations.
        E_reaction : `float`
            Kinetic energy (MeV) of incident neutron in fission event.
                    
            Needed to create file name in prescribed TALYS format of
            fission fragment files and for GEF simulations.
        MC_runs : `int`
            Given in the GEF input file as a multiplication factor
            of the default number of Monte Carlo simulations. 
                    
            The GEF default value is 10^5 MC simulations. 
        With_TALYS_flag : `boolean`
            If `True`, includes TALYS evaporation simulation. 
        path_dict : `dict` [`str`,`str`,`str`,`str`,`str`,`str`,`str`,`str`]
            Dictionary with local paths needed for simulations.

        Returns
        -------
        unperturbed_result_dict : `dict`
            Dictionary of unperturbed GEF FY results with keys:
                ``'FY_TALYS_format'``
                    FY from GEF in ".ff" file format [1]_
                    For GEF "lmd" option: (`numpy.ndarray`, (300,11)).
                    For GEF "lmd+" option: (`numpy.ndarray`, (300,19)).
                    See the `See Also` section.
                ``'GEF_data_dict'``
                    Dictionary of unperturbed GEF results (`dict`).
                ``'TALYS_data_dict'``
                    Dictionary of unperturbed TALYS results (`dict`).
                ``'unpert_ignored_events'`` 
                    Number of removed GEF fission events because they only
                    occur once (for statistical reasons) [1]_ (`int`).

        Notes
        -----
        The ``unique_unpert_thread_ID`` is important. It is used to keep
        simulation results separate when multi-threading and to 
        distinguish perturbed ".ff" files in the GEF library [1]_.

        In case ``With_TALYS_flag`` = "False", the function returns 
        ``TALYS_data_dict`` as an empty dictionary.

        See Also
        --------
        ``Reaction.create_paths_and_directories()``
        ``Reaction.read_and_clear_GEF_results()``
        ``Reaction.read_and_clear_TALYS_results()``
        ``Reaction.FY_results()``
        ``Reaction.GEF_FY_for_TALYS()``

        References
        ----------
        .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
        GEF and its influence on the nuclear evaporation in TALYS." 
        Technical report UPTEC F 23065, Uppsala University, 2023, 
        url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>, 
        accessed 1 Janauary 2024.

        Examples
        --------
        >>> Unpert_FY = Reaction.create_unperturbed_FY(92,236,
                                             2.53e-8,1e6,'True',path_dict)
        [{'FY_TALYS_format':(numpy.ndarray,(300,11)),
                   'GEF_data_dict':dict,'TALYS_data_dict':dict,
                                           'unpert_ignored_events':int}]
        """
        #-------------------------------------------------- TALYS OUTPUT DICT --------------------------------------------------#
        TALYS_data_dict = {}
        #--------------------------------------------- CREATE GEF AND TALYS PATHS ----------------------------------------------#
        unique_unpert_thread_ID = 'Unperturbed_Param'                                                       # Unique thread identifier.
        GEF_working_dir_path = str(path_dict['GEF_working_dir_path'])
        TALYS_ff_file_path = str(path_dict['TALYS_ff_file_path'])
        TALYS_input_path = str(path_dict['TALYS_input_path']) 
        #------------------------------------------------ RUN GEF SIMULATION ---------------------------------------------------# 
        GEF_data_paths = Reaction.create_GEF_workingdir_and_inputfile(unique_unpert_thread_ID,Z_target,A_compound,E_reaction,MC_runs,GEF_working_dir_path)              
        GEF_cwd_path    = GEF_data_paths['GEF_cwd_path']                                                    # Local path to individual unperturbed GEF simulation output folder.
        GEF_LMD_path    = GEF_data_paths['GEF_LMD_path']                                                    # Local path to individual unperturbed GEF simulation ".lmd" file.
        GEF_DAT_path    = GEF_data_paths['GEF_DAT_path']                                                    # Local path to individual unperturbed GEF simulation ".dat" file.
        GEF_DMP_EN_path = GEF_data_paths['GEF_DMP_EN_path']                                                 # Local path to individual unperturbed GEF simulation "EN.dat" file.
        Reaction.clear_MyParameters_dat(GEF_cwd_path)                                                       # Clears MyParameters.dat of previous values.
        Reaction.eraseFolders(GEF_cwd_path)                                                                 # Make sure folders with old results are erased before start.
        with subprocess.Popen(["GEF"],stdin=PIPE,stdout=DEVNULL,cwd=GEF_cwd_path) as GEF_process:   
            GEF_process.communicate(input=b'\n')
        FY_TALYS_format,unpert_ignored_events = Reaction.FY_results(GEF_LMD_path,A_compound,Z_target)
        GEF_data_dict = Reaction.read_and_clear_GEF_results(GEF_cwd_path,GEF_DAT_path,GEF_DMP_EN_path)      # Reads and stores data from GEF .dat file
        del GEF_cwd_path,GEF_LMD_path,GEF_DAT_path,GEF_process
        #------------------------------------------------ RUN TALYS SIMULATION ---------------------------------------------------#
        if With_TALYS_flag == True:            
            Reaction.print_TALYS_ff_files(TALYS_ff_file_path,FY_TALYS_format,Z_target,A_compound,E_reaction,unique_unpert_thread_ID)
            TALYS_unpert_input_command = 'talys' 
            TALYS_unpert_input_file    = f'{unique_unpert_thread_ID}_input.in' 
            TALYS_unpert_output_file   = f'{unique_unpert_thread_ID}_output.out' 
            TALYS_unpert_cwd           = os.path.join(TALYS_input_path+f'{unique_unpert_thread_ID}',"")     # Empty string at end add a '/' to file name.
            TALYS_Input.create_TALYS_input_file(TALYS_input_path,unique_unpert_thread_ID,TALYS_unpert_input_file,Z_target,A_compound,E_reaction)
            try: 
                subprocess.run(["bash","-c",str(f"{TALYS_unpert_input_command} < {TALYS_unpert_input_file} > {TALYS_unpert_output_file}")],stdout=DEVNULL,cwd=TALYS_unpert_cwd)
            except FileNotFoundError as e:
                print(f'TALYS could not run because it cannot find the input file.\n{e}\n')
                sys.exit(e)
            except subprocess.CalledProcessError as e:
                print(f'An error occured while running TALYS\n{e}\n')
                sys.exit()
            TALYS_data_dict = Reaction.read_and_clear_TALYS_results(TALYS_unpert_cwd,TALYS_ff_file_path,unique_unpert_thread_ID) 
            del TALYS_unpert_input_command,TALYS_unpert_input_file,TALYS_unpert_output_file,TALYS_unpert_cwd
        #----------------------------------------------- IMPORTANT INDENTATION -------------------------------------------------#
        unperturbed_result_dict = {'FY_TALYS_format':FY_TALYS_format,'GEF_data_dict':GEF_data_dict,'TALYS_data_dict':TALYS_data_dict,'unpert_ignored_events':unpert_ignored_events}
        return unperturbed_result_dict 
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#

    def delete_GEF_result_folder(GEF_results_folder_path):     
        """Delete GEF output data from simulation.
        
        After GEF output is stored in object it is deleted to save space.
        
        Parameters
        ----------
        GEF_results_folder_path : `str`
            Local path to GEF output folder.

        Returns
        -------
        Function has no return value.

        Examples
        --------
        Function deletes files and folders at specified location.

        >>> Reaction.delete_GEF_result_folder("local/path/GEF/output/")
        []
        """      
        #-----------------------------------------------------------------------------------------------------------------------#
        if os.path.isdir(GEF_results_folder_path):
            shutil.rmtree(GEF_results_folder_path)    
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#
            
    def delete_TALYS_ff_files(pth_TALYS_FY_folder,unique_thread_ID):
        """Delete perturbed ff-file from GEF library after completed TALYS
        simulation.
        
        Parameters
        ----------
        pth_TALYS_FY_folder : `str`
            Local path to GEF ff-file folder in GEF library.
        unique_thread_ID : `str`
            Used to identify individual perturbed GEF ff-file when 
            multi-threading.
        
        Returns
        -------
        Function has no return value.

        Notes
        -----
        Function singles out files ending with the unique ID + ".ff" and 
        deletes them, leaving all original files intact.
        TALYS has an internal function that converts all user input in 
        input file to lower case. Hence, when constructing the file name
        used to identify the individual perturbed GEF ff-file, the 
        "unique_thread_ID" has to be converted to lower case. 

        Examples
        --------
        Function deletes perturbed GEF ff-file in GEF library.

        >>> Reaction.delete_TALYS_ff_files("/local/path/TALYS/output/","TMC_5")
        []
        """
        #-----------------------------------------------------------------------------------------------------------------------#
        for file_name in os.scandir(pth_TALYS_FY_folder):
            if file_name.is_file() and file_name.name.endswith(str(unique_thread_ID).lower()+'.ff'):
                os.remove(os.path.join(pth_TALYS_FY_folder,file_name.name))
        del file_name
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#
        
    def delete_TALYS_result_files(path_TALYS_cwd):
        """Delete TALYS output data from simulation.
        
        After TALYS output is stored in object it is deleted to save space.
        
        Parameters
        ----------
        path_TALYS_cwd : `str`
            Local path to TALYS output folder.

        Returns
        -------
        Function has no return value.

        Examples
        --------
        Function deletes files and folders at specified location.

        >>> Reaction.delete_TALYS_result_files("/local/path/TALYS/output/")
        []
        """      
        #-----------------------------------------------------------------------------------------------------------------------#
        if os.path.isdir(path_TALYS_cwd):
            shutil.rmtree(path_TALYS_cwd)
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#
            
    def eraseFolders(GEF_cwd_path):
        """Delete GEF default output folders after completed GEF simulation.
        
        Parameters
        ----------
        GEF_cwd_path : `str`
            Local path to modified GEF software.

        Returns
        -------
        Function has no return value.

        See Also
        --------
        GEF 2023-V1.1 "ReadMe" file:<https://www.khschmidts-nuclear-web.eu/GEF_code/GEF-2023-1-1/Standalone/Readme.txt>

        Notes
        -----
        GEF produces several default output folders containg various
        simulation data. After completed GEF simulation these folders
        must be erased, otherwise GEF exits during the next simulation
        loop because the folders indicate that the simulation has already
        been performed. This behaviour is a part of the GEF softwares way
        of running multiple simulations specified in a list. (See the GEF
        documentation for more information.)
        
        Examples
        --------
        Function deletes files and folders at specified location.

        >>> Reaction.eraseFolders("/local/path/modified/GEF/")
        []
        """
        #-----------------------------------------------------------------------------------------------------------------------#
        if os.path.isdir(os.path.join(GEF_cwd_path,'ctl')):
            shutil.rmtree(os.path.join(GEF_cwd_path,'ctl',""))
        if os.path.isdir(os.path.join(GEF_cwd_path,'dmp')):
            shutil.rmtree(os.path.join(GEF_cwd_path,'dmp',""))
        if os.path.isdir(os.path.join(GEF_cwd_path,'out')):
            shutil.rmtree(os.path.join(GEF_cwd_path,'out',""))
        if os.path.isdir(os.path.join(GEF_cwd_path,'tmp')):
            shutil.rmtree(os.path.join(GEF_cwd_path,'tmp',""))
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#
            
    def FY_results(LMD_path,A_compound,Z_target):
        """Convert "raw" GEF fission fragment yield data from the GEF 
        ".lmd" or ".lmd+" file into GEF ".ff" file format used in GEF 
        library.
        
        Parameters
        ----------
        LMD_path : `str`
            Local path to GEF ".lmd" file from perturbed GEF simulation.
        A_compound : `int`
            Mass number of compound target nuclei. I.e A + 1.
                    
            Needed to be able to delete multi-chance fission events.
            (See the "See Also" section).
        Z_target : `int`
        Atomic number of target chemical element. 
                
        Needed to identify lines of fission fragment yield data in
        ".lmd+" file.

        Returns
        -------
        With GEF "lmd" option: 
        FY_TALYS_format : `numpy.ndarray` (number of GEF MC simulations),11)
            dtype = (`numpy.float32`). GEF FY in GEF library format [1]_.
        With GEF "lmd+" option: 
        FY_TALYS_format : `numpy.ndarray` (number of GEF MC simulations),19)
            dtype = (`numpy.float32`). GEF FY in GEF library format [1]_. 
                          
            The GEF "lmd+" option includes neutron and gamma 
            energies per fragment. See the `Notes` section.
        ignored_events : `int`
            Number of multi-chance fission events removed from an 
            individual GEF simulation.

        Raises
        ------
        UserWarning
            Raised if there is no ".lmd" file at specified location, or
            if there is and the file is empty.

        See Also
        --------
        ``Reaction.GEF_FY_for_TALYS()``
        ``Reaction. runs_MC``
        GEF 2023-V1.1 "ReadMe" file:<https://www.khschmidts-nuclear-web.eu/GEF_code/GEF-2023-1-1/Standalone/Readme.txt>

        Notes
        -----
        The path to the ".lmd" and the ".lmd+" files are the same. The 
        ".lmd+" version holds information about the energies of the
        prompt neutrons and prompt gammas as a function of A.

        When the user chooses the GEF option "lmd+", the ".lmd" file 
        changes format and has to be handled differently. This is checked
        for by a "try-catch" clause. If the `numpy.loadtxt()` encounters
        an uneven number of columns it is an ".lmd+" file and is 
        pre-processed by this function.

        The return value "FY_TALYS_format" is an numpy array with a number 
        of lines equal to the number of GEF Monte Carlo simulations, e.g 
        1e6 (See the "See Also" section) specified by the user. 
        For the GEF "lmd" option it has 11 columns ([0-10]). They contain: 
        [0] = Z1, [1] = A1, [2] = Z2, [3] = A2, [4] = Yield, [5] = TKE(pre),
        [6] = TXE, [7] = Eexc1, [8] = W1, [9] = Eexc2, [10] = W2, where
        Eexc = E*, W = std(E*) and 1,2 are fragments 1 and 2. 
        For the GEF "lmd+" option it has 19 columns ([0-10]). They contain:
        [0] = Z1, [1] = A1, [2] = Z2, [3] = A2, [4] = Yield, [5] = TKE(pre),
        [6] = TXE, [7] = Eexc1, [8] = W1, [9] = Eexc2, [10] = W2, 
        [11] = E(n)_light_frag, [12] = std(E(n)_light_frag), 
        [13] = E(n)_heavy_frag, [14] = std(E(n)_heavy_frag), 
        [15] = E(g)_light_frag, [16] = std(E(g)_light_frag), 
        [17] = E(g)_heavy_frag, [17] = std(E(g)_heavy_frag),  where
        Eexc = E*, W = std(E*) and 1,2 are fragments 1 and 2. 

        The ".ff" file format for the TALYS fission fragment yield
        library requires that multi-chance fission events and fission
        events that only occur once, are removed [1]_. The variable
        ``ignored_events`` represents the number of fission events
        removed because they only occured once. It is saved because with 
        it, the number of multi-chance fission events that were removed
        can be calculated together with the total number of GEF fission
        events (which is given by the variable ``runs_MC``).

        References
        ----------
        .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
        GEF and its influence on the nuclear evaporation in TALYS." 
        Technical report UPTEC F 23065, Uppsala University, 2023, accessed
        1 Janauary 2024, 
        url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>.

        Examples
        --------
        For this example, the "lmd" option is chosen and 35 fission events 
        in GEF occured only once.

        >>> FY, ignoredEvents = Reaction.FY_results("/local/LMD/path",236)
        [numpy.ndarray (number of GEF MC simulations, 11), 35]

        For this example, the "lmd+" option is chosen and 42 fission events 
        in GEF occured only once.

        >>> FY, ignoredEvents = Reaction.FY_results("/local/LMD/path",236)
        [numpy.ndarray (number of GEF MC simulations, 19), 42]
        """
        #--------------------------------------------- READ GEF OUTPUT DATA-----------------------------------------------------#
        try: 
            if os.path.isfile(LMD_path) and os.path.getsize(LMD_path)!=0:                                               # Check if ".lmd" file exists and is not empty.
                try:                                                                                                    # Check if the user chose GEF option "lmd" or "lmd+".
                    lmdData = np.loadtxt(LMD_path,dtype=np.float32,comments='*',usecols=(2,3,4,5,18,19,22))
                    FY = np.zeros((len(lmdData[:,0]),7),dtype=np.float32)                                               # Allocate array for simulation with "lmd" option (without neutron and gamma energies per fragment).
                except Exception as e:                                                                                  # Pre-processes ".lmd+" file.
                    with open(LMD_path,'r') as LMD_plus_file:
                        All_lmd_plus_data = LMD_plus_file.readlines()     
                    #------------------------------- FY DATA WITHOUT GAMMA RAY ENERGES -----------------------------------------#                
                    lmd_plus_data = [line.strip() for line in All_lmd_plus_data if line.strip().startswith(str(Z_target))]              # Singles out lines with FY data (removes energy data).
                    lmdData = np.loadtxt(lmd_plus_data ,dtype=np.float32,comments='*',usecols=(2,3,4,5,18,19,22))                       # FY data identical to when using GEF "lmd" option.
                    FY = np.zeros((len(lmdData[:,0]),11),dtype=np.float32)                                                              # Allocate array for simulation with "lmd+" option (with neutron and gamma energies per fragment).
                    #---------------------------------- NEUTRON ENERGY LIGHT FRAGMENTS -----------------------------------------#
                    E_n_light_frag       = [line.strip()[1:] for line in All_lmd_plus_data if line.strip().startswith(str(1))]          # Lines with energy of emitted n in light fragment (and remove '1').
                    #---------------------------------- NEUTRON ENERGY LIGHT FRAGMENTS -----------------------------------------#
                    E_n_heavy_frag       = [line.strip()[1:] for line in All_lmd_plus_data if line.strip().startswith(str(2))]          # Lines with energy of emitted n in heavy fragment (and remove '2').
                    #------------------------------- GAMMMA RAY ENERGY LIGHT FRAGMENTS -----------------------------------------#
                    E_competition_g_light_frag       = [line.strip()[1:] for line in All_lmd_plus_data if line.strip().startswith(str(3))]    # Lines with energy of emitted g in competition with n light fragment (and remove '3').
                    E_statistical_g_light_frag       = [line.strip()[1:] for line in All_lmd_plus_data if line.strip().startswith(str(4))]    # Lines with energy of statistical g in light fragment after neutrons (and remove '4').
                    E_prompt_collective_g_light_frag = [line.strip()[1:] for line in All_lmd_plus_data if line.strip().startswith(str(5))]    # Lines with energy of prompt collective g light fragment, final state GS (and remove '5').
                    #------------------------------- GAMMMA RAY ENERGY HEAVY FRAGMENTS -----------------------------------------#
                    E_competition_g_heavy_frag       = [line.strip()[1:] for line in All_lmd_plus_data if line.strip().startswith(str(6))]    # Lines with energy of emitted g in competition with n heavy fragment (and remove '6').
                    E_statistical_g_heavy_frag       = [line.strip()[1:] for line in All_lmd_plus_data if line.strip().startswith(str(7))]    # Lines with energy of statistical g in heavy fragment after neutrons (and remove '7').
                    E_prompt_collective_g_heavy_frag = [line.strip()[1:] for line in All_lmd_plus_data if line.strip().startswith(str(8))]    # Lines with energy of prompt collective g heavy fragment, final state GS (and remove '8').
                    #--------- SUM NEUTRON AND GAMMA RAY ENERGIES FOR DIFFERENT CONTRIBUTIONS AND COLLECT IN ARRAY -------------#
                    # Columns: [0]=E(n)_light,[1]=E(n)_heavy,[2]=E(g)_comp_light,[3]=E(g)_stat_light,[4]=E(g)_coll_light,[5]=E(g)_comp_heavy,[6]=E(g)_stat_heavy,[7]=E(g)_coll_heavy
                    E_values_n_g = np.zeros((len(lmdData[:,0]),8))                                                                                                             
                    for i,line in enumerate(E_n_light_frag):                                                                            # ================== Store data for ========================#    
                        location_n_in_string = [loc for loc in range(0,len(line.split()),4)]                                            # Picks out every fourth element in list (neutron energies)
                        E_values_n_g[i,0] = np.sum([line.split()[n] for n in location_n_in_string],dtype=np.float32)                    # [0] = E(n) light fragment.
                        del i,line,location_n_in_string
                    for i,line in enumerate(E_n_heavy_frag):
                        location_n_in_string = [loc for loc in range(0,len(line.split()),4)]                                            # Picks out every fourth element in list (neutron energies)
                        E_values_n_g[i,1] = np.sum([line.split()[n] for n in location_n_in_string],dtype=np.float32)                    # [1] = E(n) heavy fragment.
                        del i,line,location_n_in_string
                    for i,line in enumerate(E_competition_g_light_frag):
                        E_values_n_g[i,2] = np.sum([float(energy) for energy in line.split() if energy[0].isdigit()],dtype=np.float32)  # [2] = E(g) competition light fragment.
                        del i,line
                    for i,line in enumerate(E_statistical_g_light_frag):
                        E_values_n_g[i,3] = np.sum([float(energy) for energy in line.split() if energy[0].isdigit()],dtype=np.float32)  # [3] = E(g) statistical light fragment.
                        del i,line
                    for i,line in enumerate(E_prompt_collective_g_light_frag):
                        E_values_n_g[i,4] = np.sum([float(energy) for energy in line.split() if energy[0].isdigit()],dtype=np.float32)  # [4] = E(g) collective light fragment.
                        del i,line
                    for i,line in enumerate(E_competition_g_heavy_frag):
                        E_values_n_g[i,5] = np.sum([float(energy) for energy in line.split() if energy[0].isdigit()],dtype=np.float32)  # [5] = E(g) competition heavy fragment.
                        del i,line
                    for i,line in enumerate(E_statistical_g_heavy_frag):
                        E_values_n_g[i,6] = np.sum([float(energy) for energy in line.split() if energy[0].isdigit()],dtype=np.float32)  # [6] = E(g) statistical heavy fragment.
                        del i,line
                    for i,line in enumerate(E_prompt_collective_g_heavy_frag):
                        E_values_n_g[i,7] = np.sum([float(energy) for energy in line.split() if energy[0].isdigit()],dtype=np.float32)  # [7] = E(g) collective heavy fragment.
                        del i,line
                   #----------------------------------- ADD NEUTRON AND GAMMA ENERGIES TO LMD DATA ARRAY -----------------------#          
                    FY[:,7]  += E_values_n_g[:,0]                                                                                       # Add energy of neutrons from light fragment to FY array col 7.
                    FY[:,8]  += E_values_n_g[:,1]                                                                                       # Add energy of neutrons from heavy fragment to FY array col 8.
                    FY[:,9]  += E_values_n_g[:,2] + E_values_n_g[:,3] + E_values_n_g[:,4]                                               # Add energy of all gamma emissions from light fragments to FY array col 9.
                    FY[:,10] += E_values_n_g[:,5] + E_values_n_g[:,6] + E_values_n_g[:,7]                                               # Add energy of all gamma emissions from light fragments to FY array col 10.
                    del LMD_plus_file,All_lmd_plus_data,lmd_plus_data,E_n_light_frag,E_n_heavy_frag,E_competition_g_light_frag,E_statistical_g_light_frag,E_prompt_collective_g_light_frag,
                    E_competition_g_heavy_frag,E_statistical_g_heavy_frag,E_prompt_collective_g_heavy_frag,E_values_n_g 
                #-------------------------------- PICK OUT DATA FROM LMD DATA ARRAY --------------------------------------------#    
                Z1     = lmdData[:,0]       # Col:  2  in LMD/LMD+ file.     
                Z2     = lmdData[:,1]       # Col:  3  in LMD/LMD+ file.        
                A1sci  = lmdData[:,2]       # Col:  4  in LMD/LMD+ file.        
                A2sci  = lmdData[:,3]       # Col:  5  in LMD/LMD+ file. 
                Eexc1  = lmdData[:,4]       # Col: 18  in LMD/LMD+ file. 
                Eexc2  = lmdData[:,5]       # Col: 19  in LMD/LMD+ file. 
                TKEpre = lmdData[:,6]       # Col: 22  in LMD/LMD+ file.  
                #----------------------------------- POPULATE FY VECTOR WITH LMD DATA ------------------------------------------#
                FY[slice(0,len(Z1)),0]     = Z1
                FY[slice(0,len(A1sci)),1]  = A1sci 
                FY[slice(0,len(Z2)),2]     = Z2
                FY[slice(0,len(A2sci)),3]  = A2sci
                FY[slice(0,len(TKEpre)),4] = TKEpre 
                FY[slice(0,len(Eexc1)),5]  = Eexc1
                FY[slice(0,len(Eexc2)),6]  = Eexc2
                #------------------------------------- CREATE FY IN ".FF" FILE FORMAT ------------------------------------------#
                num_of_events = len(lmdData[:,0])                                                                                       # Number of GEF simulations. Needed for calculating yields. 
                FY_TALYS_format, ignored_events = Reaction.GEF_FY_for_TALYS(FY,A_compound,num_of_events)
                del lmdData,Z1,A1sci,Z2,A2sci,TKEpre,Eexc1,Eexc2
                return FY_TALYS_format, ignored_events
            else:
                raise UserWarning("UserWarning: loadtxt: Empty input file: "+LMD_path)
        except UserWarning as e:
            sys.exit(e)
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#
            
    def GEF_FY_for_TALYS(FY,A_compound,num_of_events):
        """Assemble perturbed fission fragment yield data from GEF in 
        same ".ff" file format as the default files in the GEF library.
        
        Parameters
        ----------
        FY : `numpy.ndarray` (number of GEF MC simulations,7 or 11)
            "Raw" GEF ".lmd" (7 columns) or ".lmd+" (11 columns) fission 
            fragment yield data with average values for all unique fission 
            events.

            See the "See Also" and "Notes" sections.
        A_compound : `int`
            Mass number of compound target nuclei. I.e A + 1.
                    
            Needed to be able to delete multi-chance fission events.
            See "Notes" section.
        num_of_events : `int`
            Number of GEF Monte Carlo simulations.

        Returns
        -------
        With GEF "lmd" option: 
        FY_TALYS_format : `numpy.ndarray` (300,11)
            dtype = (`numpy.float32`). GEF fission fragment yields in ".ff" 
            file library format [1]_. 

        With GEF "lmd+" option: 
        FY_TALYS_format : `numpy.ndarray` (300,19)
            dtype = (`numpy.float32`). GEF fission fragment yields in ".ff" 
            file library format [1]_. 
                          
            The GEF "lmd+" option includes neutron and gamma 
            energies per fragment. See ``Reaction.FY_results()``.
        ignored_events : `int`
            GEF  fission events that only occur once are removed due to
            statistical reasons [1]_. 

        See Also
        --------
        ``Reaction.FY_results()``
        `numpy.mean()`(url:
        <https://numpy.org/doc/stable/reference/generated/numpy.mean.html>)

        The N-dimensional array (ndarray), 
        (url: <https://numpy.org/doc/stable/reference/arrays.ndarray.html#index-2>)
        
        Notes
        -----
        Multi-chance fission events in GEF are removed. These are events
        where the nucleus has enough excitation energy to emit a particle
        prior to fission. Multi-chance fission events are identified by
        adding the mass number of the two fission fragments and checking
        that the sum matches the compound nucleus mass number. 

        Sometimes GEF swaps places for the light and heavy fission
        fragment in the ".lmd" file. This function checks for that and
        remedies it by swapping back.

        When calculating the mean value and standard deviations for the 
        FY_results, dtype = `numpy.float64`is used because `numpy.float32`
        and lower can lead to erronous results. See the Python 
        documentaion.

        The ``FY_TALYS`` array holding the ".lmd" data from GEF is 
        initialized to 300 rows. Since a numpy array is contiguous, it 
        cannot be dynamically extended during runtime. A typical GEF 
        simulation will produce, in the vicinity of, 200 unique fission 
        events. The array can be initialized to any number that is longer
        that the number of unique fission events.

        References
        ----------
        .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
        GEF and its influence on the nuclear evaporation in TALYS." 
        Technical report UPTEC F 23065, Uppsala University, 2023, accessed
        1 Janauary 2024, 
        url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>.

        Examples
        --------
        For this example, the "lmd" option was chosen and 35 fission events 
        in GEF occured only once.

        >>> Reaction.GEF_FY_for_TALYS(numpy.ndarray (1e6,7), 236, 1e6)
        [numpy.ndarray (300,11), 35]

        For this example, the "lmd+" option was chosen and 42 fission events 
        in GEF occured only once.

        >>> Reaction.GEF_FY_for_TALYS(numpy.ndarray (1e6,7), 236, 1e6)
        [numpy.ndarray (300,19), 42]
        """
        #--------------- CREATE FY_TALYS VECTOR WITH NUMBER OF COLUMNS DEPENDING ON GEF OPTION LMD OR LMD+ ---------------------#
        rows, columns = np.shape(FY)
        if columns == 7:                                    
            FY_TALYS = np.zeros((300,11),dtype=np.float32)                                                                          # columns = 7 for GEF lmd option. Add 4 for mean and std vectors.
        else:
            FY_TALYS = np.zeros((300,19),dtype=np.float32)                                                                          # columns = 11 for GEF lmd+ option. Add 8 for mean and std vectors.
        del rows
        #------------------------------------- DELETE MULTI-CHANCE FISSION EVENTS ----------------------------------------------#
        index_Multichance_fission = np.where(FY[:,1]+FY[:,3] != A_compound)                                                         # Check if mass number of light and heavy fragment matches compound nucleus mass number.
        if index_Multichance_fission[0].size > 0: 
            FY = np.delete(FY[:,:],index_Multichance_fission[0],axis=0)                                                             # Delete rows with multi-chance fission events.
        del index_Multichance_fission 
        #------------------------------- SWAP PLACES OF ELEMENTS IF AL > AH BEFORE PERFORMING AVERAGES -------------------------#
        index_Al_larger_than_Ah = np.where(FY[:,1] > FY[:,3])
        if index_Al_larger_than_Ah[0].size > 0:
            for row_number in index_Al_larger_than_Ah[0]:
                FY[row_number,0],FY[row_number,2],FY[row_number,1],FY[row_number,3] = \
                FY[row_number,2], FY[row_number,0],FY[row_number,3], FY[row_number,1]   
        del index_Al_larger_than_Ah  
        #------------------------------------ REMOVE EVENTS THAT ONLY OCCUR ONCE -----------------------------------------------#
        unique_yields ,index_unique, unique_counts = np.unique(FY[:,0:4], return_index=True, return_counts=True, axis=0) 
        index_event_occur_more_than_once = np.nonzero(unique_counts > 1) 
        index_to_pick = index_unique[index_event_occur_more_than_once[0]]                                                           # Pick out index of unique events that occur more than once.
        #----------------------------------- POPULATE [0] = Z1, [1] = A2, [2] = Z2, [3] = A2 -----------------------------------#
        FY_TALYS[slice(len(index_to_pick)),0] = FY[index_to_pick,0] # Z1
        FY_TALYS[slice(len(index_to_pick)),1] = FY[index_to_pick,1] # A1sci
        FY_TALYS[slice(len(index_to_pick)),2] = FY[index_to_pick,2] # Z2
        FY_TALYS[slice(len(index_to_pick)),3] = FY[index_to_pick,3] # A2sci
        #------------------------------------------------- POPULATE [4] = YIELDS -----------------------------------------------#
        ignored_events = np.count_nonzero(unique_counts == 1) 
        FY_TALYS[slice(len(index_to_pick)),4] = unique_counts[index_event_occur_more_than_once[0]]/(len(FY[:,0])-ignored_events)    # Divide by those events that are left to get the yield as a fraction.
        #---------------------- POPULATE [5] = TKE, [6] = TXE, [7] = Eexc1, [8] = Wl, [9] = Eexc2, [10] = Wh -------------------#
        for k, row in enumerate(unique_yields[index_event_occur_more_than_once]):
            ind = (row == FY[:,0:4]).all(axis=1).nonzero() 
            FY_TALYS[k,5]  = np.mean(FY[ind[0],4], dtype=np.float64)                    # [5]  = TKE.
            FY_TALYS[k,6]  = np.mean(FY[ind[0],5] + FY[ind[0],6], dtype=np.float64)     # [6]  = TXE.            
            FY_TALYS[k,7]  = np.mean(FY[ind[0],5], dtype=np.float64)                    # [7]  = Eexc1 (E* light fragment).
            FY_TALYS[k,8]  = np.std(FY[ind[0],5], dtype=np.float64, ddof=1)             # [8]  = Wl (std(E*), light fragment).
            FY_TALYS[k,9]  = np.mean(FY[ind[0],6], dtype=np.float64)                    # [9]  = Eexc2 (E* heavy fragment).
            FY_TALYS[k,10] = np.std(FY[ind[0],6], dtype=np.float64, ddof=1)             # [10] = Wh (std(E*), heavy fragment).
            #-------------------------- POPULATE neutron and gamma energies and std --------------------------------------------#
            if columns == 11:   # If simulation is with GEF "lmd+" option.
                FY_TALYS[k,11]  = np.mean(FY[ind[0],7], dtype=np.float64)               # [11] = Mean neutron energies of light fragments.
                FY_TALYS[k,12]  = np.std(FY[ind[0],7], dtype=np.float64, ddof=1)        # [12] = Std neutron energies of light fragments.
                FY_TALYS[k,13]  = np.mean(FY[ind[0],8], dtype=np.float64)               # [13] = Mean neutron energies of heavy fragments.
                FY_TALYS[k,14]  = np.std(FY[ind[0],8], dtype=np.float64, ddof=1)        # [14] = Std neutron energies of heavy fragments.
                FY_TALYS[k,15]  = np.mean(FY[ind[0],9], dtype=np.float64)               # [15] = Mean gamma energies of light fragments.
                FY_TALYS[k,16]  = np.std(FY[ind[0],9], dtype=np.float64, ddof=1)        # [16] = Std gamma energies of light fragments.
                FY_TALYS[k,17]  = np.mean(FY[ind[0],10], dtype=np.float64)              # [17] = Mean gamma energies of heavy fragments.
                FY_TALYS[k,18]  = np.std(FY[ind[0],10], dtype=np.float64, ddof=1)       # [18] = Std gamma energies of heavy fragments.
        del unique_yields,index_unique,unique_counts,index_event_occur_more_than_once,index_to_pick,num_of_events,k,row,ind,columns
        return FY_TALYS, ignored_events
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#

    def create_GEF_workingdir_and_inputfile(unique_param_ID,Z_target,A_compound,E_reaction,MC_runs,GEF_working_dir_path):
        """Create individual GEF output folder and individual input file.
        
        Parameters
        ----------
        unique_param_ID : `str`
            Used to create individual GEF output folder and input file.
        Z_target : `int`
            Atomic number of target chemical element. 
                
            Needed to create path to GEf ".lmd" and ".dat" files.
        A_compound : `int`
            Mass number of compound target nuclei. I.e A + 1.
                    
            Needed to create path to GEf ".lmd" and ".dat" files.
        E_reaction : `float`
            Kinetic energy (MeV) of incident neutron in fission event.
                    
            Needed to create path to GEF ".lmd" and ".dat" files.
        MC_runs : `int`
            Given in the GEF input file as a multiplication factor
            of the default number of Monte Carlo simulations. 
                    
            The GEF default value is 10^5 MC simulations. 
        GEF_working_dir_path : `str`
            Path to modified GEF output folder in McPUFF output folder
            
        Returns
        -------
        GEF_data_paths : `dict`
            GEF_data_paths dictionary with keys:
            
            ``"GEF_cwd_path"`` (`str`)
                Individual GEF output folder.
            ``"GEF_LMD_path"`` (`str`)
                Path to GEF ".lmd" file from simulation.
            ``"GEF_DAT_path"`` (`str`)
                Path to GEF ".dat" file from simulation.
            ``"GEF_DMP_EN_path"``
                Path to GEF "EN.dmp" file from simulation.
            
        See Also
        --------
        ``GEF_input.create_GEF_input_file_and_folder()``

        Notes
        -----
        GEF cannot handle a `float` with a zero decimal. It converts it 
        from `float` to `int`, causing McPUFF to crash as GEF crashes due
        to error in the incident neutron energy. To remedy this, this 
        function checks if the decimal equals zero. If it does, the 
        function converts the incident neutron energy to `int`, otherwise
        the function passes it as it is.
        
        Examples
        --------
        Function returns a dictionary of string corrensponding to the
        paths described in the "Returns" section.

        >>> GEF_paths = Reaction.create_GEF_workingdir_and_inputfile("TMC_5",
                             92, 236, 2.53e-8, 1e6, "/local/path/GEF/output/")
        [ dict [str, str, str] ]
        """
        #------------------------------------- CREATE UNIQUE FOLDER FOR INDIVIDUAL GEF SIMULATION ------------------------------#
        if not os.path.isdir(os.path.join(GEF_working_dir_path,f'GEF_{unique_param_ID}')):                                          # If folder does not exist, create it.
            GEF_cwd_path = os.path.join(GEF_working_dir_path,f'GEF_{unique_param_ID}')                                              # Declare path to folder.
            os.mkdir(GEF_cwd_path)
            os.path.join(GEF_working_dir_path,f'GEF_{unique_param_ID}',"")                                                          # Empty string at the end adds a'/'.
        else:
            GEF_cwd_path = os.path.join(GEF_working_dir_path,f'GEF_{unique_param_ID}',"")                                           # If folder does exist, declare path to folder.
        #-------------------------------- CREATE PATHS TO GEF LMD AND DAT FILES IN PUT FOLDER ----------------------------------#
        #---------------------------- GEF CANNOT HANDLE ZERO DECIMALS, IT CONVERTS THEM FROM FLOAT TO INT ----------------------#
        if str(E_reaction).endswith('.0'):                                                                                          # If e.g E_reaction = 2.0, it has to be int because GEF rounds off 2.0 to 2.
            GEF_LMD_path = os.path.join(GEF_cwd_path,'out',f'Z{str(Z_target)}_A{str(A_compound)}_n_E{int(E_reaction)}MeV.lmd')
        else:
            GEF_LMD_path = os.path.join(GEF_cwd_path,'out',f'Z{str(Z_target)}_A{str(A_compound)}_n_E{str(E_reaction)}MeV.lmd')      # Construct path based on how GEF names files.
        GEF_DAT_path = os.path.join(GEF_cwd_path,'out',f'GEF_{str(Z_target)}_{str(A_compound)}_n.dat')
        #------------------------------------- CREATE PATHS TO FILES IN DMP FOLDER ---------------------------------------------#
        GEF_DMP_EN_path = os.path.join(GEF_cwd_path,'dmp',f'Z{str(Z_target)}_A{str(A_compound)}_n_E{str(E_reaction)}MeV','EN.dmp')
        #------------------------------------------- CREATE GEF INPUT FILE -----------------------------------------------------#
        GEF_input.create_GEF_input_file_and_folder(GEF_cwd_path,Z_target,A_compound,E_reaction,MC_runs)
        #------------------------------------- COLLECT PATHS IN DICTIONARY AND RETURN ------------------------------------------#
        GEF_data_paths = {'GEF_cwd_path':GEF_cwd_path,'GEF_LMD_path':GEF_LMD_path,'GEF_DAT_path':GEF_DAT_path,'GEF_DMP_EN_path':GEF_DMP_EN_path}
        return GEF_data_paths
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#

    def perturbed_calculations_single_parameter(key,dict_unpert_param_name_val,Z_target,A_compound,E_reaction,runs_MC,num_of_rand,
                                                                        num_of_workers,With_TALYS_flag,distribution_flag,path_dict):  
        """Multi-thread simulations using "concurrent.futures" for the 
        "Single_Parameters" mode.
        
        This function creates separate threads for parallel simulations.
        It calls the "Modified_Parameter.create_perturbed_FY()" function
        which performs the individual simulations.

        Parameters
        ----------
        key : `str`
            Name of parameter that is to be perturbed during a simulation.
        dict_unpert_param_name_val : `dict` [`str`, `float`]
            Dictionary of parmeter names and default parameter values.
        Z_target : `int`
            Atomic number of target chemical element. 
        A_compound : `int`
            Mass number of compound target nuclei. I.e A + 1.
        E_reaction : `float`
            Kinetic energy (MeV) of incident neutron in fission event.
        MC_runs : `int`
            The number of Monte Carlo simulations to be performed in GEF. 
        num_of_rand : `int`
            The number of random numbers to produce, which is the same as
            the number of perturbed simulations to be performed.
        num_of_workers : `int`
            The number of CPU's to use for multi-threading.

            See the "See Also" section.
        With_TALYS_flag : `boolean`
            If `True`, includes TALYS evaporation simulation. 
        distribution_flag : `str`
            Specifies if random numbers are to be drawn from "uniform", 
            "normal" or "max-min" distributions.
        path_dict : `dict` [`str`,`str`,`str`,`str`,`str`,`str`,`str`]
            Dictionary with local paths needed for simulations.

        Returns
        -------
        mod_param_object : ``Modified_Parameter``
            Perturbed parameter object with attributes:

            ``param_name``
                GEF parameter name (`str`).
            ``unpert_param_val``
                Default GEF parameter value (`float`).
            ``list_of_Random_Parameter_Value_objects``
                List of objects containing data for a perturbation of a
                parameter using a specific random number. (`list`).
            ``list_of_rand_num``
                List of all random numbers used for perturbations. (`list`).
       
        See Also
        --------
        Headline "CREATE PERTURBED FY IN 'SINGLE_PARAMETERS' MODE" in 
        Reaction.__init__.
        Class ``Modified_Parameter``.
        Class ``Random_Parameter_value``.
        ``Modified_Parameter.create_perturbed_FY()``.

        Notes
        -----
        The function creates a ``Modified_Parameter`` object for each
        thread and passes it as an argument to the  
        ``Modified_Parameter.create_perturbed_FY()`` function. That 
        function creates a ``Random_Parameter_value`` object for each
        randomly perturbation which, when completed, is stored in the 
        ``Modified_Parameter`` object. Due to the `"ALL_COMPLETED"` 
        argument, the `concurrent.futures.ThreadPoolExecutor` then waits
        until all random perturbations for a parameter are completed. This
        procedure is necessary to ensure that a ``Random_Parameter_value``
        object is not stored in the wrong ``Modified_Parameter`` object 
        when multiple threads are running asynchronously.   

        Examples
        --------
        >>> obj = Reaction.perturbed_calculations_single_parameter("_Delta_S0",
                dict_unpert_param_name_val,92,236,2.53e-8,1e6,500,40,
                False,"normal",path_dict):
        [obj (Modified_Parameter object)]  
        """ 
        #-----------------------------------------------------------------------------------------------------------------------#    
        param_name = key
        """Name of GEF parameter to perturb."""
        unpert_param_value = dict_unpert_param_name_val[key]
        """Default GEF parameter value."""
        mod_param_object = Modified_Parameter(param_name,unpert_param_value,num_of_rand,distribution_flag)
        """Perturbed parameter object"""  
        simultaneous_threads_per_param = math.floor((os.cpu_count()-num_of_workers)/2)  
        """Determines number of cpu's that are used for each parameter. """
        #----------------------------------------START MULTI-THREAD SIMULATIONS-------------------------------------------------# 
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=simultaneous_threads_per_param) as perturbed_single_executor:  
                future_rand_param_val = {perturbed_single_executor.submit(Modified_Parameter.create_perturbed_FY,param_name,
                                            unpert_param_value,rand_num,n,Z_target,A_compound,E_reaction,
                                                runs_MC,With_TALYS_flag,distribution_flag,path_dict): rand_num for n, rand_num in enumerate(mod_param_object.list_of_rand_num)}
            future_done,future_not_done = concurrent.futures.wait(future_rand_param_val,return_when="ALL_COMPLETED")
            for rand_param_obj in future_done:
                mod_param_object.list_of_Random_Parameter_Value_objects.append(rand_param_obj.result())
            del param_name,unpert_param_value,future_rand_param_val,future_done,future_not_done,perturbed_single_executor
            return mod_param_object
        except Exception as e:
            print('Encountered a problem in perturbed_calculations_single_parameter()\n')
            print(e)
            sys.exit(e)
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#
            
    def print_TALYS_ff_files(pth_TALYS_folder,FY_TALYS_format,Z_target,A_compound,E_reaction,unique_thread_ID):
        """Create a perturbed fission fragment yield file (".ff") and
        place it in the GEF ".ff" file library.
        
        Parameters
        ----------
        pth_TALYS_folder : `str`
            Path to modified GEF library of ".ff" files.
        With GEF "lmd" option:
        FY_TALYS_format : `numpy.ndarray` (300,11)
            Array with perturbed GEF fission fragment yield data assembled
            in GEF ".ff" file format for TALYS to use for simulations [1]_
            dtype = (`numpy.float32`).
        With GEF "lmd+" option:
        FY_TALYS_format : `numpy.ndarray` (300,19)
            Array with perturbed GEF fission fragment yield data assembled
            in GEF ".ff" file format for TALYS to use for simulations [1]_
            dtype = (`numpy.float32`).
                          
            The GEF "lmd+" option includes neutron and gamma 
            energies per fragment. See the `See Also` section.
        Z_target : `int`
            Atomic number of target chemical element used to name file in 
            correct format. 
        A_compound : `int`
            Mass number of compound target nuclei, i.e (A + 1), use to
              name file in correct format.
        E_reaction : `float`
            Kinetic energy (MeV) of incident neutron in fission event, 
            used to name file in correct format.
        unique_thread_ID : `str`
            Used to name files and enable TALYS to differentiate between
            perturbed ".ff" files in GEF library during multi-threading. 

        See Also
        --------
        ``Reaction.FY_results()``

        Notes
        -----
        The unperturbed ".ff" files in the GEF library are made with 
        simulations with integer values for the kinetic energy of the 
        incident neutron (E_kinetic). When TALYS performs a simulation of 
        the evaporation process with a non-integer value for E_kinetic, it 
        interpolates using data from files with adjacent integer values. 
        In order to not interfere with this procedure during perturbed 
        simulations, McPUFF creates two identical ".ff" files with 
        adjacent integer values for E_kinetic and places them both in the 
        TALYS ".ff" file library. They only difference between them are 
        the values for E_kinetic in the file-name. They don't actually 
        have different energies. The user can choose these values to suit
        specfic simulations by entering them in the list variable 
        ``energy_val``. 
        
        References
        ----------
        .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
        GEF and its influence on the nuclear evaporation in TALYS." 
        Technical report UPTEC F 23065, Uppsala University, 2023, accessed
        1 Janauary 2024, 
        url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>.

        Examples
        --------
        Function creates two files with different names for the energies at 
        specified location.

        >>> Reaction.print_TALYS_ff_files("/pth_TALYS_folder/",
                                   FY_TALYS_format,92,236,2.53e-8,"TMC_5")
        ["U236_6.00e+00MeV_gef_tmc_5.ff", "U236_7.00e+00MeV_gef_tmc_5.ff"]
        """
        number_of_FY = len(np.nonzero(FY_TALYS_format[:,0])[0])
        """Number of unique GEF fission events."""
        #-------------------------------------- SET VALUE FOR KINETIC ENERGY IN ".ff" FILE NAME --------------------------------# 
        if str(A_compound) == '239':                                                # Different elements may have different energy values for the fission barrier. This `if` statement can be extended to incorporate more elements. 
            energy_val = ['6','7']   
            """Kinetic energy values of the incident neutron.
            
            Create file names with different values that TALYS can 
            interpolate between."""   
        else: 
            energy_val = ['6','7']
        
        elements = {'92':'U','94':'Pu'}                                                                                                             # This dictionary can be extended to incorporate more elements.
        """Dictionary of atomic number and chemical symbol of elements
        used in simulations."""
        structure = 'Zl  Al   Zh  Ah   Yield       TKE[MeV]    TXE[MeV]'\
                    '    El[MeV]     Wl[MeV]     Eh[MeV]     Wh[MeV]' 
        """Required string format of ".ff" files in GEF library.""" 
        #============================ CODE FOR CREATING CUSTOM E-FILE IN THE GEF ".ff" FILE LIBRARY. ===========================#
        #path_E_file = os.path.join(pth_TALYS_folder,f'{elements[str(Z_target)]}{str(A_compound)}_gef.E')
        #if str(A_compound) == '239':
        #    with open(path_E_file,'w') as E_file:
        #        E_file.write('5.00e+00\n6.00e+00')
        #else:
        #    with open(path_E_file,'w') as E_file:
        #        E_file.write('6.00e+00\n7.00e+00')
        #------------------------------------ CREATE ".FF" FILE WITH FISSION FRAGMENT YIELDS------------------------------------#    
        for i in range(len(energy_val)):                                                                                                            # Create ".ff" files according to E_kinetic values in `energy_val` list.
            ff_file_name = f'{elements[str(Z_target)]}{str(A_compound)}_{str(energy_val[i])}.00e+00MeV_gef_{str(unique_thread_ID).lower()}.ff'                                 
            with open(os.path.join(pth_TALYS_folder,ff_file_name), 'w') as TALYS_format:
                TALYS_format.write('# Z        =   {:>3}\n'.format(Z_target) + '# A        =   {:>3}\n'.format(A_compound) +\
                '# Ex (MeV) =   {:>3.2e}\n'.format(float(E_reaction)) +\
                '# Ntotal   =   {:>3}\n'.format(number_of_FY) + f'# {structure}\n')
                for row in FY_TALYS_format[:,:]:
                    if row[0] !=0:
                        TALYS_format.write(f'{row[0]:>4.0f} {row[1]:>3.0f} {row[2]:>4.0f} {row[3]:>4.0f}  {row[4]:.4e}  {row[5]:.4e}  {row[6]:.4e}  {row[7]:.4e}  {row[8]:.4e}  {row[9]:.4e}  {row[10]:.4e}\n')           
        del number_of_FY,unique_thread_ID,energy_val,elements,structure,TALYS_format,ff_file_name,row,Z_target,A_compound,E_reaction
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#
        
    def read_and_clear_GEF_results(path_GEF_results,path_GEF_DAT,GEF_DMP_EN_path): 
        """Store selected simulation results from different GEF output 
        files and then delete output data to save disc space and memory.
        
        Parameters
        ----------
        path_GEF_results : `str`
            Path to GEF output in Output_McPUFF folder.
        path_GEF_DAT : `str`
            Local path to GEF output ".dat" file.
        GEF_DMP_EN_path: `str`
            Local path to GEF output ".dmp" file.

        Returns
        -------
        GEF_data : `dict`
            Dictionary of GEF results. 

            The contents are subject to change. The user can add or remove
            content of dictionary. Data can be any data type. More 
            information about GEF output data can be found in the GEF 
            documentation in the "See Also" section.

        See Also
        --------
        GEF 2023/1.1 "ReadMe" file. url:
        <https://www.khschmidts-nuclear-web.eu/GEF_code/GEF-2023-1-1/Standalone/Readme.txt>

        Python built-in function "str.index()" 
        url: <https://docs.python.org/3/library/stdtypes.html?highlight=string%20index>

        ``Reaction.create_GEF_workingdir_and_inputfile()``

        ``Reaction.delete_GEF_result_folder()``

        Rough estimates of GEF parameter stability intervals: 
        url:<https://github.com/UPTEC-F-23065/McPUFF/blob/12b9f864240fe133ef727955c9e7715a7e0e20a9/Rough_limits_GEF_parameter_values.txt>

        Notes
        -----
        The GEF output data is stored here as "raw" data. The format of 
        the GEF output data can sometimes be cumbersome and in order to
        keep McPUFF fast, as little formatting as possible is done to the
        output data here. Data necessary for later formatting is saved in 
        the GEF output data dictionary, e.g. "GEF_data['A_TKE_pre_ranges']"
        and "GEF_data['E_N_A_ranges']" that holds the mass numbers that 
        correlate to the data.

        The list variable ``headlines_wanted_data`` holds rows of text
        from the GEF output files used to find line index using the 
        Python `str.index()` function. For a single value, one line index
        is enough. For multiple values, using a start index and a stop 
        index is the most suitable way. The commented section below 
        provides examples. The user can add/remove data to be saved from
        a simulation by adding/removing line indexes from 
        ``headlines_wanted_data`` and adding/removing code sections 
        reading the data. Note that the Python newline character "\\n" 
        has to be included in order for the `index()` function to find a
        match.

        If data from other GEF output files are to be read, the path to 
        the file has to be passed to this function. The path must be added
        in the function ``Reaction.create_GEF_workingdir_and_inputfile()``.

        This function calls ``Reaction.delete_GEF_result_folder()`` in 
        order to delete GEF simulation data after selected data is stored.
        This function can be turned off by commenting out the code line.
        It is useful for inspecting GEF output files. Be warned that the 
        GEF output data can become quite large (GB) when performing many 
        perturbed simulations.

        For some values of the GEF parameters, the GEF software produces 
        non-physical fission event results (e.g violate energy principle).
        GEF marks these events by adding a "%" in the output file, 
        causing McPUFF to crash (non-numerical data). In order to avoid
        crashing, McPUFF removes the line containing the "%" symbol and 
        stores the rest of the data. Most likely, these data sets are 
        the ones that do not result in a TALYS simulation. Hence, 
        proposedly, these non-physical simulations can be identified by 
        the lack of TALYS results in the `pickle` file. A rough estimate 
        of the parameter value stability intervals can be found in the 
        McPUFF GitHub repository (See the "See Also" section).

        Note that the GEF data for prompt neutron multiplicity as a 
        function of mass A, has two columns for mean values. One for the 
        compound nucleus (CN) and one for the fragments. 

        Examples
        --------
        >>> GEF_data_dict = Reaction.read_and_clear_GEF_results("/local/path/GEF/output/","/local/path/GEF/.dat/file")
        [GEF_data_dict dict ]
        """
        #-----------------------------------------------------------------------------------------------------------------------#
                                #==================================================================================================================#
                                #                Explanation for index of lines in ``headlines_wanted_data`` used to read GEF data.          
                                #                ----------------------------------------------------------------------------------
                                #   Example:
                                #   -------
                                #   index_gamma_multi_func_of_A = [dat_data.index(headlines_wanted_data[0])+4, dat_data.index(headlines_wanted_data[1])-4],
                                #   means that the index() function finds the string matching headlines_wanted_data[0] in the GEF file and starts reading from [0]+4 lines.
                                #   The index() function finds the string matching headlines_wanted_data[1] in the GEF file and stops reading at [1]-4.
                                #   "col: 0-22" means that data from columns 0 to 22 in the file are to be read. "col: 3" means a single value is read.
                                #==================================================================================================================#
                                #                                               DAT FILE DATA
                                #------------------------------------------------------------------------------------------------------------------#
                                # [0]->[1]  = Gamma multiplicity as function of A :             From: [0]+4 to [1]-4,   col: 0-22
                                # [2]       = Mean gamma multiplicity:                          From: [2]-1,            col: 3 single value
                                # [3]       = Mean gamma energy :                               From: [3]-2,            col: 3 single value
                                # [4]-[5]   = Prompt-neutron multiplicity as func of pre A:     From: [4]+2 to [5]-1,   col: 0-16 
                                # [5]-[6]   = Prompt-neutron multiplicity as func of post A:    From: [5]+2 to [6]-4,   col: 0-16 
                                # [7]       = Mean prompt neutron multiplicity from fragments:  From: [7]-9,            col: 4 single value
                                # [7]       = Standard deviation prompt neutron from fragments: From: [7]-8,            col: 3 single value
                                # [8]       = Mean neutron energy:                              From: [8]+3,            col: 4 single value
                                # [9]       = Mean value TKE-pre:                               From: [9]-5,            col: 4 single value
                                # [9]       = Mean value TKE-post:                              From: [9]-5,            col: 8 single value
                                # [10]-[11] = A-TKE spectrum (pre-neutron):                     From: [10]+8 to [11]-2, col: 0-50 last col<50 elements
                                # [11]-[12] = A-TKE spectrum (post-neutron):                    From: [11]+8 to [12],   col: 0-50 last col<50 elements 
                                # [13]      = Mean value Q-bar                                  From: [13]-4,           col: 4 single value
                                # [14]      = Mean value TXE                                    From: [14]-4,           col: 4 single value
                                #------------------------------------------------------------------------------------------------------------------#
                                #                                               DMP FILE DATA
                                #------------------------------------------------------------------------------------------------------------------#
                                # [15]      = Mass numbers of fragments (lowest and highest):   From: [15]+6            col: 3 and 5. 
                                # [15]-[16] = Mean neutron energy as function of A:             From: [15]+7 to [16]-3  col: csv lines of different length. 
                                #==================================================================================================================#
        headlines_wanted_data = ['--- Mass-dependent gamma multiplicity (from fragments) ---\n',                            # [0]
                                '--- Total gamma-multiplicity distribution (emission from fragments) ---\n',                # [1]
                                '      </Gamma_multiplicity>\n',                                                            # [2]
                                '      </E_gammas>\n',                                                                      # [3]
                                ' Apre  Nmean  Nmean          Multiplicity distribution, unnormalized (0 to 16)\n',         # [4]
                                ' Apost  Nmean  Nmean         Multiplicity distribution, unnormalized (0 to 16)\n',         # [5]
                                '--- Multiplicity distribution of prompt neutrons (only fission events considered ---\n',   # [6]
                                '--- Multiplicity distribution of prompt neutrons (light fragment) ---\n',                  # [7]    
                                '      </Nspectrum>\n',                                                                     # [8]
                                '--- TXE spectrum (bins with zero content suppressed) ---\n',                               # [9]
                                '--- A-TKE spectrum (pre-neutron)---\n',                                                    # [10]             
                                '--- A-TKE spectrum (post-neutron)---\n',                                                   # [11]
                                '    </A_TKE>\n',                                                                           # [12]
                                '--- TKE spectrum (pre- and post-neutron) (bins with zero content suppressed) ---\n',       # [13]
                                '--- A-Ekin spectrum (pre-neutron)---\n',                                                   # [14]
                                'S: TITLE(Mean neutron energy over pre-neutron mass (from fragments in fragment frame))\n', # [15]
                                'S: ANALYZER(ENApostfs)\n']                                                                 # [16]
        """List of rows in file used to find index of lines to read."""
        GEF_data = {}
        """Dictionary used to store perturbed GEF simulation data."""
        if os.path.isfile(path_GEF_DAT):
            try:
                #-------------------------------------- READ ".dat"-FILE AND CREATE INDEXES  -----------------------------------#
                with open(path_GEF_DAT,'r') as GEF_dat_file:
                    dat_data = GEF_dat_file.readlines()
                index_gamma_multi_func_of_A = [dat_data.index(headlines_wanted_data[0])+4, dat_data.index(headlines_wanted_data[1])-4] 
                #----------------------------------  CHECK FOR AND REMOVE ERRONOUS DATA IN FILE  -------------------------------#       
                #-------------- CREATE INDEXES FOR PROMPT GAMMA MULTIPLICITY AS A FUNCTION OF A (PRE-NEUTRON)-------------------#
                for line in dat_data[index_gamma_multi_func_of_A[0]:index_gamma_multi_func_of_A[1]]:                    
                    if '%' in line:                                                                     # Check if there are non-numerical entries.
                        del dat_data[dat_data.index(line)]                                              # If so, remove line and update indexes.
                        index_gamma_multi_func_of_A = [dat_data.index(headlines_wanted_data[0])+4, dat_data.index(headlines_wanted_data[1])-4] 
                #-------------- CREATE INDEXES FOR PROMPT NEUTRON MULTIPLICITY AS A FUNCTION OF A (PRE-NEUTRON) ----------------#
                index_prompt_neutron_multi_func_of_A_pre = [dat_data.index(headlines_wanted_data[4])+2, dat_data.index(headlines_wanted_data[5])-1] 
                for line in dat_data[index_prompt_neutron_multi_func_of_A_pre[0]:index_prompt_neutron_multi_func_of_A_pre[1]]:    
                    if '%' in line:                                                                     # Check if there are non-numerical entries.
                        del dat_data[dat_data.index(line)]                                              # If so, remove line and update indexes.
                        index_prompt_neutron_multi_func_of_A_pre = [dat_data.index(headlines_wanted_data[4])+2, dat_data.index(headlines_wanted_data[5])-1]
                #-------------- CREATE INDEXES FOR PROMPT NEUTRON MULTIPLICITY AS A FUNCTION OF A (POST-NEUTRON) ---------------#
                #index_prompt_neutron_multi_func_of_A_post = [dat_data.index(headlines_wanted_data[5])+2, dat_data.index(headlines_wanted_data[6])-4]
                #for line in dat_data[index_prompt_neutron_multi_func_of_A_post[0]:index_prompt_neutron_multi_func_of_A_post[1]]:    
                #    if '%' in line:                                                                    # Check if there are non-numerical entries.
                #        del dat_data[dat_data.index(line)]                                             # If so, remove line and update indexes.
                #        index_prompt_neutron_multi_func_of_A_post = [dat_data.index(headlines_wanted_data[5])+2, dat_data.index(headlines_wanted_data[6])-4]
                #----------------------------------------------- STORE GEF GAMMA DATA  -----------------------------------------#
                #---------------------------------------------------------------------------------------------------------------#
                #   [0]-[1]: Prompt gamma multiplicity as a function on mass A.   [2]: Mean gamma multiplicity as a function on mass A.
                #   [3]: Mean gamma energy.
                #---------------------------------------------------------------------------------------------------------------#
                GEF_data['gamma_multi_func_of_A'] =  [np.loadtxt(dat_data,dtype=np.float32, skiprows=index_gamma_multi_func_of_A[0], 
                                                                                                    max_rows=(index_gamma_multi_func_of_A[1]-index_gamma_multi_func_of_A[0]))]
                GEF_data['mean_gamma_multi'] = [dat_data[dat_data.index(headlines_wanted_data[2])-1].split()[3]]   
                GEF_data['Mean gamma energy'] = [dat_data[dat_data.index(headlines_wanted_data[3])-2].split()[3]]  
                #--------------------------------------------- STORE GEF NEUTRON DATA  -----------------------------------------#
                #---------------------------------------------------------------------------------------------------------------#
                #   [4]-[5]: Prompt neutron multiplicity as function of A (pre-neutron).    [5]-[6]: Prompt neutron multiplicity as function of A (post-neutron).   
                #   [7]: mean prompt neutron multiplicity as function of A from fragments.  [7]: Standard deviation of prompt n from fragments.
                #   [8]: Mean neutron energy. 
                #---------------------------------------------------------------------------------------------------------------#
                GEF_data['prompt_neutron_multi_func_of_A_pre'] = [np.loadtxt(dat_data,dtype=np.float32,skiprows=index_prompt_neutron_multi_func_of_A_pre[0],
                                                                                        max_rows=(index_prompt_neutron_multi_func_of_A_pre[1]-index_prompt_neutron_multi_func_of_A_pre[0]))]
                #GEF_data['prompt_neutron_multi_func_of_A_post'] = [np.loadtxt(dat_data,dtype=np.float32,skiprows=index_prompt_neutron_multi_func_of_A_post[0],
                #                                                                        max_rows=(index_prompt_neutron_multi_func_of_A_post[1]-index_prompt_neutron_multi_func_of_A_post[0]))]
                GEF_data['mean_prompt_neutron_multi_from_frag'] = [dat_data[dat_data.index(headlines_wanted_data[7])-9].split()[4]]   
                GEF_data['st_dev_prompt_neutron_multi_from_frag'] = [dat_data[dat_data.index(headlines_wanted_data[7])-8].split()[3]] 
                GEF_data['mean_neutron_E'] = [dat_data[dat_data.index(headlines_wanted_data[8])+3].split()[4]]  
                #--------------------------------------------- A VS TKE PRENEUTRON DATA ----------------------------------------#
                #---------------------------------------------------------------------------------------------------------------#
                #   [9]: Mean TKE pre-neutron.   [10]-[11]: A_TKE value ranges and Average TKE pre-neutron as a function of A. 
                #---------------------------------------------------------------------------------------------------------------#
                GEF_data['mean_value_TKE_pre'] = [dat_data[dat_data.index(headlines_wanted_data[9])-5].split()[4]]  
                index_A_TKE_preneutron = [dat_data.index(headlines_wanted_data[10])+8, dat_data.index(headlines_wanted_data[11])-2]     
                Astart = dat_data[index_A_TKE_preneutron[0]-5].split()[4]   # Index from '--- A-TKE spectrum (pre-neutron)---\n' +8 -5. (Reuse index for A-TKE spectrum (pre-neutron)).
                Astop  = dat_data[index_A_TKE_preneutron[0]-5].split()[6]
                Estart = dat_data[index_A_TKE_preneutron[0]-5].split()[11]
                Estop  = dat_data[index_A_TKE_preneutron[0]-5].split()[13]
                GEF_data['A_TKE_pre_ranges'] = [f'A_TKE_pre_A_{Astart}_to_{Astop}_E_{Estart}_to_{Estop}']
                GEF_data['A_TKE_pre_neutron'] = [dat_data[index_A_TKE_preneutron[0]:index_A_TKE_preneutron[1]]]
                del Astart, Astop, Estart, Estop
                #----------------------------------- A VS TKE POST-NEUTRON DATA ------------------------------------------------#
                #---------------------------------------------------------------------------------------------------------------#
                #    [9]: Mean TKE post-neutron.     [11]-[12]: A_TKE value ranges and Average TKE post-neutron as a function of A. 
                #---------------------------------------------------------------------------------------------------------------#
                #GEF_data['mean_value_TKE_post'] = [dat_data[dat_data.index(headlines_wanted_data[9])-5].split()[8]]  
                #index_A_TKE_postneutron = [dat_data.index(headlines_wanted_data[11])+8, dat_data.index(headlines_wanted_data[12])]
                #Astart = dat_data[index_A_TKE_postneutron[0]-5].split()[4]  #Index from '--- A-TKE spectrum (post-neutron)---\n' +8 -5. Reuse index for A-TKE spectrum (post-neutron)).
                #Astop  = dat_data[index_A_TKE_postneutron[0]-5].split()[6]
                #Estart = dat_data[index_A_TKE_postneutron[0]-5].split()[11]
                #Estop  = dat_data[index_A_TKE_postneutron[0]-5].split()[13]
                #GEF_data['A_TKE_post_ranges'] = [f'A_TKE_post_A_{Astart}_to_{Astop}_E_{Estart}_to_{Estop}']
                #GEF_data['A_TKE_post_neutron'] = [dat_data[index_A_TKE_postneutron[0]:index_A_TKE_postneutron[1]]] 
                #--------------------------------------------- MEAN VALUE Q-BAR ------------------------------------------------#
                #---------------------------------------------------------------------------------------------------------------#
                #   [13]: Mean value Q-bar (MeV)
                #---------------------------------------------------------------------------------------------------------------#
                GEF_data['mean_value_Q_bar'] = [dat_data[dat_data.index(headlines_wanted_data[13])-4].split()[4]] 
                #--------------------------------------------- MEAN VALUE TXE ------------------------------------------------#
                #---------------------------------------------------------------------------------------------------------------#
                #   [14]: Mean value TXE (MeV)
                #---------------------------------------------------------------------------------------------------------------#
                GEF_data['mean_value_TXE'] = [dat_data[dat_data.index(headlines_wanted_data[14])-4].split()[4]]             
                #---------------------------------------------------------------------------------------------------------------#  
            except Exception as e:
                sys.exit(e)
            del GEF_dat_file,dat_data,index_gamma_multi_func_of_A,index_prompt_neutron_multi_func_of_A_pre,
            index_A_TKE_preneutron# ,Astart,Astop,Estart,Estop,index_A_TKE_postneutron,index_prompt_neutron_multi_func_of_A_post
            #-------------------------------------- READ ".dmp"-FILE AND CREATE INDEXES  -----------------------------------#
            try:
                with open(GEF_DMP_EN_path,'r') as GEF_dmp_file:
                    dmp_data = GEF_dmp_file.readlines()
                index_E_N_func_of_A = [dmp_data.index(headlines_wanted_data[15])+7, dmp_data.index(headlines_wanted_data[16])-3] 
                #---------------------------------------------------------------------------------------------------------------#
                #   [15]: E_N_A_ranges  (Mass number ranges for E(n) as a function af A).
                #   [15]-[16]: E_N_vs_A (Neutron energy values as a function of A).
                #---------------------------------------------------------------------------------------------------------------#
                Astart = dmp_data[dmp_data.index(headlines_wanted_data[15])+6].split()[3]
                Astop = dmp_data[dmp_data.index(headlines_wanted_data[15])+6].split()[5] 
                GEF_data['E_N_A_ranges'] = f'A_from: {Astart} to: {Astop}, increment by 1.'
                GEF_data['E_N_vs_A'] = dmp_data[index_E_N_func_of_A[0]:index_E_N_func_of_A[1]]
                del headlines_wanted_data,GEF_dmp_file,dmp_data,index_E_N_func_of_A,Astart,Astop
            except Exception as e:
                sys.exit(e)
            #--------------------------------------------------------------------------------------------------------------------------------------------------#
            #   Function for deleting GEF simulation files and folders after a simulation is finished. The data from the files are stored in the McPUFF output. #
            #---------------------------------------------------------------------------------------------------------------------------------------------------#
            Reaction.delete_GEF_result_folder(path_GEF_results)     # Comment out this line in order to study the complete GEF simulation output.      
            #-------------------------------------------------------------------------------------------------------------------#
            return GEF_data
        else:
            print('No such file can be found')
            sys.exit()
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#
            
    def read_and_clear_TALYS_results(pth_TALYS_CWD,pth_TALYS_FY_folder,unique_thread_ID): 
        """Store selected simulation results from different TALYS output 
        files and then delete output data to save disc space and memory.
        
        Parameters
        ----------
        pth_TALYS_CWD : `str`
            Local path to individual TALYS simulation output folder.
        pth_TALYS_FY_folder : `str`
            Local path to modified GEF library of ".ff" files.
        unique_thread_ID : `str`
            Used to name files and enable TALYS to differentiate between
            perturbed ".ff" files in GEF library during multi-threading. 

        Returns
        -------
        TALYS_data : `dict`
            Dictionary of TALYS results. 

            The contents are subject to change. The user can add or remove
            content of dictionary. Data can be any data type. More 
            information about the TALYS output data can be found in the 
            TALYS documentation in the "See Also" section.

        See Also
        --------
        "TALYS-1.96/2.0. Simulation of nuclear reactions". 
        url:<https://www-nds.iaea.org/talys/tutorials/talys_v1.96.pdf>

        `numpy.loadtxt()`. 
        url: <https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html>

        ``Reaction.delete_TALYS_result_files()``

        ``Reaction.delete_TALYS_ff_files()``

        Notes
        -----
        The TALYS output files are identified by the start of the file
        name string. The user can add/remove data to be stored by 
        following the procedure for the existing data storage below.

        The TALYS output data for the prompt gamma multiplicity as a 
        function of mass A, is given "per fragment". Hence, the values 
        should be multiplied by two to get the total prompt gamma 
        multiplicity.  

        This function calls ``Reaction.delete_TALYS_result_files()`` and
        ``Reaction.delete_TALYS_ff_files()`` in order to delete the TALYS
        output files and perturbed ".ff" file in the GEF library after 
        each simulation is completed. To retain the output data for 
        analysis, these functions can be commented out. Be warned that 
        the TALYS output data can become quite large (GB) when performing 
        many perturbed simulations. 

        Examples
        --------
        >>> Reaction.read_and_clear_TALYS_results("/local/path/TALYS/output/folder","/local/path/GEF/ff-file/library","TMC_5")
        [ dict ]
        """
        #-----------------------------------------------------------------------------------------------------------------------#
        data_to_read = ['EavA',         # [0]: Average emission energy per A.
                        'nugA',         # [1]: Mean number and multiplicity of prompt gamma as a function of A.
                        'nunA',         # [2]: Mean number and multiplicity of prompt neutrons as a function of A.
                        'pfgs',         # [3]: Average energy of PFGS and PFGS (prompt fission gamma spectrum)
                        'pfns',         # [4]: Average energy of PFNS and PFNS (prompt fission neutron spectrum)
                        'Pnug',         # [5]: Mean number prompt gammas
                        'Pnun',         # [6]: Mean number prompt neutrons and distributions of number of prompt neutrons
                        'yieldA']       # [7]: Fission products (yield after evaporation) and Fission fragments (data from GEF run) as a function of A
        """List of strings to match when storing data."""
        TALYS_data = {}   
        """Dictionary for storing TALYS output data from simulations."""                  
        try:
            for file_name in os.scandir(pth_TALYS_CWD):
                #---------------------------------------------------------------------------------------------------------------#
                #  data_to_read[0] = EavA: Average emission energy per A.                                                       #
                #---------------------------------------------------------------------------------------------------------------#
                if file_name.is_file() and file_name.name.startswith(data_to_read[0]):
                    path_file = os.path.join(pth_TALYS_CWD,file_name.name)
                    with open(path_file,'r') as file:   
                        data = file.readlines()
                    TALYS_data[f'{data_to_read[0]}_Avg_emission_E_for_n_g_func_of_A'] = [np.loadtxt(path_file,dtype=np.float32)]
                    del path_file,file,data
                #--------------------------------------------------------------------------------------------------------------------------------#
                # data_to_read[1] = nugA: Mean value of prompt gamma multiplicity and Average prompt gamma multiplicity as a function on mass A. #
                #--------------------------------------------------------------------------------------------------------------------------------#
                if file_name.is_file() and file_name.name.startswith(data_to_read[1]):
                    path_file = os.path.join(pth_TALYS_CWD,file_name.name)
                    with open(path_file,'r') as file:   
                        data = file.readlines()
                    TALYS_data[f'{data_to_read[1]}_Mean_value_(nubar-prompt)'] = [data[2].split()[5]]
                    TALYS_data[f'{data_to_read[1]}_Avg_prompt_gamma_multiplicity_func_of_A'] = [np.loadtxt(path_file,dtype=np.float32)]
                    del path_file,file,data
                #-------------------------------------------------------------------------------------------------------------------------------------#
                #  data_to_read[2] = nunA: Mean value of prompt neutron multiplicity and Average prompt neutron multiplicity as a function on mass A. # 
                #-------------------------------------------------------------------------------------------------------------------------------------#
                elif file_name.is_file() and file_name.name.startswith(data_to_read[2]):
                    path_file = os.path.join(pth_TALYS_CWD,file_name.name)
                    with open(path_file,'r') as file:   
                        data = file.readlines()
                    TALYS_data[f'{data_to_read[2]}_Mean_value_(nubar_prompt)'] = [data[2].split()[5]]
                    TALYS_data[f'{data_to_read[2]}_Avg_prompt_neutron_multiplicity_func_of_A'] = [np.loadtxt(path_file,dtype=np.float32)]
                    del path_file,file,data
                #---------------------------------------------------------------------------------------------------------------#
                #  data_to_read[3] = pfgs: Average energy of PFGS and PFGS (prompt fission gamma spectrum).                     #
                #---------------------------------------------------------------------------------------------------------------#
                elif file_name.is_file() and file_name.name.startswith(data_to_read[3]):
                    path_file = os.path.join(pth_TALYS_CWD,file_name.name)
                    with open(path_file,'r') as file:   
                        data = file.readlines()
                    TALYS_data[f'{data_to_read[3]}_E_average_gamma_MeV'] = [data[3].split()[3]]
                    TALYS_data[f'{data_to_read[3]}'] = [np.loadtxt(path_file,dtype=np.float32,usecols=(0,1))] 
                    del path_file,file,data
                #---------------------------------------------------------------------------------------------------------------#
                # data_to_read[4] = pfns: Average energy of PFNS and PFNS (prompt fission neutron spectrum).                    #
                #---------------------------------------------------------------------------------------------------------------#
                elif file_name.is_file() and file_name.name.startswith(data_to_read[4]):
                    path_file = os.path.join(pth_TALYS_CWD,file_name.name)
                    with open(path_file,'r') as file:   
                        data = file.readlines()
                    TALYS_data[f'{data_to_read[4]}_E_average_neutron_MeV'] = [data[3].split()[3]]
                    TALYS_data[f'{data_to_read[4]}'] = [np.loadtxt(path_file,dtype=np.float32,usecols=(0,1,2))]
                    del path_file,file,data
                #---------------------------------------------------------------------------------------------------------------#
                # data_to_read[5] = Pnug: Mean number prompt gammas.                                                            #
                #---------------------------------------------------------------------------------------------------------------#
                elif file_name.is_file() and file_name.name.startswith(data_to_read[5]):
                    path_file = os.path.join(pth_TALYS_CWD,file_name.name)
                    with open(path_file,'r') as file:   
                        data = file.readlines()
                    TALYS_data[f'{data_to_read[5]}_Mean_value_(nubar-prompt)'] = [data[2].split()[5]]
                    del path_file,file,data
                #---------------------------------------------------------------------------------------------------------------#
                # data_to_read[6] = Pnun: Mean number prompt neutrons, Distributions of number of prompt neutrons.              #
                #---------------------------------------------------------------------------------------------------------------#
                elif file_name.is_file() and file_name.name.startswith(data_to_read[6]):
                    path_file = os.path.join(pth_TALYS_CWD,file_name.name)
                    with open(path_file,'r') as file:   
                        data = file.readlines()
                    TALYS_data[f'{data_to_read[6]}_Mean_value_(nubar_prompt)'] = [data[2].split()[5]]
                    TALYS_data[f'{data_to_read[6]}'] = [np.loadtxt(path_file,dtype=np.float32)]
                    del path_file,file,data
                #---------------------------------------------------------------------------------------------------------------#
                # data_to_read[7] = yieldA: Number of fission fragments.                                                        #
                #   Col 1: Mass number A.                                                                                       #
                #   Col 2: Fission products (yield after evaporation).                                                          #
                #   Col 3: Fission fragments (data from GEF run).                                                               #
                #---------------------------------------------------------------------------------------------------------------#
                elif file_name.is_file() and file_name.name.startswith(data_to_read[7]):
                    path_file = os.path.join(pth_TALYS_CWD,file_name.name)
                    with open(path_file,'r') as file:   
                        data = file.readlines()
                    TALYS_data[f'{data_to_read[7]}_Number_of_nuclides'] = [data[2].split()[4]]
                    TALYS_data[f'{data_to_read[7]}'] = [np.loadtxt(path_file,dtype=np.float32,usecols=(0,1,2))]
                    del path_file,file,data
        except FileNotFoundError as e:
                sys.exit(e)
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------------#
        #   Functions for deleting simulation files and perturbed ff-files after a TALYS simulation is finished. The data from the files are stored in the McPUFF output. #
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------------#
        Reaction.delete_TALYS_result_files(pth_TALYS_CWD)                       # Comment out this line in order to retain TALYS output files and folders.
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------------#
        Reaction.delete_TALYS_ff_files(pth_TALYS_FY_folder,unique_thread_ID)    # Comment out this line in order to retain perturbed ".ff" files in GEF ".ff" library.
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------------#
        return TALYS_data
    #------------------------------------------------------- END OF METHOD ---------------------------------------------------------#

#----------------------------------------------------------- END OF CLASS ----------------------------------------------------------#
    
class Custom_Thread(Thread):
    """Custom thread class with capabilities to store output for
    unperturbed simulations.

    Parameters
    ----------
    target : callable
        ``Reaction.create_unperturbed_FY()`` function with arguments:
        - ``Z_target`` : Atomic number of target chemical element (`int`).
        - ``A_compound`` : Mass number of compound target nuclei (`int`).
        - ``E_reaction`` : Kinetic energy (MeV) of incident neutron in 
                           fission event (`float`).
        - ``MC_runs`` : The number of Monte Carlo simulations to be 
                        performed in GEF (`int`).
        - ``With_TALYS_flag`` : If `True`, includes TALYS evaporation 
                                simulation (`boolean`).
        - ``path_dict`` : Dictionary with local paths needed for 
                          simulations `dict` [`str`, `str`, `str`,
                          `str`, `str`, `str`, `str`, `str`].

    See Also
    --------
    ``Reaction.create_unperturbed_FY()``
    ``Reaction.read_and_clear_GEF_results()``
    ``Reaction.read_and_clear_TALYS_results()``
    ``Reaction.FY_results()``
    ``Reaction.GEF_FY_for_TALYS`()``

    Notes
    -----
    Simulation for unperturbed simulations performed separately because
    it only has to do one simulation and its faster/easier than to 
    include it in the `concurrent.futures`-multithread simulations. Since
    a separate thread cannot communicate with the main thread, it is not
    possible to store output data from the thread. Hence, a custom class
    that overloads the built-in Python `Thread` class is used.
    
    References
    ----------
    .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
    GEF and its influence on the nuclear evaporation in TALYS." 
    Technical report UPTEC F 23065, Uppsala University, 2023, accessed
    1 Janauary 2024, 
    url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>.

    Examples
    --------
    To perform an unperturbed thermal fission simulation for U235 + n 
    with 1e6 Monte Carlo fission events in GEF and no TALYS simulation.

    >>> unperturbed_thread = 
                    Custom_Thread(target=Reaction.create_unperturbed_FY, 
                        args=[92,236,2.53e-8,1e6,False,path_dict])
    [unperturbed_thread (Custom_Thread object)]     
    """
    #--------------------------------------------------- CLASS ATTRIBUTES ------------------------------------------------------#
    unperturbed_FY = None
    """Array with unperturbed GEF fission fragment yield data in ".ff" 
    file format 
    
    For GEF "lmd" option:  (`numpy.ndarray` (300,11)).
    For GEF "lmd+" option: (`numpy.ndarray` (300,19)).
    See the `See Also` section."""

    unperturbed_ignored_events = None
    """Number of GEF fission events removed due to statistical reasons 
    because they only occur once (`int`) [1]_."""

    unperturbed_GEF_results = None
    """Dictionary of unperturbed GEF simulation results (`dict`) 
    (See the `See Also` section.)"""

    unperturbed_TALYS_results = None
    """Dictionary of unperturbed TALYS simulation results (`dict`) 
    (See the `See Also` section.)"""
    #------------------------------------------------------ CONSTRUCTOR --------------------------------------------------------#
    def __init__(self,target=None,args=[]):
        Thread.__init__(self,target=None,args=[])
        self.target     = target
        self.Z_target        = args[0]
        self.A_compound      = args[1]
        self.E_reaction      = args[2]
        self.runs_MC         = args[3]
        self.With_TALYS_flag = args[4]
        self.path_dict       = args[5]
        self.unperturbed_FY = None
        self.unperturbed_ignored_events = None
        self.unperturbed_GEF_results = {}
        self.unperturbed_TALYS_results = {}
    #------------------------------------------------------ CLASS METHODS ------------------------------------------------------#
    def run(self):
        """Custom run method to store output data from thread target.

        Overloads `threading.Thread.run()` method.
        
        Parameters
        ----------
        Function takes no arguments.

        Returns
        -------
        For GEF "lmd" option:
        result_dict : `dict` [`numpy.ndarray` (300,11), `dict`, `dict`,`int`]
        For GEF "lmd+" option:
        result_dict : `dict` [`numpy.ndarray` (300,19), `dict`, `dict`,`int`]

        See Also
        --------
        `threading.Thread.run()`. 
        url: <https://docs.python.org/3/library/threading.html?highlight=thread#module-threading>

        Notes
        -----
        In case ``With_TALYS_flag`` = "False", 
        Reaction.create_unperturbed_FY() returns ``TALYS_data_dict`` 
        as an empty dictionary.
        
        Examples
        --------
        Custom thread stores output from thread in a dictionary.

        For GEF "lmd" option:
        >>> result_dict = self.target(92,236,2,1e5,True,path_dict)
        [result_dict dict [numpy.ndarray (300,11), dict, dict, int] ]

        For GEF "lmd+" option:
        >>> result_dict = self.target(92,236,2,1e5,True,path_dict)
        [result_dict dict [numpy.ndarray (300,19), dict, dict, int] ]
        """
        #-----------------------------------------------------------------------------------------------------------------------#
        result_dict = self.target(self.Z_target,self.A_compound,self.E_reaction,self.runs_MC,self.With_TALYS_flag,self.path_dict)
        self.unperturbed_FY = result_dict['FY_TALYS_format']
        self.unperturbed_GEF_results = result_dict['GEF_data_dict']
        self.unperturbed_TALYS_results = result_dict['TALYS_data_dict'] 
        self.unperturbed_ignored_events = result_dict['unpert_ignored_events']
#------------------------------------------------------- END OF CLASS ----------------------------------------------------------#
        
class Modified_Parameter(Reaction):
    """Create object that holds all simulation information about a 
    specific GEF parameter for ``Single_Parameter`` mode.

    The object holds the name and default value of a GEF parameter when 
    performing simulations using the ``Single_Parameters`` mode. For each 
    random number, a ``Random_Parameter_value`` object is created. This 
    holds all information about a parameter for a simulation with a
    specific random number. These objects are stored in a list in the
    ``Modified_Parameter`` object. Simulations for several parameters
    may be performed asynchronously and the results are kept separate
    by storing them in their own ``Modified_Parameter`` object.

    Parameters
    ----------
    param_name : `str`
        Name of GEF parameter. 
    unpert_param_val : `float`
        Unperturbed default GEF parameter value.
    num_of_rands : `int`
        Number of random numbers to produce using the specified 
        distribution.

        This number equals the number of simulations performed in
        McPUFF.
    distribution_flag : `str`
            Specifies if random numbers are to be drawn from "uniform", 
            "normal" or "max-min" distributions.

    See Also
    --------
    Rough estimates of GEF parameter stability limits: 
    url:<https://github.com/UPTEC-F-23065/McPUFF/blob/12b9f864240fe133ef727955c9e7715a7e0e20a9/Rough_limits_GEF_parameter_values.txt>

    ``package_McPUFF.Gaussian_GEF_Param.gaussian_st_dev_for_parameter()``.
    ``number_of_randoms`` in `Reaction.__init__`.
    
    Notes
    -----
    A list of the GEF parameter names and default values can be found in 
    appendix A in reference [1]_.

    When using the "max-min" distribution, the Max/Min values are given as 
    "random numbers" in a list, e.g self.list_of_rand_num = [1,1.5]. The 
    numbers are used to scale the unperturbed parameter values. I.e. 1.5 
    equals 1.5 * unperturbed parameter value. It is possible to run 
    simulations for more values, e.g. self.list_of_rand_num = 
    [0.95,0.97,0.99,1.01,1.03,1.05]. When doing so, the number of values 
    in the list must be matched by the value of the "number_of_randoms"-
    parameter given in the main file (See the "See also" section). I.e. 4 
    values in the "max-min" list means that the "number_of_randoms" 
    parameter must equal 4. The method uses the framework of random 
    perturbations for TMC-simulations, hence the name "random_values" for 
    the parameter values.

    References
    ----------
    .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
    GEF and its influence on the nuclear evaporation in TALYS." 
    Technical report UPTEC F 23065, Uppsala University, 2023, accessed
    1 Janauary 2024, 
    url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>.

    Examples
    --------    
    Example produces an object that holds information about a parameter
    as well as a list of objects with results for perturbed simulations.

    >>> mod_param_object = Modified_Parameter("_P_DZ_Mean_S1",- 0.1710014,
                            500,"normal")
    [mod_param_object (Modified_Parameter object)]
    """
    #--------------------------------------------------- CLASS ATTRIBUTES ------------------------------------------------------#
    list_of_Random_Parameter_Value_objects = None
    """List that holds objects with results for a perturbed simulation of
     a parameter using a specific random number. """
    
    scaling_number_uniform = None
    """Scaling number for uniform distribution.
    
    Lets user adjust span for the uniform distribution. (See the 
    "See Also" section for parameter stability limits). Creates the upper
    and lower limits for the span of the uniform distribution. The 
    scaling is given as a percentage of the unperturbed parameter value. 
    E.g. scaling_number_uniform = 0.5 sets the span of the uniform 
    distribution as: [0.5*unperturbed_value, 1.5*unperturbed_value].
    """

    upper_lim = None
    """Upper limit of uniform distribution."""

    lower_lim = None
    """Lower limit of uniform distribution."""

    list_of_rand_num = None
    """List of random numbers is created here because it is unique to
    each parameter."""

    GEF_st_dev = None
    """Standard deviation of the Gaussian (normal) distribution is given 
    as a percentage of the unperturbed parameter value in the file: 
    ``Gaussian_GEF_Param.py`` (See the "See Also" section)."""
    #------------------------------------------------------ CONSTRUCTOR --------------------------------------------------------#
    def __init__(self, param_name,unpert_param_val,num_of_rands,distribution_flag): 
        self.param_name = param_name                                
        self.unpert_param_val = np.float32(unpert_param_val)                                                                        # Unperturbed value used for centering distributions on unperturbed parameter value.
        self.list_of_Random_Parameter_Value_objects = []            
        #---------------------------------- LIST OF RANDOM NUMBERS FOR "UNIFORM" DISTRIBUTION ----------------------------------#
        if distribution_flag == 'uniform':
            scaling_number_uniform = np.float32(0.5)                                                    # Scaling value for limits of distribution. User can adjust the span of the uniform distribution here.
            upper_lim = abs(self.unpert_param_val) + scaling_number_uniform*abs(self.unpert_param_val)                              # Upper numerical limit of distribution.
            lower_lim = abs(self.unpert_param_val) - scaling_number_uniform*abs(self.unpert_param_val)                              # Lower numerical limit of distribution.
            self.list_of_rand_num = (np.random.default_rng().uniform(low=lower_lim,high=upper_lim,size=num_of_rands)).astype(dtype=np.float32) # Creation of list of random numbers.
            del distribution_flag,scaling_number_uniform,upper_lim,lower_lim
        #----------------------------------- LIST OF RANDOM NUMBERS FOR "NORMAL" DISTRIBUTION ----------------------------------#
        elif distribution_flag == 'normal':
            GEF_st_dev = abs(Gaussian_GEF_Param.gaussian_st_dev_for_parameter(param_name,unpert_param_val))                         # Retrieval of std from external input file.
            self.list_of_rand_num = (np.random.default_rng().normal(self.unpert_param_val,scale=float(GEF_st_dev),size=num_of_rands)).astype(np.float32)    # Creation of list of random numbers.
            del distribution_flag,GEF_st_dev
        #----------------------------------- LIST OF RANDOM NUMBERS FOR "MAX-MIN" DISTRIBUTION ---------------------------------#
        elif distribution_flag == 'max-min':
            self.list_of_rand_num = [1,1.5]                                                             # List of values used to scale unperturbed parameter values in Max/Min-mode. (See description of variable).   
            del distribution_flag     
        else:
            print(f'The "{distribution_flag}"-distribution does not exist')
            sys.exit('Exiting program')
    #------------------------------------------------------ CLASS METHODS ------------------------------------------------------#
    def create_perturbed_parameter_value(param_name,rand_num,unpert_param_val,distribution_flag):
        """Perturb GEF parameter by adding a random number to the GEF
        default value.
        
        Parameters
        ----------
        param_name : `str`
            Name of GEF parameter.
        rand_num : `float`
            Random number added to default GEF parameter value in order to
            perturb the parameter.
        unpert_param_val : `float`
            Default GEF parameter value. 
        distribution_flag : `str`
            Name of distribution from which the random number was drawn.
            
        Returns
        -------
        pert_param_val : `float`
            Perturbed GEF parameter value.

        See Also
        --------
        Class ``Modified_Parameter``
        ``Gaussian_GEF_Param.gaussian_st_dev_for_parameter()``

        Notes
        -----
        For parameters whose default values are not zero, the "uniform", 
        "normal" and "max-min" distributions are centered on the GEF parameter
        default value. 

        The parameters in the list ``special_case_parameters`` are the GEF
        parameters whose default values are zero. Since it is meaningless
        to create a perturbed parameter value as a percentage of the 
        default value zero, these parameters are appointed a random value
        chosen by the user. The procedure depends on which distribution 
        that is used, as follows (See code for each distribution below):

        **Perturbed value for "uniform" distribution in range:**
        [-scaling_special_case_parameters,+scaling_special_case_parameters]
        (Centered on zero).

        **Perturbed value for "normal" distribution in range:**
        [-st_dev_special_param, +st_dev_special_param] (centered on zero).
          
        **Perturbed value for "max-min" distribution in range:**
        [(min-1), (max-1)]
        ``rand_num`` in the code below is not random but a value chosen by 
        the user. See class ``Modified_Parameter`` for more information.
             
        These ranges and scaling parameters can be changed by the user.
        
        Examples
        --------
        Example for "uniform" distribution and "special_case_parameter":
        >>> pert_val = Reaction.create_perturbed_parameter_value("ZPOL1",
                                                            0,0,"uniform")
        [-0.15] (for random number = 0.35)
        
        Example for "normal" distribution:
        >>> pert_val = Reaction.create_perturbed_parameter_value("_betaH0",
                                                45.6133,45.67577,"normal")
        [45.6133]

        Example for "max-min" distribution:
        >>> pert_val = Reaction.create_perturbed_parameter_value(
                            "_dE_Defo_S1",1.5,-4.254451,"max-min")
        [-6.3816765]

        Example for "max-min" distribution and "special_case_parameter":
        >>> pert_val = Reaction.create_perturbed_parameter_value("kappa",
                                                1.5,0,"max-min")
        [0.5]
        """
        #-----------------------------------------------------------------------------------------------------------------------#
        special_case_parameters = ['_Delta_S0','P_att_pol2','P_att_pol3','kappa','Epot_shift','SIGDEFO_slope','ZPOL1','P_n_x','T_orbital']
        """List of all GEF parameters with default value = 0."""
        param_name = str(param_name)
        """Name of GEF parameter."""
        rand_num = np.float32(rand_num)
        """Random number used to perturb default GEF parameter value."""
        unpert_param_val = np.float32(unpert_param_val)
        """Default GEF parameter value."""
        distribution_flag = str(distribution_flag)
        """Distribution from which random number is drawn."""
        #-------------------------- CREATION OF PERTURBED PARAMETER VALUES FOR "UNIFORM" DISTRIBUTION---------------------------#
        if distribution_flag == 'uniform':
            if param_name in special_case_parameters:                                                                       # Cannot produce perturbation as percentage if default value = 0. Set value using special case scaling parameter.   
                scaling_special_case_parameters = np.float32(0.5)                                                           # The user can change this value to set the range of the special case random number distributions.     
                special_case_rand_num = np.random.rand(1).astype('float32')                                                 # One random number in range [0,1[
                pert_param_val = -scaling_special_case_parameters + special_case_rand_num*2*scaling_special_case_parameters
                pert_param_val = round(pert_param_val[0],9)                                                                 # [0] Important, otherwise pert_param_val is list.
                del special_case_parameters,param_name,rand_num,unpert_param_val,distribution_flag,scaling_special_case_parameters,special_case_rand_num
                return pert_param_val
            else:                         
                pert_param_val = round((unpert_param_val/abs(unpert_param_val))*rand_num,9)                                 # Multiply random number by normalized unperturbed parameter value to get correct sign.             
                del special_case_parameters,param_name,rand_num,unpert_param_val,distribution_flag
                return pert_param_val
        #-------------------------- CREATION OF PERTURBED PARAMETER VALUES FOR "NORMAL" DISTRIBUTION----------------------------#
        elif distribution_flag == 'normal':                  
            if param_name in special_case_parameters:                                                                       # Cannot produce perturbation as percentage if default value = 0. Set value using special case scaling parameter.
                st_dev_special_param = np.float32(0.03)                                                                     # The user can change this value to set the range of the special case random number distributions.
                special_case_rand_num = (np.random.default_rng().normal(loc=0.0,scale=st_dev_special_param,size=1)).astype(dtype=np.float32)                                                                                                                  
                pert_param_val = round(special_case_rand_num[0],9)                                                          # [0] Important, otherwise pert_param_val is list.
                del special_case_parameters,param_name,rand_num,unpert_param_val,distribution_flag,st_dev_special_param,special_case_rand_num
                return pert_param_val
            else:
                pert_param_val = round(rand_num,9)                                                                          # Random number created using ``Gaussian_GEF_Param.gaussian_st_dev_for_parameter()``. See the "See Also" section.
                del special_case_parameters,param_name,rand_num,unpert_param_val,distribution_flag
                return pert_param_val
        #-------------------------- CREATION OF PERTURBED PARAMETER VALUES FOR "MAX-MIN" DISTRIBUTION---------------------------#
        elif distribution_flag == 'max-min':
            if param_name in special_case_parameters:                                                                       # Cannot produce perturbation as percentage if default value = 0. Set value using special case scaling parameter.
                pert_param_val = round(rand_num-1,9)                                                                        # The user can change this value to set the range of the special case random number distributions.        
                del special_case_parameters,param_name,rand_num,unpert_param_val,distribution_flag
                return pert_param_val
            else:               
                pert_param_val = round(unpert_param_val*rand_num,9)             
                del special_case_parameters,param_name,rand_num,unpert_param_val,distribution_flag
                return pert_param_val
        else:
            print(f'The distribution named {distribution_flag} is not available')
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#
    
    def create_perturbed_FY(param_name,unpert_param_val,rand_num,enumeration_rand_num,Z_target,A_compound,E_reaction,MC_runs,With_TALYS_flag,distribution_flag,path_dict):
        """Perform simulation with perturbed GEF parameter values using
        ``Single_Parameter`` mode.

        This function performs the actual perturbed simulations. It calls
        GEF and TALYS during execution and collects all results.
        
        Parameters
        ----------
        param_name : `str`
            Name of GEF parameter.
        unpert_param_val : `float`
            Default value of GEF parameter.
        rand_num : `float`
            Random number used to perturb the default GEF parameter value.
        enumeration_rand_num : `int`
            Simulation count used to create unique ID for simulation.

            The total number of simulations is determined at the outset by
            the user as the number of random numbers. The simulation count 
            is unique and is appended to the simulation files and folders 
            to keep them separate when multi-threading.
        Z_target : `int`
            Atomic number of target chemical element. 
        A_compound : `int`
            Mass number of compound target nuclei. I.e A + 1.
        E_reaction : `float`
            Kinetic energy (MeV) of incident neutron in fission event.
        MC_runs : `int`
            The number of Monte Carlo simulations to be performed in GEF. 
        With_TALYS_flag : `boolean`
            If `True`, includes TALYS evaporation simulation. 
        distribution_flag : `str`
            Specifies if random numbers are to be drawn from "uniform", 
            "normal" or "max-min" distributions.
        path_dict : `dict` [`str`,`str`,`str`,`str`,`str`,`str`,`str`]
            Dictionary with local paths needed for simulations.
        
        Returns
        -------
        rand_param_val_obj : ``Random_Parameter_value``
            Parameter object with attributes:

            ``param_name``
                Name of GEF parameter (`str`).
            ``unpert_param_val``
                Default value of GEF parameter (`float`). 
            ``random_value``
                Random number used to perturb the default GEF parameter
                value (`float`).
            ``enumeration_rand_num``
                Simulation count used to create unique ID for 
                simulation (`int`).
            ``ignored_events``
                Number of GEF fission events removed due to statistical
                reasons because they only occur once (`int`).
            ``pert_param_val``
                Perturbed GEF parameter value used in simulation (`float`).
            ``perturbed_FY``
                Array of GEF fission fragment yields in GEF ".ff" library
                format. `dtype` = numpy.float(32). 
                
                For GEF "lmd" option:  (`numpy.ndarray`, (300,11)).
                For GEF "lmd+" option: (`numpy.ndarray`, (300,19)).
            ``perturbed_GEF_results``
                Dictionary with GEF simulation results (`dict`).
            ``perturbed_TALYS_results``
                Dictionary with TALYS simulation results (`dict`).

        See Also
        --------
        Class ``Random_Parameter_value``
        ``package_McPUFF.TALYS_Input``
        ``Reaction.FY_results()``
        ``Reaction.GEF_FY_for_TALYS()``

        Notes
        -----
        The ``unique_pert_thread_ID`` is very important. It is needed to
        keep simulation results separate during multi-thread simulations
        and is also used to differentiate between perturbed and 
        unperturbed fission fragment yield files (".ff") in the GEF 
        library. It is passed to TALYS as the value of a specially created
        keyword called "geffissionfileid" [1]_. See the "See Also" section for
        more information.

        References
        ---------
        .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
        GEF and its influence on the nuclear evaporation in TALYS." 
        Technical report UPTEC F 23065, Uppsala University, 2023, accessed
        1 Janauary 2024, 
        url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>.

        Examples
        --------
        Performs perturbed simulations with GEF and TALYS for 
        the ``Single_Parameters`` mode.
        >>> Sim_obj = Reaction.create_perturbed_FY("T_low_S11",0.36,
                0.360927,342,92,236,2.53e-8,1e6,True,"normal",path_dict)
        [Sim_obj (Random_Parameter_value object)]
        """
        #------------------------------------- CREATE UNIQUE THREAD NAME FOR MULTI-THREADING -----------------------------------#
        unique_pert_thread_ID = f'Param_{param_name}_{enumeration_rand_num}'                                        # Used to separate simulations when multi-threading.
        #---------------------------------- CREATE ``RANDOM_PARAMETER_VALUE`` OBJECT TO STORE RESULTS --------------------------#
        rand_param_val_obj = Random_Parameter_value(param_name,unpert_param_val,rand_num,enumeration_rand_num)      # Holds all perturbed simulation results for simulation with specific random number.
        #----------------------------------------------- CREATE PERTURBED PARAMETER VALUES -------------------------------------#
        rand_param_val_obj.pert_param_val = Modified_Parameter.create_perturbed_parameter_value(param_name,rand_num,unpert_param_val,distribution_flag)     
        #----------------------------------------------- CREATE GEF PATHS FOR SIMULATION ---------------------------------------#     
        dict_of_GEF_paths = Reaction.create_GEF_workingdir_and_inputfile(unique_pert_thread_ID,Z_target,A_compound,E_reaction,MC_runs,path_dict['GEF_working_dir_path'])    # Add more paths if other files are to be read.
        GEF_cwd_path = dict_of_GEF_paths['GEF_cwd_path']                                                            # Local path to individual perturbed GEF simulation output folder.
        GEF_LMD_path = dict_of_GEF_paths['GEF_LMD_path']                                                            # Local path to individual perturbed GEF simulation ".lmd" file.
        GEF_DAT_path = dict_of_GEF_paths['GEF_DAT_path']                                                            # Local path to individual perturbed GEF simulation ".dat" file.
        GEF_DMP_EN_path = dict_of_GEF_paths['GEF_DMP_EN_path']                                                      # Local path to individual perturbed GEF simulation "EN.dat" file.
        Reaction.eraseFolders(GEF_cwd_path)                                                                         # Make sure folders are erased before start. Otherwise GEF exits because it thinks the simulation has already been performed.
        #-------------------------------- WRITE PERTURBED PARAMETER VALUES IN MY_PARAMETERS.DAT --------------------------------#
        with open(os.path.join(GEF_cwd_path,'MyParameters.dat'), 'w') as inputMyParameters_dat:
            row = param_name+' = '+str(rand_param_val_obj.pert_param_val)                                           # Enter perturbed GEF parameter value to be used in simulation.
            row = ''.join(map(str,row))
            inputMyParameters_dat.write(row) 
        #----------------------------------------------- RUN PERTURBED GEF SIMULATION ------------------------------------------#
        with subprocess.Popen(["GEF"],stdin=PIPE,stdout=DEVNULL,cwd=GEF_cwd_path) as GEF_process:                   # Run perturbed GEF simulation. Make sure modified GEF is on path.
            GEF_process.communicate(input=b'\n')
        rand_param_val_obj.perturbed_FY, rand_param_val_obj.ignored_events = Reaction.FY_results(GEF_LMD_path,A_compound,Z_target)
        rand_param_val_obj.perturbed_GEF_results = Reaction.read_and_clear_GEF_results(GEF_cwd_path,GEF_DAT_path,GEF_DMP_EN_path)   # Reads and stores data from ".dat" file. Deletes files and folders after data is stored in object.
        del GEF_cwd_path,GEF_LMD_path,GEF_DAT_path,GEF_process
        if With_TALYS_flag == True:                                                                                 # If flag set to TRUE, run TALYS.
            #----------------------------------- CREATE PERTURBED .FF FILE AND TALYS PATHS FOR SIMULATION ----------------------#                                                                          
            Reaction.print_TALYS_ff_files(path_dict['TALYS_ff_file_path'],rand_param_val_obj.perturbed_FY,Z_target,A_compound,E_reaction,unique_pert_thread_ID)
            TALYS_pert_input_command = 'talys'  
            TALYS_pert_input_file    = f'{unique_pert_thread_ID}_input.in'                                 
            TALYS_pert_output_file   = f'{unique_pert_thread_ID}_output.out'
            TALYS_pert_cwd           = os.path.join(path_dict['TALYS_input_path']+f'{unique_pert_thread_ID}',"")    # Unique TALYS output folder for each perturbed simulation.
            TALYS_Input.create_TALYS_input_file(path_dict['TALYS_input_path'],unique_pert_thread_ID,TALYS_pert_input_file,Z_target,A_compound,E_reaction)
            #---------------------------------------------- RUN PERTURBED TALYS SIMULATION -------------------------------------#
            try:
                subprocess.run(["bash","-c",str(f"{TALYS_pert_input_command} < {TALYS_pert_input_file} > {TALYS_pert_output_file}")],cwd=TALYS_pert_cwd)    # Run perturbed TALYS simulation. Make sure modified TALYS is on path.
            except FileNotFoundError as e:
                print(f'TALYS could not run because it cannot find the input file.\n{e}')
                sys.exit(e)
            except subprocess.CalledProcessError as e:
                print(f'An error occured while running TALYS\n{e}')
                sys.exit()
            #------------------------- STORE PERTURBED TALYS DATA AND DELETE OUTPUT FILES AND FOLDERS --------------------------#
            rand_param_val_obj.perturbed_TALYS_results = Reaction.read_and_clear_TALYS_results(TALYS_pert_cwd,path_dict['TALYS_ff_file_path'],unique_pert_thread_ID)    # Store selected TALYS output and then delete files and folders.
            del TALYS_pert_input_command,TALYS_pert_input_file,TALYS_pert_output_file,TALYS_pert_cwd
        #------------------------------------------- IMPORTANT INDENTATION -----------------------------------------------------#
        del unique_pert_thread_ID
        return rand_param_val_obj       
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#

    def create_perturbed_TMC_FY(dict_unpert_param_name_val,enumeration_rand_num,Z_target,A_compound,E_reaction,MC_runs,path_dict,With_TALYS_flag,distribution_flag):
        """Perform simulation with perturbed GEF parameter values using
        ``TMC`` mode.

        This function performs the actual perturbed simulations. It calls
        GEF and TALYS during execution and collects all results.
        
        Parameters
        ----------
        dict_unpert_param_name_val : `dict` [`str`, `float`]
            The keys are the GEF parameter names and the values are the
            default GEF parameter values.
        enumeration_rand_num : `int`
            Simulation count used to create unique ID for simulation.

            The total number of simulations is determined at the outset by
            the user as the number of random numbers. The simulation count
            is unique and is appended to the simulation files and folders 
            to keep them separate when multi-threading.
        Z_target : `int`
            Atomic number of target chemical element. 
        A_compound : `int`
            Mass number of compound target nuclei. I.e A + 1.
        E_reaction : `float`
            Kinetic energy (MeV) of incident neutron in fission event.
        MC_runs : `int`
            The number of Monte Carlo simulations to be performed in GEF. 
        path_dict : `dict` [`str`,`str`,`str`,`str`,`str`,`str`,`str`]
            Dictionary with local paths needed for simulations.
        With_TALYS_flag : `boolean`
            If `True`, includes TALYS evaporation simulation. 
        distribution_flag : `str`
            Specifies if random numbers are to be drawn from "uniform", 
            "normal" or "max-min" distributions.
        
        Returns
        -------
        tmc_obj : ``TMC_Object``
            Simulation object with attributes:

            ``dictionary_unpert_param_name_val`` 
                Dictionary with GEF parameter names and default parameter 
                values (`dict` [`str`, `float`])
            ``enumeration_rand_val``
                 Simulation count used to create unique ID for simulation
                 (`int`).
            ``list_of_TMC_Mod_Param_objects`` 
                List of objects holding information about individual 
                parameters (`list` [``TMC_Mod_Param_object``]).
            ``ignored_events`` 
                Number of GEF fission events removed for statistical 
                reasons because they only occur once (`int`).
            ``perturbed_FY``
                Array with perturbed GEF fission fragment yield data from
                simulations. `dtype` = numpy.float(32). 
                
                For GEF "lmd" option:  (`numpy.ndarray`, (300,11)).
                For GEF "lmd+" option: (`numpy.ndarray`, (300,19)).
            ``perturbed_GEF_results``
                Dictionary with GEF simulation results (`dict`).
            ``perturbed_TALYS_results``
                Dictionary with TALYS simulation results (`dict`).
        See Also
        --------
        Class ``TMC_Mod_Param_object``
        Class ``TMC_Object``
        ``package_McPUFF.TALYS_Input``
        ``Reaction.FY_results()``
        ``Reaction.GEF_FY_for_TALYS()``

        Notes
        -----
        The ``unique_pert_thread_ID`` is very important. It is needed to
        keep simulation results separate during multi-thread simulations
        and is also used to differentiate between perturbed and 
        unperturbed fission fragment yield files (".ff") in the GEF 
        library. It is passed to TALYS as the value of a specially created
        keyword called "geffissionfileid" [1]_. See the "See Also" section for
        more information.

        References
        ----------
        .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
        GEF and its influence on the nuclear evaporation in TALYS." 
        Technical report UPTEC F 23065, Uppsala University, 2023, accessed
        1 Janauary 2024, 
        url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>.

        Examples
        --------
        Performs perturbed simulations with GEF and TALYS for the ``TMC`` mode.
        >>> Sim_obj = Reaction.create_perturbed_TMC_FY(dictParamNameVal,
                         572,92,236,2.53e-8,1e6,path_dict,True,"normal")
        [Sim_obj (TMC_Object object)]  
        """
        #------------------------------------- CREATE UNIQUE THREAD NAME FOR MULTI-THREADING -----------------------------------#
        unique_pert_thread_ID = f'TMC_{enumeration_rand_num}'                                                                                                   # Used to separate simulations when multi-threading.
        #------------------------------------------ CREATE ``TMC`` OBJECT TO STORE RESULTS -------------------------------------#
        tmc_obj = TMC_Object(dict_unpert_param_name_val, enumeration_rand_num)     
        #------------------------ CREATE INDIVIDUAL PARAMETER OBJECTS AND PERTURBED PARAMETER VALUES ---------------------------#  
        for param_name, unpert_param_val in dict_unpert_param_name_val.items():                    
            tmc_mod_param_obj = TMC_Mod_Param_object(param_name,unpert_param_val,distribution_flag)                                                             # Holds all perturbed simulation results for simulation with specific random number.
            tmc_mod_param_obj.pert_param_val = Modified_Parameter.create_perturbed_parameter_value(param_name,tmc_mod_param_obj.random_value,unpert_param_val,distribution_flag)
            tmc_obj.list_of_TMC_Mod_Param_objects.append(tmc_mod_param_obj)                                                                                     # Must append here to loop over all chosen parameters in GEF run.
        #----------------------------------------------- CREATE GEF PATHS FOR SIMULATION ---------------------------------------#
        dict_of_GEF_paths = Reaction.create_GEF_workingdir_and_inputfile(unique_pert_thread_ID,Z_target,A_compound,E_reaction,MC_runs,path_dict['GEF_working_dir_path'])
        GEF_cwd_path = dict_of_GEF_paths['GEF_cwd_path']                                                                                                        # Local path to individual perturbed GEF simulation output folder.
        GEF_LMD_path = dict_of_GEF_paths['GEF_LMD_path']                                                                                                        # Local path to individual perturbed GEF simulation ".lmd" file.
        GEF_DAT_path = dict_of_GEF_paths['GEF_DAT_path']                                                                                                        # Local path to individual perturbed GEF simulation ".dat" file.
        GEF_DMP_EN_path = dict_of_GEF_paths['GEF_DMP_EN_path']                                                                                                  # Local path to individual perturbed GEF simulation "EN.dat" file.
        #-------------------------------- WRITE PERTURBED PARAMETER VALUES IN MY_PARAMETERS.DAT --------------------------------#
        Reaction.eraseFolders(GEF_cwd_path)                                                                                                                     # Make sure folders are erased before start. Otherwise GEF exits because it thinks the simulation has already been performed.
        Reaction.clear_MyParameters_dat(GEF_cwd_path)                                       
        for tmc_mod_param_obj in tmc_obj.list_of_TMC_Mod_Param_objects:
            with open(os.path.join(GEF_cwd_path,'MyParameters.dat'), 'a') as inputMyParameters_dat:
                row = tmc_mod_param_obj.param_name+' = '+str(tmc_mod_param_obj.pert_param_val)                                                                  # Enter perturbed GEF parameter value to be used in simulation.
                row = ''.join(map(str,row))
                inputMyParameters_dat.write(row+'\n')
        #----------------------------------------------- RUN PERTURBED GEF SIMULATION ------------------------------------------#
        with subprocess.Popen(["GEF"],stdin=PIPE,stdout=DEVNULL,cwd=GEF_cwd_path) as GEF_process:                                                               # Run perturbed GEF simulation. Make sure modified GEF is on path.
            GEF_process.communicate(input=b'\n')
        tmc_obj.perturbed_FY, tmc_obj.ignored_events = Reaction.FY_results(GEF_LMD_path,A_compound,Z_target)
        tmc_obj.perturbed_GEF_results = Reaction.read_and_clear_GEF_results(GEF_cwd_path,GEF_DAT_path,GEF_DMP_EN_path)                                          # Reads and stores data from ".dat" file. Deletes files and folders after data is stored in object.
        if With_TALYS_flag == True:                                                                                                                             # If flag set to TRUE, run TALYS      
            #------------------------- CREATE PERTURBED .FF FILE AND TALYS PATHS FOR SIMULATION --------------------------------#
            Reaction.print_TALYS_ff_files(path_dict['TALYS_ff_file_path'],tmc_obj.perturbed_FY,Z_target,A_compound,E_reaction,unique_pert_thread_ID)
            TALYS_pert_input_command = 'talys'  
            TALYS_pert_input_file    = f'{unique_pert_thread_ID}_input.in'                       
            TALYS_pert_output_file   = f'{unique_pert_thread_ID}_output.out'
            TALYS_pert_cwd           = os.path.join(path_dict['TALYS_input_path']+f'{unique_pert_thread_ID}',"")                                                # Unique TALYS output folder for each perturbed simulation.
            TALYS_Input.create_TALYS_input_file(path_dict['TALYS_input_path'],unique_pert_thread_ID,TALYS_pert_input_file,Z_target,A_compound,E_reaction)
            #---------------------------------------------- RUN PERTURBED TALYS SIMULATION -------------------------------------#
            try:
                subprocess.run(["bash","-c",str(f"{TALYS_pert_input_command} < {TALYS_pert_input_file} > {TALYS_pert_output_file}")],cwd=TALYS_pert_cwd)        # Run perturbed TALYS simulation. Make sure modified TALYS is on path.
            except FileNotFoundError as e:
                print(f'TALYS could not run because it cannot find the input file.\n{e}')
                sys.exit(e)
            except subprocess.CalledProcessError as e:
                print(f'An error occured while running TALYS\n{e}')
                sys.exit()
            #------------------------- STORE PERTURBED TALYS DATA AND DELETE OUTPUT FILES AND FOLDERS --------------------------#
            tmc_obj.perturbed_TALYS_results = Reaction.read_and_clear_TALYS_results(TALYS_pert_cwd,path_dict['TALYS_ff_file_path'],unique_pert_thread_ID)       # Store selected TALYS output and then delete files and folders.
            del TALYS_pert_input_command,TALYS_pert_input_file,TALYS_pert_output_file,TALYS_pert_cwd
        #------------------------------------------- IMPORTANT INDENTATION -----------------------------------------------------#
        del unique_pert_thread_ID,GEF_cwd_path,GEF_LMD_path,GEF_DAT_path,GEF_process
        return tmc_obj
    #----------------------------------------------------- END OF METHOD -------------------------------------------------------#

#--------------------------------------------------------- END OF CLASS --------------------------------------------------------#
    
class Random_Parameter_value(Modified_Parameter):
    """Create object to hold perturbed simulation results for a specific
     random number when using the ``Single_Parameters`` mode. 
     
    Parameters
    ----------
    parameter_name : `str`
        Name of GEF parameter.
    unperturbed_param_val : `float`
        Default GEF parameter value.
    rand_val : `float`
        Random number used to perturb default GEF parameter value.
    enum_rand_num : `int`
        Simulation count used to create unique ID for simulation.

    See Also
    --------
    Class ``Modified_Parameter``
    ``Modified_Parameter.create_perturbed_parameter_value()``
    ``Gaussian_GEF_Param.gaussian_st_dev_for_parameter()``

    Notes
    -----
    When performing a "Total Monte Carlo" simulation [1]_, possibly 
    thousands of GEF and TALYS simulations are performed. In order to 
    retain data from each individual simulation, the data is stored in
    separate objects (``Random_Parameter_value`` objects for 
    ``Single_Parameters`` mode). This makes it possible to relate specific
    data for fission observables to specific parameter values and to 
    perform a correlation analysis between the fission observables and
    the GEF parameter values.

    References
    ----------
    .. [1] P. Karlsson, "Total Monte Carlo of the fission model in
        GEF and its influence on the nuclear evaporation in TALYS." 
        Technical report UPTEC F 23065, Uppsala University, 2023, accessed
        1 Janauary 2024, 
        url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>.

    Examples
    --------
    Instantiates a ``Random_Parameter_value`` object:
    >>> rand_param_val_obj = Random_Parameter_value("_P_Z_Curv_S1",
                                               0.3814453,-5.28966e-3,1924)
    [rand_param_val_obj (Random_Parameter_value object)]
    """
    #--------------------------------------------------- CLASS ATTRIBUTES ------------------------------------------------------#
    param_name = None
    """Name of GEF parameter (`str`)."""

    unpert_param_val = None
    """Default GEF parameter value (`float`)."""

    random_value = None
    """Random value used to perturb GEF parameters (`float`)."""

    enumeration_rand_num = None
    """Simulation count used to create unique ID (`str`)."""

    ignored_events = None
    """Number of GEF fission events removed for statistical reasons 
    (`int`)."""

    pert_param_val = None
    """Perturbed GEF parameter value (`float`)."""

    perturbed_FY = None
    """Array of GEF fission fragment yields in GEF ".ff" library
    format. `dtype` = numpy.float(32). 
    
    For GEF "lmd" option:  (`numpy.ndarray`, (300,11)).
    For GEF "lmd+" option: (`numpy.ndarray`, (300,19))."""

    perturbed_GEF_results = None
    """Dictionary with GEF simulation results (`dict`)."""

    perturbed_TALYS_results = None
    """Dictionary with TALYS simulation results (`dict`)."""
    #------------------------------------------------------ CONSTRUCTOR --------------------------------------------------------#
    def __init__(self, parameter_name, unperturbed_param_val, rand_val,enum_rand_num):
        self.param_name = parameter_name
        self.unpert_param_val = unperturbed_param_val
        self.random_value = rand_val
        self.enumeration_rand_num = enum_rand_num
        self.ignored_events = None
        self.pert_param_val = None
        self.perturbed_FY = None
        self.perturbed_GEF_results = {}
        self.perturbed_TALYS_results = {}
#------------------------------------------------------- END OF CLASS ----------------------------------------------------------#

class TMC_Object(Modified_Parameter):
    """Create object to hold perturbed simulation results for a specific
     random number when using the ``TMC`` mode. 
     
    Parameters
    ----------
    dict_unpert_param_name_val : `dict` [`str`, `float`]
        Dictionary of GEF parameter names and default parameter values.
    enum_rand_num : `int`
        Simulation count used to create unique ID for simulation.

    See Also
    --------
    ``TMC_Mod_Param_object``
    ``Random_Parameter_value``
    ``Modified_Parameter``

    Notes
    -----
    The ``TMC_Object`` object has the same role in a ``TMC`` simulation as 
    a ``Modified_Parameter`` object has in a ``Single_Parameters`` 
    simulation. In the ``Single_Parameters`` mode, the parameters are 
    perturbed one at a time and ``Random_Parameter_value`` objects store 
    all perturbed results for each perturbation. In the ``TMC`` mode many, 
    possibly all, parameters are perturbed simultaneously in the same 
    simulation and the information about each parameter is stored in a
    ``TMC_Mod_Param_object`` object. To avoid redundancy, all perturbed 
    simulation data is stored once in the ``TMC_object`` instead of 
    storing it in each ``TMC_Mod_Param_object`` object.
    
    Examples
    --------
    Instantiates a ``TMC_Object`` object:
    >>> tmc_obj = TMC_Object(dict_unpert_param_name_val, 325) 
    [tmc_obj (TMC_Object object)]
    """
    #--------------------------------------------------- CLASS ATTRIBUTES ------------------------------------------------------#
    dictionary_unpert_param_name_val = None
    """Dictionary of GEF parmeter names and default parameter values 
    (`dict` [`str`, `float`])."""

    enumeration_rand_val = None
    """Simulation count used to create unique ID for simulation (`int`)."""

    list_of_TMC_Mod_Param_objects = None
    """List for storing parameter objects 
    (`list` [`TMC_Mod_Param_object`])."""

    ignored_events = None
    """Number of GEF fission events removed for statistical reasons 
    (`int`)."""

    perturbed_FY = None
    """Array of GEF fission fragment yields in GEF ".ff" library
    format. `dtype` = numpy.float(32). 
    
    For GEF "lmd" option:  (`numpy.ndarray`, (300,11)).
    For GEF "lmd+" option: (`numpy.ndarray`, (300,19))."""

    perturbed_GEF_results = None
    """Dictionary with GEF simulation results (`dict`)."""

    perturbed_TALYS_results = None
    """Dictionary with TALYS simulation results (`dict`)."""
    #------------------------------------------------------ CONSTRUCTOR --------------------------------------------------------#
    def __init__(self, dict_unpert_param_name_val, enum_rand_val):
        self.dictionary_unpert_param_name_val = dict_unpert_param_name_val
        self.enumeration_rand_val = enum_rand_val
        self.list_of_TMC_Mod_Param_objects = []    
        self.ignored_events = None
        self.perturbed_FY = None
        self.perturbed_GEF_results = {}
        self.perturbed_TALYS_results = {}
#------------------------------------------------------- END OF CLASS ----------------------------------------------------------#
        
class TMC_Mod_Param_object(Modified_Parameter):
    """Create object to store information about individual parameter 
    perturbation during ``TMC`` simulation.
    
    Parameters
    ----------
    parameter_name : `str`
        Name of GEF parameter.
    unperturbed_param_val : `float`
        Default GEF parameter value.
    distribution_flag : `str`
            Specifies if random numbers are to be drawn from "uniform", 
            "normal" or "max-min" distributions.

    See Also
    --------
    ``Modified_Parameter``
    ``Random_Parameter_value``
    ``Gaussian_GEF_Param.gaussian_st_dev_for_parameter()``
    Rough estimates of GEF parameter stability limits: 
    url:<https://github.com/UPTEC-F-23065/McPUFF/blob/12b9f864240fe133ef727955c9e7715a7e0e20a9/Rough_limits_GEF_parameter_values.txt>

    Notes
    -----
    The ``TMC_Mod_Param_object`` has the same role in a ``TMC`` mode 
    simulation as a ``Random_Parameter_value`` object has in a 
    ``Single_Parameters`` mode simulation. It stores simulation 
    information about a specific parameter. The difference between them is
    that the ``Random_Parameter_value`` object stores the perturbed 
    simulation results for each parameter, whereas the 
    ``TMC_Mod_Param_object`` object does not. Since many, or all, 
    parameters are perturbed at once during a ``TMC`` simulation they all
    share the same perturbed results for the fission observables. In order
    to avoid redundancy, the perturbed results are stored only once in 
    the ``TMC_Object`` object together with all the 
    ``TMC_Mod_Param_object`` objects.

    Examples
    --------
    Instantiates a ``TMC_Object`` object:
    >>> tmc_mod_param_obj = TMC_Mod_Param_object("P_Shell_SL4",-0.4,"normal")
    [tmc_mod_param_obj (TMC_Mod_Param_object object)]
    """
    #--------------------------------------------------- CLASS ATTRIBUTES ------------------------------------------------------#
    param_name = None
    """Name of GEF parameter (`str`)."""

    unpert_param_val = None
    """Default GEF parameter value (`float`)."""

    pert_param_val = None
    """Perturbed GEF parameter value (`float`)."""

    random_value = None
    """Random value used to perturb GEF defaults parameter value (`float`)."""

    scaling_number_uniform = None
    """Scaling number for uniform distribution.
    
    Lets user adjust span for the uniform distribution. (See the 
    "See Also" section for parameter stability limits). Creates the upper
    and lower limits for the span of the uniform distribution. The 
    scaling is given as a percentage of the unperturbed parameter value. 
    E.g. scaling_number_uniform = 0.5 sets the span of the uniform 
    distribution as: [0.5*unperturbed_value, 1.5*unperturbed_value].
    """
    upper_lim = None
    """Upper limit of uniform distribution."""

    lower_lim = None
    """Lower limit of uniform distribution."""

    GEF_st_dev = None
    """Standard deviation of the Gaussian (normal) distribution is given 
    as a percentage of the unperturbed parameter value in the file: 
    ``Gaussian_GEF_Param.py`` (See the "See Also" section)."""
    #------------------------------------------------------ CONSTRUCTOR --------------------------------------------------------#
    def __init__(self, parameter_name, unperturbed_param_val,distribution_flag):
        self.param_name = str(parameter_name)
        self.unpert_param_val = np.float32(unperturbed_param_val)
        self.pert_param_val = None
        self.random_value = None
        #------------------------------------ CREATE RANDOM NUMBER FOR "UNIFORM" DISTRIBUTION ----------------------------------#
        if distribution_flag == 'uniform':
            scaling_number_uniform = np.float32(0.5)                                                                                    # Scaling value for limits of distribution. User can adjust the span of the uniform distribution here.
            upper_lim = abs(self.unpert_param_val) + scaling_number_uniform*abs(self.unpert_param_val)                                  # Upper numerical limit of distribution.
            lower_lim = abs(self.unpert_param_val) - scaling_number_uniform*abs(self.unpert_param_val)                                  # Lower numerical limit of distribution.
            self.random_value = ((np.random.default_rng().uniform(low=lower_lim,high=upper_lim,size=1)).astype(dtype=np.float32))[0]    #Important [0]. Removes [ndarray] that prohibits MyParameters from reading
            del scaling_number_uniform,upper_lim,lower_lim,distribution_flag
        #------------------------------------ CREATE RANDOM NUMBER FOR "NORMAL" DISTRIBUTION -----------------------------------#
        elif distribution_flag == 'normal':
            GEF_st_dev = abs(np.float32(Gaussian_GEF_Param.gaussian_st_dev_for_parameter(self.param_name,self.unpert_param_val)))       # Retrieval of std from external input file.
            self.random_value = (np.random.default_rng().normal(self.unpert_param_val,scale=GEF_st_dev,size=1)).astype(np.float32)[0]   # Important [0]. Goes into [ndarray] that otherwise prohibits MyParameters from reading
            del distribution_flag,GEF_st_dev
        else:
            print(f'The "{distribution_flag}"-distribution does not exist')
            sys.exit('Exiting program')
#------------------------------------------------------- END OF CLASS ------------------------------------------------------#
            
#------------------------------------------------------- END OF SCRIPT -----------------------------------------------------#