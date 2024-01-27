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
"""Initialize a `McPUFF_Perturbed_Data` object by deserializing
a `pickle` file with simulation data."""
from package_McPUFF import McPUFF_program as main_program
from matplotlib import pyplot as plt
import numpy as np
import pickle
import sys
import os
class McPUFF_Perturbed_Data(main_program.TMC_Mod_Param_object):
    """Recreate a `Reaction` object from a `pickle` file by initializing 
    the instance attributes of a `McPUFF_Perturbed_Data` object.

    If several `pickle` files are present in the target folder, objects
    from each `pickle` file are appended to the list of objects. This 
    enables simulations to be made in batches if there are time 
    constraints for the simulations. This is done differently for the
    McPUFF ``Single_Parameters`` and ``TMC`` modes (See the 'Notes' 
    section).
    
    Parameters
    ----------
    pth_PKL : `str`
        Local path to `pickle` file(s) with serialized simulation data.
    
    See Also
    --------
    ``McPUFF_program.Reaction``
    ``McPUFF_program.Modified_Parameter``
    ``McPUFF_program.TMC_Mod_Param_object``
    ``package_McPUFF.Parameters_to_vary``
    ``McPUFF_program.GEF_FY_for_TALYS()``
    ``McPUFF_program.FY_results()``
    ``McPUFF_program.read_and_clear_GEF_results()``
    ``McPUFF_program.read_and_clear_TALYS_results()``

    Notes
    -----
    This Class inherits class ``McPUFF_program.TMC_Mod_Param_object`` from
    the main program, which inherits ``McPUFF_program.Modified_Parameter``, 
    which inherits ``McPUFF_program.Reaction``. Hence this class gets 
    access to all class and instance attributes of the main program. 
    For clarity, this class does not call the constructors of the classes
    in the main program, but initializes the instance attributes 
    explicitly.

    The `pickle` file(s) needs to be initialized differently depending on 
    whether the ``Single_Parameters`` or ``TMC`` mode is used. When several 
    `pickle` files from McPUFF simulations with the ``Single_Parameters` 
    mode are loaded, in order to keep the results separated, the results 
    from each `pickle` file has to be stored as a ``Reaction`` object in 
    the ``list_of_Reac_objects``-list. When `pickle` files from McPUFF 
    simulations with the ``TMC`` mode are loaded, as long as the `pickle`
    files are for the same fission reaction there is no risk of 
    mixing the results since all relevant information for a simulation is 
    contained within its own ``TMC_Object``. The results from each `pickle`
    file can then simply be added to the ``list_of_TMC_Objects``-list.
    The mode used is determined by checking if the file name starts with a 
    specific string sequence. The user can change this string in the 
    `if`-statements in the code to suit the chosen file name. In this file
    the modes are identified with the strings 'Single' or 'TMC' at the 
    beginning. The user should rename the finished `pickle` files 
    accordingly. The reason for this setup is that TMC simulations can be 
    performed in batches if there are time constraints.

    For each simulation loop in McPUFF, following the procedure 
    developed in [1]_, fission events that only occur once are 
    removed. For more information, see [2]_. Hence, the actual 
    simulation results are fewer than specified by user. The 
    actual number of simulations in a loop can be determined using
    ``unperturbed_ignored_events``.

    References
    ----------
    .. [1] F. Nordström. Benchmark of the fission channels in TALYS. 
    Technical Report UPTEC ES 21016, Uppsala university, 2021.

    .. [2] P. Karlsson, "Total Monte Carlo of the fission model in
    GEF and its influence on the nuclear evaporation in TALYS." 
    Technical report UPTEC F 23065, Uppsala University, 2023, 
    url: <http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598>, 
    accessed 1 Janauary 2024.

    .. [3] Python Software Foundation. Built-in Functions---Python 3.12.1 
    documentation, 2024, 
    url: <https://docs.python.org/3/library/functions.html#vars>, 
    accessed 4 January 2024.

    Examples
    --------
    Creates and initializes an instance of the class using a pickle file.
    >>> MPD_object = McPUFF_Perturbed_Data("/local/path/to/pickle/file(s)/")
    [MPD_object (McPUFF_Perturbed_Data object)]
    """
    #--------------------------------------------------- CLASS ATTRIBUTES ------------------------------------------------------#
    list_of_Mod_Param_objects  = None 
    """Holds objects created in the ``Single_Parameters`` mode 
    (`list` [``McPUFF_program.Modified_Parameter``])."""
    list_of_TMC_Objects = None 
    """Holds objects created in the ``TMC`` mode 
    (`list` [``McPUFF_program.TMC_Object``])."""
    list_of_Reac_objects = None 
    """Holds ``Reaction`` objects from ``Single_Parameters`` simulations
    when using data from multiple `pickle` files 
    (`list` [``McPUFF_program.Reaction``])."""
    reaction_info = None 
    """Holds information about:``Mc_runs``,``Z_target``,``A_compound``
    ,``E_reaction``,``num_of_rand`` 
    (`dict` [`int`,`int`,`int`,`float`, `int`])."""
    dictionary_unpert_param_name_val = None 
    """Dictionary of parameter name(s) and default parameter value(s) 
    (`dict` [`str`, `float`]).

    Parameters to perturb are chosen by user in 
    ``package_McPUFF.Parameters_to_vary``."""
    program_flag = None 
    """Flag for determining which program mode to use in simulations 
    (`str`)."""
    With_TALYS_flag = None
    """`True` for including simulation of evaporation process with TALYS
    (`boolean`)."""
    unperturbed_FY = None 
    """dtype = `float32`. Holds unperturbed fission fragment yield GEF
    results.  
    
    For GEF 'lmd'-option:  (`numpy.ndarray`, (300,11)). 
    Holds: Z, A, Yff, avg(TXE), avg(TKE), E* and std(E*) for light and 
    heavy fission fragment per fission event.
    
    For GEF 'lmd+'-option: (`numpy.ndarray`, (300,19)).
    Holds: Z, A, Yff, avg(TXE), avg(TKE), E* and std(E*), E(n), std(E(n)), 
    E(g), std(E(g)) for light and heavy fission fragment per fission event.
     
    For more information, see the `See Also` section and [2]_."""
    unperturbed_ignored_events = None
    """Fission events in GEF 'lmd' file that are removed from results 
    (`int`).

    See `Notes` for more information."""
    unperturbed_GEF_results = None 
    """The keys of this dictionary are subject to change, use the Python
    command `vars(object)` [3]_ to see the contents (`dict`). 

    See also ``McPUFF_program.read_and_clear_GEF_results()``."""
    unperturbed_TALYS_results = None 
    """The keys of this dictionary are subject to change, use the Python
    command `vars(object)` [3]_ to see the contents (`dict`). 

    See also ``McPUFF_program.read_and_clear_TALYS_results()``."""
    #------------------------------------------------------ CONSTRUCTOR --------------------------------------------------------#
    def __init__(self,pth_PKL):
        self.list_of_Mod_Param_objects = []    
        self.list_of_TMC_Objects = []
        self.list_of_Reac_objects = [] 
        for file_name in os.scandir(pth_PKL):
            #--------------------------------------------------------------------------------#
            #               For pickle files from simulations with 'TMC' mode                #
            #--------------------------------------------------------------------------------#
            if file_name.is_file() and file_name.name.startswith('TMC'): 
                print(f'File name: {file_name.name}')
                with open(str(os.path.join(pth_PKL,file_name)),'rb') as f:
                    reac_obj = pickle.load(f)
                if isinstance(reac_obj, main_program.Reaction):
                    # All serialized data is retrieved as `str` type. Conversion must be done when using data. For data type, see main program McPUFF_program.
                    self.reaction_info = getattr(reac_obj, 'reaction_info')                                         # Same for all pickle files 
                    self.dictionary_unpert_param_name_val = getattr(reac_obj, 'dictionary_unpert_param_name_val')   # Same for all pickle files
                    self.list_of_Mod_Param_objects.extend(getattr(reac_obj, 'list_of_Mod_Param_objects'))           # Extend list with obj from each pickle file
                    self.list_of_TMC_Objects.extend(getattr(reac_obj,'list_of_TMC_Objects'))                        # Extend list with obj from each pickle file
                    self.program_flag = getattr(reac_obj, 'program_flag')                                           # Same for all pickle files
                    self.With_TALYS_flag = getattr(reac_obj, 'With_TALYS_flag')                                     # Same for all pickle files   
                    self.unperturbed_FY = getattr(reac_obj, 'unperturbed_FY')                                       # Same for all pickle files
                    self.unperturbed_ignored_events = getattr(reac_obj,'unperturbed_ignored_events')                # One value set because one run
                    self.unperturbed_GEF_results = getattr(reac_obj, 'unperturbed_GEF_results')                     # Same for all pickle files
                    self.unperturbed_TALYS_results = getattr(reac_obj, 'unperturbed_TALYS_results')                 # Same for all pickle files
                    del f, reac_obj
            #-------------------------------------------------------------------------------------------------------------------#
            #                           For pickle files from simulations with 'Single_Parameter' mode                          #
            #-------------------------------------------------------------------------------------------------------------------#
            elif file_name.is_file() and file_name.name.startswith('Single'): 
                print(f'File name: {file_name.name}')
                with open(str(os.path.join(pth_PKL,file_name)),'rb') as f:
                    reac_obj = pickle.load(f)
                if isinstance(reac_obj, main_program.Reaction):
                    self.list_of_Reac_objects.append(reac_obj)
                    del f, reac_obj 
            else:
                print('Incorrect input for Perturbed_Parameter_Data. No such pickle file can be found')
                sys.exit()
    #--------------------------------------------------- CLASS METHODS ---------------------------------------------------------#
                
    #***************************************************************************************************************************#    
    #                                                       PLOTS                                                               #
    #***************************************************************************************************************************#  
                           
    #======================================================================================================================================================#
    #               EXAMPLE OF LOADING A PICKLE FILE FOR 'SINGLE_PARAMETERS' MODE AND CALCULATING THE CORRELATION COEFFICIENTS                             #
    #               FOR DIVISION OF EXCITATION ENERGY BETWEEN THE LIGHT AND HEAVY OBJECT FOR VARIOUS FISSION OBSERVABLES.                                  #                            
    #======================================================================================================================================================#  
                          
    def Correlation_division_Eexc_light_heavy_Single_Parameters(MPD_obj):
        """Example of how to load data from a `pickle` file after a McPUFF
        ``Single_Parameters`` mode simulation.

        Parameters
        ----------
        MPD_obj : `object` [``McPUFF_Perturbed_Data``]
            A McPUFF_Perturbed_Data object instantiated from a ``McPUFF``
            simulation `pickle` file.

        Returns
        -------
        Function has no return value.

        See Also
        --------
        ``McPUFF_data_analysis``
        McPUFF_object_structure: url: <https://github.com/UPTEC-F-23065/McPUFF/blob/ac8c29ce13729ff85d598f8165555891a9e56ea4/ObjectStructureMcPUFF.png>
        
        Notes
        -----
        This function serves as an example of how to navigate the McPUFF
        object structure for the ``Single_Parameters`` mode. During 
        exucution, McPUFF stores all data in an object structure in order
        to keep the results from the asynchronous, parallel multi-thread 
        simulations separated. The data can be retrieved at a later stage
        by entering the lists in the objects. For more information about
        the object structure, see the 'McPUFF_object_structure' in the 
        `See Also` section.

        Examples
        --------
        This function prints a .txt file with correlation coefficients for 
        various ratios and plots a 'least-squares-fit' of the simulation
        results for a chosen parameter.

        >>> MPD.Correlation_division_Eexc_light_heavy_Single_Parameters(MPD_object)
        []
        """
        #----------------------------- CREATE LIST WITH ALL MODIFIED PARAMETER OBJECTS ----------------------------#
        Mod_param_obj_list = []                                                                                                         # Collect 'Modified_Parameter' objects from all simulations in one list.
        list_param_names = []                                                                                                           # List of selected GEF parameter names for indexing.
        #----------------------------------------------------------------------------------------------------------#
        #   Loop through all 'Reaction' objects in 'list_of_reac_objects'.                                         #
        #       Loop through all 'Modified_Parameter' objects in each 'Reaction' object.                           #
        #           Loop through all 'Random_Parameter_value' objects in each 'Modified_Parameter' object.         #
        #----------------------------------------------------------------------------------------------------------#
        #------------------------ LOOP THROUGH ALL 'REACTION' OBJECTS IN 'LIST_OF_REAC_OBJECTS' -------------------#
        for reac_obj in MPD_obj.list_of_Reac_objects:                                                                                   # Read data from each 'Reaction' object.
            if isinstance(reac_obj,main_program.Reaction):
                #------------ LOOP THROUGH ALL 'MODIFIED_PARAMETER' OBJECTS IN EACH 'REACTION' OBJECT -------------#
                for mod_param_obj in reac_obj.list_of_Mod_Param_objects:                                                                # Read data from each 'Modified_Parameter' object
                    if isinstance(mod_param_obj,main_program.Modified_Parameter):
                        Mod_param_obj_list.append(mod_param_obj)
                        list_param_names.append(mod_param_obj.param_name)
        #----------------------------- CREATE EMPTY DICTIONARY WITH PARAMETER NAMES -------------------------------#
        param_data_dict = {key: {} for key in list_param_names}                                                                         # Holds simulation info for each modified parameter.
        correlation_data_dict = {key: {} for key in list_param_names}                                                                   # Holds correlation info for each modified parameter.
        #----------------------------- CREATE VECTOR WITH VALUES FOR RATIOS ---------------------------------------#
        #   Columns: [0]=pert_param_val,[1]=E1/E2,[2]=TXE/TKE,[3]=Q_bar,[4]=E1(n)/E2(n),[5]=E1(g)/E2(g)
        for modified_param in Mod_param_obj_list: 
            if isinstance(modified_param,main_program.Modified_Parameter):
                #----------------------------------- CREATE EMPTY DATA ARRAY --------------------------------------#
                number_of_randoms = len(modified_param.list_of_Random_Parameter_Value_objects)
                param_data_vec = np.zeros((number_of_randoms,6))                                                                        # Holds simulation info for each random number.
                #----- LOOP THROUGH ALL 'RANDOM_PARAMETER_VALUE' OBJECTS IN EACH 'MODIFIED_PARAMETER' OBJECT ------#
                for n,rand_param_obj in enumerate(modified_param.list_of_Random_Parameter_Value_objects):                               # Read data from each 'Random_Parameter_value' object.
                    if isinstance(rand_param_obj,main_program.Random_Parameter_value): 
                        number_of_fission_events = len(rand_param_obj.perturbed_FY[np.nonzero(rand_param_obj.perturbed_FY[:,0]),0][0])  # Number of nonzero elements in array.
                        #------------------------------------- E1/E2 ----------------------------------------------#
                        E1 = np.mean(rand_param_obj.perturbed_FY[0:number_of_fission_events,7][0])   
                        E2 = np.mean(rand_param_obj.perturbed_FY[0:number_of_fission_events,9][0])
                        #------------------------------------ TXE/TKE ---------------------------------------------#
                        TXE   = float(rand_param_obj.perturbed_GEF_results['mean_value_TXE'][0])                                        # Name of dictionary key can be looked up in McPUFF_program.read_and_clear_GEF_results()
                        TKE   = float(rand_param_obj.perturbed_GEF_results['mean_value_TKE_pre'][0]) 
                        #---------------------------------- avg(Q_bar) --------------------------------------------#     
                        Q_bar = float(rand_param_obj.perturbed_GEF_results['mean_value_Q_bar'][0])     
                        #---------------------------------- E1(n)/E2(n) -------------------------------------------#                                                                    
                        E1_n  = np.mean(rand_param_obj.perturbed_FY[0:number_of_fission_events,11])
                        E2_n  = np.mean(rand_param_obj.perturbed_FY[0:number_of_fission_events,13])
                        #---------------------------------- E1(g)/E2(g) -------------------------------------------# 
                        E1_g  = np.mean(rand_param_obj.perturbed_FY[0:number_of_fission_events,15])
                        E2_g  = np.mean(rand_param_obj.perturbed_FY[0:number_of_fission_events,17])                                                                                                                   
                        #-------------------------- ASSIGN DATA TO DATA VECTOR ------------------------------------#
                        param_data_vec[n,0] = float(rand_param_obj.pert_param_val)  # Perturbed parameter value
                        param_data_vec[n,1] = E1/E2                                 # average E1/E2
                        param_data_vec[n,2] = TXE/TKE                               # average TXE/TKE
                        param_data_vec[n,3] = Q_bar                                 # average Q_bar
                        param_data_vec[n,4] = E1_n/E2_n                             # average E1(n)/E2(n)
                        param_data_vec[n,5] = E1_g/E2_g                             # average E1(g)/E2(g) 
                #---------------------------------- END RANDOM PARAMETER LOOP -------------------------------------#
                #---------------------------------- ASSIGN DATA TO DICTIONARY -------------------------------------#
                param_data_dict[f'{rand_param_obj.param_name}'] = param_data_vec
                #----------------------------- CALCULATE CORRELATION COEFFICIENTS ---------------------------------#
                correlation_data_vec = np.zeros((5,2,2))                                                                                # Multidimensional array. Correlation data shape is 2x2-matrix. Array is 5 x (2x2).
                corr_E1_E2   = np.corrcoef(param_data_vec[:,0],param_data_vec[:,1])
                corr_TXE_TKE = np.corrcoef(param_data_vec[:,0],param_data_vec[:,2])
                corr_Q_bar   = np.corrcoef(param_data_vec[:,0],param_data_vec[:,3])
                corr_E1n_E2n = np.corrcoef(param_data_vec[:,0],param_data_vec[:,4])
                corr_E1g_E2g = np.corrcoef(param_data_vec[:,0],param_data_vec[:,5])
                correlation_data_vec[0,:,:] = corr_E1_E2
                correlation_data_vec[1,:,:] = corr_TXE_TKE
                correlation_data_vec[2,:,:] = corr_Q_bar
                correlation_data_vec[3,:,:] = corr_E1n_E2n
                correlation_data_vec[4,:,:] = corr_E1g_E2g
                correlation_data_dict[f'{rand_param_obj.param_name}'] = correlation_data_vec
                del E1,E2,TXE,TKE,E1_n,E2_n,E1_g,E2_g,n,correlation_data_vec,corr_E1_E2,corr_TXE_TKE,corr_Q_bar,corr_E1n_E2n,corr_E1g_E2g,param_data_vec
        #----------------------------------------- PRINT RESULTS TO FILE ------------------------------------------#                 
        with open("/__local_path__/correlation_coefficients.txt", 'w') as outputfile:                                               # Add local path to save position.
            outputfile.write('Parameter name'.ljust(19," ")+'E1/E2'.ljust(12," ")+'TXE/TKE'.ljust(12," ")+'Q_bar'.ljust(10," ")+'E1(n)/E2(n)'.ljust(15," ")+'E1(g)/E2(g)'+'\n')
            outputfile.write(78*'-'+'\n')
            for para_name in list_param_names:
                string = str(para_name).ljust(18," ")
                for i in range(5):
                    if i == 3:
                        distance = 8
                    else:
                        distance = 5
                    if correlation_data_dict[para_name][i,0,1]<0:
                        string += str('{:05.4f}'.format(correlation_data_dict[para_name][i,0,1]))+distance*' '
                    else:
                        string += ' '+str('{:05.4f}'.format(correlation_data_dict[para_name][i,0,1]))+distance*' '
                outputfile.write(string+'\n')
        #--------------------------------------------- PLOT RESULTS -----------------------------------------------# 
        #--------------------------------- Plot scatter data with regression line  --------------------------------#      
        fig, ax = plt.subplots(figsize = (12,8))
        x_data = param_data_dict['_Delta_S0'][:,0]                                                                                 # Add name of parameter to plot.
        y_data = param_data_dict['_Delta_S0'][:,4]                                                                                 # Add name of parameter to plot.
        ax.scatter(x_data, y_data, s=60, alpha=0.7, edgecolors="k",label='E1(n)/E2(n) for parameter')                              # Add label text.
        b, a = np.polyfit(x_data, y_data, deg=1)
        xseq = np.linspace(min(x_data), max(x_data), num=100)
        ax.plot(xseq, a + b * xseq, color="r",ls='--', lw=2.5,label='least squares fit')
        ax.set_xlabel('Parameter values'); ax.set_ylabel('average E1(n)/E2(n)')
        ax.legend()
        plt.title('Least squares fit avg(E1(n)/E2(n)) vs parameter 100 sims')                                                      # Add title text.
        plt.savefig('/__local_path__/Scatterplot_parameter.png',dpi=300,format='png')                                              # Add local path to save position.
        plt.show()               
#------------------------------------------------------- END OF METHOD --------------------------------------------#
#------------------------------------------------------- END OF EXAMPLE -------------------------------------------#
        
#======================================================================================================================================================#
#               EXAMPLE OF LOADING A PICKLE FILE FOR 'TMC' MODE AND CALCULATING THE CORRELATION COEFFICIENTS                                           #
#               FOR DIVISION OF EXCITATION ENERGY BETWEEN THE LIGHT AND HEAVY OBJECT FOR VARIOUS FISSION OBSERVABLES.                                  #                            
#======================================================================================================================================================#  
                        
    def Correlation_division_Eexc_light_heavy_TMC(MPD_obj):
        """Example of how to load data from a `pickle` file after a McPUFF
        ``TMC`` mode simulation.

        Parameters
        ----------
        MPD_obj : `object` [``McPUFF_Perturbed_Data``]
            A McPUFF_Perturbed_Data object instantiated from a ``McPUFF``
            simulation `pickle` file.

        Returns
        -------
        Function has no return value.

        See Also
        --------
        ``McPUFF_data_analysis``
        McPUFF_object_structure: url: <https://github.com/UPTEC-F-23065/McPUFF/blob/ac8c29ce13729ff85d598f8165555891a9e56ea4/ObjectStructureMcPUFF.png>
        
        Notes
        -----
        This function serves as an example of how to navigate the McPUFF
        object structure for the ``TMC`` mode. During 
        exucution, McPUFF stores all data in an object structure in order
        to keep the results from the asynchronous, parallel multi-thread 
        simulations separated. The data can be retrieved at a later stage
        by entering the lists in the objects. For more information about
        the object structure, see the 'McPUFF_object_structure' in the 
        `See Also` section.

        The results of this example should not be taken at face value. 
        During a ``TMC`` simulation, all selected parameters are perturbed 
        at once and influence the shared results for all parameters. Hence, 
        if one parameter inflicts a change, an unrelated parameter may also 
        show a correlation. The example is presented like this to aid in a 
        comparison with the ``Single_Parameters`` example when it comes to
        the navigation of the object structure.

        Examples
        --------
        This function prints a .txt file with correlation coefficients for 
        various ratios and plots a 'least-squares-fit' of the simulation
        results for a chosen parameter.

        >>> MPD.Correlation_division_Eexc_light_heavy_TMC(MPD_object)
        []
        """
        #----------------------------- CREATE LIST WITH ALL "TMC_OBJECT" OBJECTS ----------------------------------#
        #   When using the "TMC" mode, all of the selected GEF parameters are individually perturbed at once.      #
        #   The "TMC_Object" holds all information and results from one loop of McPUFF, i.e one perturbed run.     #
        #   All perturbed parameters share the same simulation result for each perturbed run, which is stored      #
        #   in the 'TMC_Object'. But each individual perturbed parameter value is stored in a unique               #
        #   'TMC_Mod_Param_object' object for each run.                                                            #
        #----------------------------------------------------------------------------------------------------------#
        #----------------------------------------------------------------------------------------------------------#
        #   Loop through all 'TMC_Object' objects in 'list_of_TMC_Objects'.                                        #
        #       Loop through all 'TMC_Mod_Param_object' objects in each 'TMC_Object' object.                       #
        #----------------------------------------------------------------------------------------------------------#

        #----------------------------- CREATE EMPTY DICTIONARY WITH PARAMETER NAMES -------------------------------#
        if isinstance(MPD_obj.list_of_TMC_Objects[0],main_program.TMC_Object):
            param_data_dict = {key: [] for key in MPD_obj.list_of_TMC_Objects[0].dictionary_unpert_param_name_val}                      # Holds simulation info for each modified parameter.
            correlation_data_dict = {key: {} for key in MPD_obj.list_of_TMC_Objects[0].dictionary_unpert_param_name_val}                # Holds correlation info for each modified parameter.
        # ----------------------- LOOP THROUGH ALL 'TMC_OBJECTS' IN 'LIST_OF_TMC_OBJECTS' -------------------------#
        for k,tmc_obj in enumerate(MPD_obj.list_of_TMC_Objects):
            if isinstance(tmc_obj,main_program.TMC_Object):
                #------------ CREATE VECTOR WITH VALUES FOR RATIOS FOR EACH TMC_OBJECT (EACH PERTURBED RUN) -------#
                #   Columns: [0]=pert_param_val,[1]=E1/E2,[2]=TXE/TKE,[3]=Q_bar,[4]=E1(n)/E2(n),[5]=E1(g)/E2(g)
                number_of_fission_events = len(tmc_obj.perturbed_FY[np.nonzero(tmc_obj.perturbed_FY[:,0]),0][0])                        # Number of nonzero elements in array.
                #------------------------------------- E1/E2 ------------------------------------------------------#
                E1 = np.mean(tmc_obj.perturbed_FY[0:number_of_fission_events,7][0])   
                E2 = np.mean(tmc_obj.perturbed_FY[0:number_of_fission_events,9][0])
                #------------------------------------ TXE/TKE -----------------------------------------------------#
                TXE   = float(tmc_obj.perturbed_GEF_results['mean_value_TXE'][0])                                                       # Name of dictionary key can be looked up in McPUFF_program.read_and_clear_GEF_results()
                TKE   = float(tmc_obj.perturbed_GEF_results['mean_value_TKE_pre'][0]) 
                #---------------------------------- avg(Q_bar) ----------------------------------------------------#     
                Q_bar = float(tmc_obj.perturbed_GEF_results['mean_value_Q_bar'][0])     
                #---------------------------------- E1(n)/E2(n) ---------------------------------------------------#                                                                    
                E1_n  = np.mean(tmc_obj.perturbed_FY[0:number_of_fission_events,11])
                E2_n  = np.mean(tmc_obj.perturbed_FY[0:number_of_fission_events,13])
                #---------------------------------- E1(g)/E2(g) ---------------------------------------------------# 
                E1_g  = np.mean(tmc_obj.perturbed_FY[0:number_of_fission_events,15])
                E2_g  = np.mean(tmc_obj.perturbed_FY[0:number_of_fission_events,17])      
                #--- LOOP THROUGH ALL 'TMC_MOD_PARAM' OBJECTS IN 'LIST_OF_TMC_MOD_PARAM_OBJECT' OBJECTS -----------#                                                                                                             
                for tmc_mod_param_obj in tmc_obj.list_of_TMC_Mod_Param_objects:
                    if isinstance(tmc_mod_param_obj,main_program.TMC_Mod_Param_object):
                        #-------------------- CREATE EMPTY DATA ARRAY ---------------------------------------------#
                        TMC_Mod_Param_data_vec = np.zeros((1,6))
                        #-------------------------- ASSIGN DATA TO DATA VECTOR ------------------------------------#
                        TMC_Mod_Param_data_vec[0,0] = float(tmc_mod_param_obj.pert_param_val)                                           # Perturbed parameter value
                        TMC_Mod_Param_data_vec[0,1] = E1/E2                                                                             # average E1/E2
                        TMC_Mod_Param_data_vec[0,2] = TXE/TKE                                                                           # average TXE/TKE
                        TMC_Mod_Param_data_vec[0,3] = Q_bar                                                                             # average Q_bar
                        TMC_Mod_Param_data_vec[0,4] = E1_n/E2_n                                                                         # average E1(n)/E2(n)
                        TMC_Mod_Param_data_vec[0,5] = E1_g/E2_g                                                                         # average E1(g)/E2(g) 
                        #-------------------------- ASSIGN DATA TO DICTIONARY -------------------------------------#
                        param_data_dict[f'{tmc_mod_param_obj.param_name}'].append(TMC_Mod_Param_data_vec)
                        del TMC_Mod_Param_data_vec
                        #------------------- END 'TMC_MOD_PARAM_OBJECT' OBJECTS LOOP ------------------------------#
                #---------------------------------- END 'TMC_OBJECT' LOOP -----------------------------------------#
        #------------------------------------ CALCULATE CORRELATION COEFFICIENTS ----------------------------------#
        for parameter_name in param_data_dict:
            stacked_data = np.vstack(param_data_dict[parameter_name])                                                                   # Stack data from arrays in list of arrays into one array.
            correlation_data_vec = np.zeros((5,2,2))                                                                                    # Multidimensional array. Correlation data shape is 2x2-matrix. Array is 5 x (2x2).
            corr_E1_E2   = np.corrcoef(stacked_data[:,0],stacked_data[:,1])
            corr_TXE_TKE = np.corrcoef(stacked_data[:,0],stacked_data[:,2])
            corr_Q_bar   = np.corrcoef(stacked_data[:,0],stacked_data[:,3])
            corr_E1n_E2n = np.corrcoef(stacked_data[:,0],stacked_data[:,4])
            corr_E1g_E2g = np.corrcoef(stacked_data[:,0],stacked_data[:,5])
            correlation_data_vec[0,:,:] = corr_E1_E2
            correlation_data_vec[1,:,:] = corr_TXE_TKE
            correlation_data_vec[2,:,:] = corr_Q_bar
            correlation_data_vec[3,:,:] = corr_E1n_E2n
            correlation_data_vec[4,:,:] = corr_E1g_E2g
            correlation_data_dict[parameter_name] = correlation_data_vec        
        del E1,E2,TXE,TKE,E1_n,E2_n,E1_g,E2_g,k,correlation_data_vec,corr_E1_E2,corr_TXE_TKE,corr_Q_bar,corr_E1n_E2n,corr_E1g_E2g
        #----------------------------------------- PRINT RESULTS TO FILE ------------------------------------------#                 
        with open("/__local_path__/correlation_coefficients.txt", 'w') as outputfile:                                                  # Add local path to save position.
            outputfile.write('Parameter name'.ljust(19," ")+'E1/E2'.ljust(12," ")+'TXE/TKE'.ljust(12," ")+'Q_bar'.ljust(10," ")+'E1(n)/E2(n)'.ljust(15," ")+'E1(g)/E2(g)'+'\n')
            outputfile.write(78*'-'+'\n')
            for para_name in param_data_dict:
                string = str(para_name).ljust(18," ")
                for i in range(5):
                    if i == 3:
                        distance = 8
                    else:
                        distance = 5
                    if correlation_data_dict[para_name][i,0,1]<0:
                        string += str('{:05.4f}'.format(correlation_data_dict[para_name][i,0,1]))+distance*' '
                    else:
                        string += ' '+str('{:05.4f}'.format(correlation_data_dict[para_name][i,0,1]))+distance*' '
                outputfile.write(string+'\n')
        #--------------------------------------------- PLOT RESULTS -----------------------------------------------# 
        #--------------------------------- Plot scatter data with regression line  --------------------------------#      
        fig, ax = plt.subplots(figsize = (12,8))        
        x_data = np.vstack(param_data_dict['T_orbital'])[:,0]                                                                           # Add name of parameter to plot.
        y_data = np.vstack(param_data_dict['T_orbital'])[:,4]                                                                           # Add name of parameter to plot.
        ax.scatter(x_data, y_data, s=60, alpha=0.7, edgecolors="k",label='E1(n)/E2(n) for parameter')                                   # Add label text.
        b, a = np.polyfit(x_data, y_data, deg=1)
        xseq = np.linspace(min(x_data), max(x_data), num=100)
        ax.plot(xseq, a + b * xseq, color="r",ls='--', lw=2.5,label='least squares fit')
        ax.set_xlabel('Parameter values'); ax.set_ylabel('average E1(n)/E2(n)')
        ax.legend()
        plt.title('Least squares fit avg(E1(n)/E2(n)) vs parameter 100 sims')                                                           # Add title text.
        plt.savefig('/__local_path__/Scatterplot_parameter.png',dpi=300,format='png')                                                  # Add local path to save position.
        plt.show()    
#------------------------------------------------------- END OF METHOD --------------------------------------------#
#------------------------------------------------------- END OF EXAMPLE -------------------------------------------#
#------------------------------------------------------- END OF CLASS ---------------------------------------------#
        