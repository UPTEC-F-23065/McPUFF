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
"""Set the standard deviation for the Gaussian distribution used for
creating random parameter values.
"""
def gaussian_st_dev_for_parameter(parameter_name,param_val):
    """Set the standard deviation for the Gaussian distribution.
     
    The standard deviation is created as a percentage of the of the
    default parameter value. The percentage is chosen by the user by
    changing the value for the parameter `perturbation_in_percent`. 
    
    Parameters
    ----------
    parameter_name : `str`
        The parameter name is used as the key to create an entry
        in the dictionary of parameter names and values of std.
    param_val : `float`
        The default GEF parameter value. Used to calculate the std as a
        percentage of the default value. 
    
    Returns
    -------
    st_dev_param : `dict` [`str`, `float`]
        The keys are the parameter names and the values are the
        std used for the Gaussian distribution.

    See Also
    --------
    `McPUFF_program.create_perturbed_parameter_value()`

    Notes
    -----
    Some GEF parameters have a suggested std in the GEF source code.
    Those parameters that do not are commented by "No value in GEF".
    Some GEF parameters have a default value of zero. These are commented 
    by "Default parameter value in GEF is zero". The std of these are set 
    in the method `McPUFF_program.create_perturbed_parameter_value`.

    When this function is called, it sets all parameter values equal to 
    the `param_val` argument. This is redundant, but the code is 
    written in this way because the structure makes it possible, with some
    changes in the code and function calls, to set the standard deviation
    for each parameter individually, should it be desired. (The std for a 
    parameter can be hard coded).

    Examples
    --------
    (If `perturbation_in_percent` = 0.03.)

    >>> std = Gaussian_GEF_Param.gaussian_st_dev_for_parameter(name_of_parameter,5)
    >>> print(std[name_of_parameter])
    [0.15]
    """
    #--------------------------------------------------- METHOD VARIABLES ------------------------------------------------------#
    parameter_value = float(param_val)    
    """Default GEF parameter value (`float`)."""     
    perturbation_in_percent = float(0.03)                   # CHANGE THIS VALUE TO SET STANDARD DEVIATION FOR THE DISTRIBUTION.
    """Width of std for Gaussian distribution (`float`)."""
    st_dev_param = {
                    '_Delta_S0':parameter_value*perturbation_in_percent,            # Default parameter value in GEF is zero.
                    'EOscale':parameter_value*perturbation_in_percent,
                    '_P_DZ_Mean_S1':parameter_value*perturbation_in_percent,
                    '_P_corr_S1':parameter_value*perturbation_in_percent,
                    '_P_DZ_Mean_S2':parameter_value*perturbation_in_percent,
                    '_P_DZ_Mean_S3':parameter_value*perturbation_in_percent,
                    '_P_DZ_Mean_S4':parameter_value*perturbation_in_percent,
                    'ZC_Mode_4L':parameter_value*perturbation_in_percent,           # No value in GEF.
                    '_P_Z_Curv_S1':parameter_value*perturbation_in_percent,
                    'P_Z_Curvmod_S1':parameter_value*perturbation_in_percent,       # No value in GEF.
                    '_P_Z_Curv_S2':parameter_value*perturbation_in_percent,
                    '_S2leftmod':parameter_value*perturbation_in_percent,
                    'P_Z_Curvmod_S2':parameter_value*perturbation_in_percent,       # No value in GEF.
                    '_P_A_Width_S2':parameter_value*perturbation_in_percent,
                    '_P_Z_Curv_S3':parameter_value*perturbation_in_percent,
                    'P_Z_Curvmod_S3':parameter_value*perturbation_in_percent,       # No value in GEF.
                    'P_Z_Curv_SL4':parameter_value*perturbation_in_percent,         # No value in GEF.
                    'P_Z_Sigma_SL4':parameter_value*perturbation_in_percent,        # No value in GEF.
                    '_P_Z_Curv_S4':parameter_value*perturbation_in_percent,
                    'P_Z_Curvmod_S4':parameter_value*perturbation_in_percent,       # No value in GEF.
                    '_P_Shell_S1':parameter_value*perturbation_in_percent,
                    '_P_Shell_S2':parameter_value*perturbation_in_percent,
                    '_P_Shell_S3':parameter_value*perturbation_in_percent,
                    'P_Shell_SL4':parameter_value*perturbation_in_percent,          # No value in GEF.
                    '_P_Shell_S4':parameter_value*perturbation_in_percent,
                    'P_S4_mod':parameter_value*perturbation_in_percent,             # No value in GEF.
                    '_PZ_S3_olap_pos':parameter_value*perturbation_in_percent,
                    '_PZ_S3_olap_curv':parameter_value*perturbation_in_percent,
                    'ETHRESHSUPPS1':parameter_value*perturbation_in_percent,        # No value in GEF.
                    'ESIGSUPPS1':parameter_value*perturbation_in_percent,           # No value in GEF.
                    'Level_S11':parameter_value*perturbation_in_percent,            # No value in GEF.
                    'Shell_fading':parameter_value*perturbation_in_percent,         # No value in GEF.
                    '_T_low_S1':parameter_value*perturbation_in_percent,
                    '_T_low_S2':parameter_value*perturbation_in_percent,
                    '_T_low_S3':parameter_value*perturbation_in_percent,
                    '_T_low_S4':parameter_value*perturbation_in_percent,
                    '_T_low_SL':parameter_value*perturbation_in_percent,
                    'T_low_S11':parameter_value*perturbation_in_percent,            # No value in GEF.
                    '_P_att_pol':parameter_value*perturbation_in_percent,
                    'P_att_pol2':parameter_value*perturbation_in_percent,           # No value in GEF.  Default parameter value in GEF is zero.
                    'P_att_pol3':parameter_value*perturbation_in_percent,           # No value in GEF.  Default parameter value in GEF is zero.
                    '_P_att_rel':parameter_value*perturbation_in_percent,           # No value in GEF.
                    '_dE_Defo_S1':parameter_value*perturbation_in_percent,
                    '_dE_Defo_S2':parameter_value*perturbation_in_percent,
                    '_dE_Defo_S3':parameter_value*perturbation_in_percent,
                    '_dE_Defo_S4':parameter_value*perturbation_in_percent,
                    '_betaL0':parameter_value*perturbation_in_percent,
                    '_betaL1':parameter_value*perturbation_in_percent,
                    '_betaH0':parameter_value*perturbation_in_percent,
                    '_betaH1':parameter_value*perturbation_in_percent,
                    '_dbeta_S3':parameter_value*perturbation_in_percent,
                    'kappa':parameter_value*perturbation_in_percent,                # No value in GEF.  Default parameter value in GEF is zero.
                    'TCOLLFRAC':parameter_value*perturbation_in_percent,            # No value in GEF.
                    '_ECOLLFRAC':parameter_value*perturbation_in_percent,
                    'TFCOLL':parameter_value*perturbation_in_percent,               # No value in GEF.
                    'TCOLLMIN':parameter_value*perturbation_in_percent,             # No value in GEF.
                    'ESHIFTSASCI_intr':parameter_value*perturbation_in_percent,     # No value in GEF.
                    'ESHIFTSASCI_coll':parameter_value*perturbation_in_percent,     # No value in GEF.
                    '_EDISSFRAC':parameter_value*perturbation_in_percent,
                    'Epot_shift':parameter_value*perturbation_in_percent,           # No value in GEF.  Default parameter value in GEF is zero.
                    'SIGDEFO':parameter_value*perturbation_in_percent,              # No value in GEF.
                    'SIGDEFO_0':parameter_value*perturbation_in_percent,            # No value in GEF.
                    'SIGDEFO_slope':parameter_value*perturbation_in_percent,        # No value in GEF.  Default parameter value in GEF is zero.
                    'SIGENECK':parameter_value*perturbation_in_percent,             # No value in GEF.
                    'EexcSIGrel':parameter_value*perturbation_in_percent,           # No value in GEF.
                    'DNECK':parameter_value*perturbation_in_percent,                # No value in GEF.
                    'FTRUNC50':parameter_value*perturbation_in_percent,             # No value in GEF.
                    'ZTRUNC50':parameter_value*perturbation_in_percent,             # No value in GEF.
                    'FTRUNC28':parameter_value*perturbation_in_percent,             # No value in GEF.
                    'ZTRUNC28':parameter_value*perturbation_in_percent,             # No value in GEF.
                    'ZMAX_S2':parameter_value*perturbation_in_percent,              # No value in GEF.
                    'NTRANSFEREO':parameter_value*perturbation_in_percent,          # No value in GEF.
                    'NTRANSFERE':parameter_value*perturbation_in_percent,           # No value in GEF.
                    'Csort':parameter_value*perturbation_in_percent,                # No value in GEF.
                    'PZ_EO_symm':parameter_value*perturbation_in_percent,           # No value in GEF.
                    'PN_EO_Symm':parameter_value*perturbation_in_percent,           # No value in GEF.
                    'R_EO_THRESH':parameter_value*perturbation_in_percent,          # No value in GEF.
                    'R_EO_SIGMA':parameter_value*perturbation_in_percent,           # No value in GEF.
                    'R_EO_Max':parameter_value*perturbation_in_percent,             # No value in GEF.
                    '_POLARadd':parameter_value*perturbation_in_percent,
                    'POLARfac':parameter_value*perturbation_in_percent,             # No value in GEF.
                    'T_POL_RED':parameter_value*perturbation_in_percent,            # No value in GEF.
                    '_HOMPOL':parameter_value*perturbation_in_percent,                            
                    'ZPOL1':parameter_value*perturbation_in_percent,                # No value in GEF.  Default parameter value in GEF is zero.
                    'P_n_x':parameter_value*perturbation_in_percent,                # No value in GEF.  Default parameter value in GEF is zero.
                    'Tscale':parameter_value*perturbation_in_percent,               # No value in GEF.           
                    'Econd':parameter_value*perturbation_in_percent,                # No value in GEF.
                    'Etrans':parameter_value*perturbation_in_percent,               # No value in GEF.
                    'T_orbital':parameter_value*perturbation_in_percent,            # No value in GEF.  Default parameter value in GEF is zero.
                    '_Jscaling':parameter_value*perturbation_in_percent,
                    'Spin_odd':parameter_value*perturbation_in_percent,             # No value in GEF.
                    'Esort_extend':parameter_value*perturbation_in_percent,         # No value in GEF.
                    'Esort_slope':parameter_value*perturbation_in_percent,          # No value in GEF.
                    'Esort_slope_S0':parameter_value*perturbation_in_percent        # No value in GEF.
                    }
    """Dictionary of GEF parameter names and std for Gaussian distribution 
    (`dict`)."""
    del parameter_value, perturbation_in_percent
    return st_dev_param[parameter_name]
    #---------------------------------------------------------- END MODULE -----------------------------------------------------#