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
"""Select the GEF parameters that are to be perturbed in a simulation.

Notes
-----
The parameters in this list are taken from the GEF version 
"GEF-Y2023/V1.1". More information about the parameters can be found in
the GEF documentation and source code [1]_.

References
----------
.. [1] Schmidt K.H. GEF 2023/1.1, 2023. 
<https://www.khschmidts-nuclear-web.eu/GEF-2023-1-1.html>, 
accessed 4 January 2024.
"""
#-------------------------------------------------------- MODULE METHODS -------------------------------------------------------#
def parameters_to_vary():
    """Select GEF parameters to perturb in simulation from dictionary.
    
    This module contains all parameters found in the GEF source file
    "Parameters.bas". If a parameter is to be left unperturbed in a 
    simulation, then the line of the parameter should be commented out
    (See Examples).

    Parameters
    ----------
    None

    Returns
    -------
    parameter_list : `list` [`int`]

    Notes
    -----
    The function returns the position of the parameter in "Parameters.bas"
    represented by the line number counted from the top. The position is
    used to read the default value of the parameter in the constructor of
    the ``Reaction`` class in ``McPUFF_program.py``. 

    Parameter 55 ("_betaH0") has two values in "Parameters.bas", one which
    is commented out. The source code file has an empty line at line 98. 
    This dictionary mimics those features to facilitate comparison with 
    the source code file. The comments in the GEF file "Parameters.bas" 
    are also included here.  

    See Also
    --------
    ``McPUFF_program``

    Examples
    --------
    For this example, all but parameters 26,27,28,29 are commented out.
    
    >>> McPUFF_program.parameters_to_vary()
    [26, 27, 28, 29]
    """
    parameters = {
#Do not use-#'Emax_valid':2,       # /' Maximum allowed excitation energy '/
#Do not use-#'Eexc_min_multi':3,   #/' Threshold for calc. of multi-chance fission '/
        #'_Delta_S0':4,            #/' Shell effect for SL, for individual systems '/               
        #'EOscale':5,              #/' Scaling factor for even-odd structure in yields '/
#Do not use-#'Emode':6,            #/' 0: E over BF_B; 1: E over gs; 2: E_neutron; 12: E_proton '/ 
        #'_P_DZ_Mean_S1':8,
        #'_P_corr_S1':9,
        #'_P_DZ_Mean_S2':10,                                              
        #'_P_DZ_Mean_S3':11,       #/' Shift of mean Z of Mode 3 '/           
        #'_P_DZ_Mean_S4':12,       #/' Shell for structure at A around 190 '/
        #'ZC_Mode_4L':13,          #/' Maximum at Pu, higher ZC_Mode_4L -> towards Am etc. '/ 
        #'_P_Z_Curv_S1':14,                                                                      
        #'P_Z_Curvmod_S1':15,      #/' Scales energy-dependent shift '/ 
        #'_P_Z_Curv_S2':16,                                                   
        #'_S2leftmod':17,          #/' Asymmetry in diffuseness of S2 mass peak '/ 
        #'P_Z_Curvmod_S2':18,      #/' Scales energy-dependent shift '/
        #'_P_A_Width_S2':19,       #/' A width of Mode 2 (box) '/                                                
        #'_P_Z_Curv_S3':20,
        #'P_Z_Curvmod_S3':21,      #/' Scales energy-dependent shift '/
        #'P_Z_Curv_SL4':22,
        #'P_Z_Sigma_SL4':23,   
        #'_P_Z_Curv_S4':24,
        #'P_Z_Curvmod_S4':25,      #/' Scales energy-dependent shift '/
        #'_P_Shell_S1':26,         #/' Shell effect for Mode 1 (S1) '/             
        #'_P_Shell_S2':27,         #/' Shell effect for Mode 2 (S2) '/             
        #'_P_Shell_S3':28,         #/' Shell effect for Mode 3 (SA) '/             
        #'P_Shell_SL4':29,
        #'_P_Shell_S4':30,
        #'P_S4_mod':31,  
        #'_PZ_S3_olap_pos':32,     #/' Pos. of S3 shell in light fragment (in Z) '/
        #'_PZ_S3_olap_curv':33,                                                    
        #'ETHRESHSUPPS1':34,   
        #'ESIGSUPPS1':35,     
        #'Level_S11':36,           #/' Level for mode S11 (higher values means more symmetric fission) '/
        #'Shell_fading':37,        #/' fading of shell effect with E* '/                             
        #'_T_low_S1':38,                                                           
        #'_T_low_S2':39,           #/' Slope parameter for tunneling '/            
        #'_T_low_S3':40,           #/' Slope parameter for tunneling '/            
        #'_T_low_S4':41,           #/' Slope parameter for tunneling '/
        #'_T_low_SL':42,           #/' Slope parameter for tunneling '/            
        #'T_low_S11':43,           #/' Slope parameter for tunneling '/
        #'_P_att_pol':44,                                                          
        #'P_att_pol2':45,
        #'P_att_pol3':46,
        #'_P_att_rel':47,          #/' Relative portion of attenuation '/          
        #'_dE_Defo_S1':48,         #/' Deformation energy expense for Mode 1 '/    
        #'_dE_Defo_S2':49,         #/' Deformation energy expense for Mode 2 '/    
        #'_dE_Defo_S3':50,         #/' Deformation energy expense for Mode 3 '/      
        #'_dE_Defo_S4':51,         #/' Deformation energy expense for Mode 4 '/    
        #'_betaL0':52,
        #'_betaL1':53, 
        #'_betaH0':54,             #/' Offset for deformation of heavy fragment '/
#Do not use.    #Only in list for for completeness.   #'_betaH0':55,            #' This value in combination with EDISSFRAC = 0.35 gets nu-bar correct for U235T and Cf252SF
        #'_betaH1':56,
        #'_dbeta_S3':57,
        #'kappa':58,               #/' N/Z dedendence of A-asym. potential '/
        #'TCOLLFRAC':59,           #/' Tcoll per energy gain from saddle to scission '/
        #'_ECOLLFRAC':60,          #/' Fraction of pot. energy gain from saddle to scission, going into coll. excitations '/
        #'TFCOLL':61,  
        #'TCOLLMIN':62,
        #'ESHIFTSASCI_intr':63,    #/' Shift of saddle-scission energy '/ 
        #'ESHIFTSASCI_coll':64,    #/' Shift of saddle-scission energy '/
        #'_EDISSFRAC':65,
        #'Epot_shift':66,
        #'SIGDEFO':67,  
        #'SIGDEFO_0':68,
        #'SIGDEFO_slope':69,
        #'SIGENECK':70,            #/' Width of TXE by fluctuation of neck length '/   
        #'EexcSIGrel':71,
        #'DNECK':72,               #/' Tip distance at scission / fm '/
        #'FTRUNC50':73,            #/' Truncation near Z = 50 '/
        #'ZTRUNC50':74,            #/' Z value for truncation '/
        #'FTRUNC28':75,            #/' Truncation near Z = 28 '/
        #'ZTRUNC28':76,            #/' Z value for truncation '/
        #'ZMAX_S2':77,             #/' Maximum Z of S2 channel in light fragment '/
        #'NTRANSFEREO':78,         #/' Steps for E sorting for even-odd effect '/
        #'NTRANSFERE':79,          #/' Steps for E sorting for energy division '/
        #'Csort':80,               #/' Smoothing of energy sorting '/
        #'PZ_EO_symm':81,          #/' Even-odd effect in Z at symmetry '/
        #'PN_EO_Symm':82,          #/' Even-odd effect in N at symmetry '/
        #'R_EO_THRESH':83,         #/' Threshold for asymmetry-driven even-odd effect'/
        #'R_EO_SIGMA':84,
        #'R_EO_Max':85,    
        #'_POLARadd':86,           #/' Offset for enhanced polarization '/
        'POLARfac':87,            #/' Enhancement of polarization of liqu. drop '/
        'T_POL_RED':88,           #/' Reduction of temperature for sigma(Z) '/
        '_HOMPOL':89,             #/' hbar omega of polarization oscillation '/                      
        'ZPOL1':90,               #/' Extra charge polarization of S1 '/
        'P_n_x':91,               #/' Enhanced inverse neutron x section '/
        'Tscale':92,                                                                        
        'Econd':93,               #/' Pairing condenstation energy '/
        'Etrans':94,              #/' Matching energy (const. temp to Fermi gas level denstiy '/
        'T_orbital':95,           #/' From orbital ang. momentum '/
        '_Jscaling':96,           #/' General scaling of fragment angular momenta '/
        'Spin_odd':97,            #/' RMS Spin enhancement for odd Z '/ #/' Value of 0.4 adjusted to data. In conflict with Naik! '/                         
        'Esort_extend':99,        #/' Extension of energy range for E-sorting '/
        'Esort_slope':100,        #/' Onset of E-sorting around symmetry '/    
        'Esort_slope_S0':101     #/' Onset of E-sorting around symmetry for S0 channel '/
        } 
    parameter_list = []             # Create list and append position number in "Parameters.bas" for chosen parameters. Return list.                 
    """List of parameter positions in "Parameters.bas" (`int`)."""
    for param in parameters:
        parameter_list.append(parameters[param])
    return parameter_list
#------------------------------------------------------- END OF METHOD ---------------------------------------------------------#
                        
#------------------------------------------------------- END OF MODULE ---------------------------------------------------------#