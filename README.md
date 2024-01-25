# McPUFF 
$\underline{M}onte$ $\underline{C}arlo$ $\underline{P}ropagation$ $of$ $\underline{U}ncertainties$ $in$ $\underline{F}ission$ $\underline{F}ragments$

## Description


GEF [[1]](#1) and TALYS [[2]](#2) are two computer software programs used to simulate nuclear fission. There is an option called ' fymodel 4 (Okumura) ' in TALYS which makes use of GEF as a fission fragment generator and then simulates the evaporation process using a Hauser-Feshbach method. McPUFF was written to enable the user to perform simulations using the ' fymodel 4 ' option in TALYS but with perturbed parameter values in GEF. The perturbation of the parameter values in GEF, using the built-in function ' MyParameters', introduces uncertainties in the fission observables, which will be passed on to TALYS via the input data. This provides a measure of the sensitivity of the TALYS model to uncertainties in nuclear data. McPUFF performs its simulations using the ' Total Monte Carlo (TMC) ' method, described in reference [[3]](#3), which is a method for handling the propagation of uncertainties in calculations. The implementation of the TMC method requires modifications of the GEF and TALYS softwares, and these modifications are discussed in the ' Overview of the modifications of GEF and TALYS ' section. When performing a simulation using McPUFF, the user provides the necessary information about the simulated fission reaction and can choose the GEF and TALYS input data, which parameters to perturb, the magnitude of the perturbation and which distribution to draw perturbed parameter values from by making choices in a set of separate input files. All results and information about a simulation with McPUFF are stored in a python ' pickle ' file. After a simulation is completed, the pickle file can be loaded and the original object structure of McPUFF can be navigated in order to analyze the results. 

## Table of contents


- Requirements
- Overview of the GEF and TALYS modifications
- Installation
- McPUFF simulation flow chart
- McPUFF object structure
- Features
- License
- References

## Requirements


- numPy
- Matplotlib
- Modified version of TALYS 1.96. TALYS is available for download at available at [TALYS](https://tendl.web.psi.ch/tendl_2021/talys.html "TALYS software at IAEA"). 
- Modified version of GEF version 2023/1.1. GEF is available for download at [GEF](https://www.khschmidts-nuclear-web.eu/GEF.html "GEF software homepage").

## Overview of the necessary modifications of GEF and TALYS 


- Only one GEF source file needs to be modified, ' ReadParameters.mac '. The list of parameters needs to be extended to encompass all parameters in the GEF source file ' parameters.bas '. This repository contains a modified copy of the GEF source code file [ReadParameters.mac](https://github.com/UPTEC-F-23065/McPUFF/blob/0aba22b49c58036d0ab31036d52f9ad9972be772/ReadParameters.mac "Modified copy of GEF source code file."), and a pdf-file showing the necessary [GEF Modifications](https://github.com/UPTEC-F-23065/McPUFF/blob/a0ae153531fa13ed95af9f5d99c407b0fbdb05f6/Modifications_of_GEF_for_TMC_simulations.pdf "PDF with necessary modifications highlighted in yellow") if the user wants to implement the modifications on their own.
- The necessary modifications for TALYS are described in detail in the GitHub repository [Modification-of-TALYS-for-TMC-simulations ](https://github.com/UPTEC-F-23065/Modification-of-TALYS-for-TMC-simulations.git "Repository with modified copies of TALYS source files and information about necessary modifications").

## Installation


McPUFF works as a stand-alone script in Linux. To run McPUFF there are a few things that needs to be taken care of manually: 
- Modified versions of the GEF and TALYS softwares have to be compiled and be ' on path ' (for instance using the ' export PATH=.. ' command).
- Local paths to the ' McPUFF_source_files ' folder and the modified versions of GEF and TALYS have to be entered in the main script.
- There is a specific folder structure required for the program to work. When placing the McPUFF files and folders on the local system, they should be organized    as described below. (The ' Output_McPUFF ' folder and contents are included for informational purposes. They are created automatically by McPUFF).  

<pre>
McPUFF                                    # Can be named anything.
├── McPUFF_source_files
│   ├── McPUFF.py                         # Main script. Lets user set data for which fission reaction to simulate.
│   └── package_McPUFF
│       ├── GEF_input.py                  # Lets user choose GEF input data.
│       ├── Gaussian_GEF_Param.py         # Lets user set standard deviation for Gaussian distribution used for perturbed parameter values.
│       ├── McPUFF_Perturbed_Data.py      # Loads saved simulation data from pickle file and reconstructs original object structure for data analysis.
│       ├── McPUFF_program.py             # Holds all McPUFF classes and functions.
│       ├── Parameters_to_vary.py         # Lets user choose which GEF parameters to perturb in simulation. Those not chosen retain default value.
│       ├── TALYS_Input.py                # Lets user choose TALYS input data.
│       └── __init__.py                   # Includes files in main script.
└── Output_McPUFF                         # Simulation results. This folder and its contents are created automatically by McPUFF.
    ├── GEF_working_directory             # Holds GEF files and folders during McPUFF execution. Contents deleted after each loop. Can be saved for analysis.
    ├── TALYS_working_directory           # Holds TALYS files and folders during McPUFF execution. Contents deleted after each loop. Can be saved for analysis.
    └── pickle_results                    # Saves all fission reaction and results in object structure for later analysis. Created at end of McPUFF simulation.
        └── Z92_A236_n_E2.53e-08MeV.pkl   # Example of pickle file.
</pre>

## McPUFF simulation flow chart


A TMC simulation has two sources of uncertainties, one statistical and one systematical. The systematical uncertainty is a measure of the model response to uncertainties in the input data, which is the desired information. In order to reduce the statistical uncertainty, a large number (potentially thousands) of simulations are performed and hence the McPUFF program performs parallel, multi-threaded, simulations in loops until the desired number of simulations are performed. The user provides the information about the simulated fission reaction at the outset and the rest of the simulation is automated. The figure below shows a schematic of the McPUFF simulations.

<img src=FlowChartMcPUFF.png width="464" height="537" />

## McPUFF object structure


In order to keep the results from potentially thousands of simulations separated and structured, McPUFF stores all information about a fission reaction and results from a simulation in a set of objects. Each simulation with GEF and TALYS is performed in its own thread with a specific perturbed parameter value and the results are stored in a unique object. The purpose for this is to be able to relate a specific perturbed parameter value to a specific GEF and/or TALYS result. The McPUFF program has two simulation modes, a "Single_Parameters" mode and a "TMC" mode. The object structure is different depending on which mode that is used. The object structure for the two program modes are shown in the figure below.

<img src=ObjectStructureMcPUFF.png width="498" height="503" />

## Features


The TMC method requires a large number of simulations to be performed in order to reduce the statistical error and in order to decrease the computation time, McPUFF performs multi-thread calculations in parallel. Results are stored in an object structure and, in order to free up disc space, files and folders used for the calculation are deleted after each calculation is complete.

Saves data for param vals etc

## License


[License](https://github.com/UPTEC-F-23065/Modification_of_TALYS_for_TMC_simulation/blob/0e362615d513a9d40d9e6bad77ce465fc0009aed/LICENSE)

## References


<a id="1">[1]</a>
K.-H. Schmidt, B. Jurado, C. Amouroux, and C. Schmitt. General Description of Fission
Observables: GEF Model Code. Nuclear data sheets, 131:107–221, 2016.

<a id="2">[2]</a>
Arjan Koning, Stephane Hilaire, and Stephane Goriely. TALYS: modeling of nuclear reactions.
The European physical journal. A, Hadrons and nuclei, 59(6), 2023.

<a id="3">[3]</a> 
P. Karlsson. Total Monte Carlo of the fission model in GEF
and its influence on the nuclear evaporation in TALYS.
Technical Report UPTEC F 23065, Uppsala university,
2023, [Link](http://urn.kb.se/resolve?urn=urn:nbn:se:uu:diva-517598 "Thesis work at DIVA"). 
Accessed 1 January 2024.
