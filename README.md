# Reservoir, Rock Physics, Time-lapse Seismic, & Geomechanical Modelling of Carbon Capture & Storage

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

![ccs_pic](https://user-images.githubusercontent.com/51282928/72623302-767cc000-3977-11ea-963b-060339c3eeba.jpg)

This is my BSc thesis which I defended in 9 August 2019. The repo includes my BSc works with few additions post-BSc that I made to enhance my methodology. You are recommended to **read first my presented papers** to understand what is this about, in the link given at the end of this material. 

## Introduction

![ccs_workflow](https://user-images.githubusercontent.com/51282928/72629166-58688d00-3982-11ea-98a0-3c5004b3c796.png)

This is a complete modelling workflow that integrates reservoir simulation data, rock physics, time-lapse seismic, and geomechanics of CO<sub>2</sub> injection in carbonate rock. The result of reservoir simulation is **CO<sub>2</sub> saturation** and **pressure** changes. These information are used in **rock physics modelling** to (1) characterize **pore structure** of the carbonate rock based on Xu and Payne (2009) and (2) compute seismic response (**Vp**, **Vs**, and **density**) of carbonate rock as **CO<sub>2</sub> saturation** and **pressure** increases, based on Kuster and Toksöz (1974) using the previous result of **pore structure**. After seismic response is computed, synthetic seismic record is simulated (**seismic modelling**). However, this part is worked separately using NORSAR software, not in Python, so that is not covered in this repo. Last part is geomechanical modelling, that simulates **Mohr-Coulomb stress path** due to change in pressure. 

Case study is used from an Indonesian CCS pilot project in Gundih. Several key points of this CCS are:
* Injection is done in carbonate reservoir,
* Depth 3 km,
* High temperature (158<sup>o</sup>C), high pressure (273 bar),
* 10 years of injection,
* Rate of injection 800 ton per day.

## Content of this repository

Notebooks are in the *main* folder:
* `01_reservoir_simulation.ipynb`: displays the result of reservoir simulation done in Petrel Eclipse, compute the reservoir pressure increase due to CO2 injection
* `02_rock_physics_modelling.ipynb`: computes Vp, Vs, and density change due to CO2 saturation and increasing pressure
* `03_xu_payne.ipynb`: Xu-Payne (2009) carbonate pore characterization 
* `geomechanics.m`: geomechanical modelling, built in Matlab language. Will be re-written in Python language when having leisure time :)

Folder *lib* includes all functions used in each of the notebooks.

Folder *data* includes the input data.

## Next task and others in mind
* Re-writing geomechanical modelling from MATLAB code to Python.
* Computational modelling of increasing porosity effect to **change of porosity and pore compressibility**.
* Exploration of [PyFEHM](https://github.com/lanl/PyFEHMp), [PyTOUGH](https://github.com/acroucher/PyTOUGH), [Open Porous Media](https://github.com/OPM/opm-utilities) applied in CO<sub>2</sub> injection.

## Dependencies
[CoolProp](https://github.com/CoolProp/CoolProp): to calculate CO<sub>2</sub> property

## Read my papers

My BSc research had been presented, see the extended abstracts in my ResearchGate here:
* [Paper 1: Time-lapse Seismic Modelling](https://www.researchgate.net/publication/338644464_Time-lapse_Seismic_Modelling_as_Tool_for_Evaluation_of_Monitoring_Feasibility_for_CO2_Sequestration_in_Gundih_Gas_Field)
* [Paper 2: Geomechanical Modelling](https://www.researchgate.net/publication/338644468_Geomechanical_Modelling_of_CO2_Sequestration_Process_in_Gundih_Field_for_Evaluation_of_Reservoir_Integrity)
