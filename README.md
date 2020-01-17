# Reservoir, Rock Physics, Time-lapse Seismic, & Geomechanical Modelling of Carbon Capture & Storage

![ccs_pic](https://user-images.githubusercontent.com/51282928/72623302-767cc000-3977-11ea-963b-060339c3eeba.jpg)

This is my BSc thesis which I defended in 9 August 2019. The repo includes my BSc works with few additions post-BSc that I made to enhance my methodology. 

## Introduction

This is a complete modelling workflow that integrates reservoir simulation data, rock physics, time-lapse seismic, and geomechanics of CO<sub>2</sub> injection in carbonate rock. The result of reservoir simulation is **CO<sub>2</sub> saturation** and **pressure** changes. These information are used in **rock physics modelling** to (1) characterize **pore structure** of the carbonate rock based on Xu and Payne (2009) and (2) compute seismic response (**Vp**, **Vs**, and **density**) of carbonate rock as **CO<sub>2</sub> saturation** and **pressure** increases, based on Kuster and Toksöz (1974) using the previous result of **pore structure**. After seismic response is computed, 

## About this repository

Notebooks are in the *main* folder:
* `01_reservoir_simulation.ipynb`: displays the result of reservoir simulation done in Petrel Eclipse, compute the reservoir pressure increase due to CO2 injection
* `02_rock_physics_modelling.ipynb`: computes Vp, Vs, and density change due to CO2 saturation and increasing pressure
* `03_xu_payne.ipynb`: Xu-Payne (2009) carbonate pore characterization 
* `geomechanics.m`: geomechanical modelling, built in Matlab language. Will be re-written in Python language when having leisure time :)

Folder *lib* includes all functions used in each of the notebooks.

Folder *data* includes the input data.

My BSc research had been published, see the extended abstracts in my ResearchGate here:
* [Paper 1: Time-lapse Seismic Modelling](https://www.researchgate.net/publication/338644464_Time-lapse_Seismic_Modelling_as_Tool_for_Evaluation_of_Monitoring_Feasibility_for_CO2_Sequestration_in_Gundih_Gas_Field)
* [Paper 2: Geomechanical Modelling](https://www.researchgate.net/publication/338644468_Geomechanical_Modelling_of_CO2_Sequestration_Process_in_Gundih_Field_for_Evaluation_of_Reservoir_Integrity)
