# A123_26650_StateOfPowerEstimation_LPV_MPC

This repository contains the linear parameter varying (LPV) MPC-based lithium-ion battery state of power estimation algorithm for the A123 26650 m1b cell developed and used in the following publication

<a href="https://ieeexplore.ieee.org/document/9483433">M. A. Xavier, A. K. de Souza and M. S. Trimboli, "An LPV-MPC Inspired Battery SOP Estimation Algorithm Using a Coupled Electro-Thermal Model," 2021 American Control Conference (ACC), 2021, pp. 4421-4426, doi: 10.23919/ACC50511.2021.9483433.</a>

<a href="https://www.researchgate.net/publication/345630376_A_Predictive_Modeling_and_Control_Approach_to_Improving_Lithium-ion_Battery_Performance_in_Cells_Exhibiting_Large_Voltage_Hysteresis?channel=doi&linkId=5fa96cbc458515157bf7485d&showFulltext=true">A. K. de Souza, G. Plett and M. S. Trimboli, "A Predictive Modeling and Control Approach to Improving Lithium-ion Battery Performance in Cells Exhibiting Large Voltage Hysteresis," 20th Advanced Automotive Battery Conference, 2020, DOI: 10.13140/RG.2.2.31324.92804.</a>

- This software uses an LPV inspired Model Predictive Control(MPC) to computed the maximum discharge/charge power limit while enforcing constraints on current, voltage and temperature. 
- The underlying model of this MPC-based algorithm is a coupled electro-thermal (CET) model developed in thesis above. The details of the CET model is presented below. The CET model can be parameterized using the <a href="https://data.mendeley.com/datasets/p8kf893yv3/1">A123 26650 dataset</a> . For details about the model parameterization see the <a href="https://mountainscholar.org/handle/10976/167269">thesis</a>.<br/>
- mainSOP.m is the main file to run the MPC algorithm
