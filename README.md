# A123_26650_StateOfPowerEstimation_LPV_MPC

This repository contains the linear parameter varying (LPV) MPC-based lithium-ion battery state of power estimation algorithm for the A123 26650 m1b cell developed and used in the following publication

<a href="https://mountainscholar.org/handle/10976/167269">Kawakita de Souza, A. (2020). Advanced Predictive Control Strategies for Lithium-Ion Battery Management Using a Coupled Electro-Thermal Model [Master thesis, University of Colorado, Colorado Springs]. ProQuest Dissertations Publishing.</a>

- This software uses an LPV inspired Model Predictive Control(MPC) to computed the maximum discharge/charge power limit while enforcing constraints on current, voltage and temperature. 
- The underlying model of this MPC-based algorithm is a coupled electro-thermal (CET) model developed in thesis above. The details of the CET model is presented below. The CET model can be parameterized using the <a href="https://data.mendeley.com/datasets/p8kf893yv3/1">A123 26650 dataset</a> . For details about the model parameterization see the <a href="https://mountainscholar.org/handle/10976/167269">thesis</a>.<br/>
- mainSOP.m is the main file to run the MPC algorithm
