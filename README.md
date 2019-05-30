# MPC_Cancer
This repository contains the MATLAB scripts for the implementation of the Model Predictive Control for the optimisation of the cancer treatment through the combination of chemotherapy and immunotherapy.

This approach combines the linearization of the nonlinear dynamical system, which is subsequently combined with the application
of the condensed approach to the construction and solution of the Quadratic Programming problem over the length of the
prediction horizon, as described in
Jerez, J. L., Kerrigan, E. C., & Constantinides, G. A. (2011). A Condensed and Sparse QP Formulation for Predictive Control. https://doi.org/10.0/Linux-x86_64

Of course, the pendulum implementation can be easily amended so that different systems can be implemented. Code is well commented
and you should be able to read through that, given that you have the understanding of the basic Quadratic Programming,
linear algebra, and dynamical systems.

The results of the simulation are displayed below: 
![alt text](https://github.com/miroslavgasparek/MPC_Cancer/blob/master/Cancer_Treatment_MPC.jpg)
