# MPC_Cancer
This repository contains the MATLAB scripts for the implementation of the Model Predictive Control for the optimisation of the cancer treatment through the combination of chemotherapy and immunotherapy.

This approach combines the linearization of the nonlinear dynamical system, which is subsequently combined with the application
of the condensed approach to the construction and solution of the Quadratic Programming problem over the length of the
prediction horizon, as described in
Jerez, J. L., Kerrigan, E. C., & Constantinides, G. A. (2011). A Condensed and Sparse QP Formulation for Predictive Control. https://doi.org/10.0/Linux-x86_64

Of course, the cancer treatment implementation can be easily amended so that different systems can be implemented. Code is well commented
and you should be able to read through that, given that you have the understanding of the basic Quadratic Programming,
linear algebra, and dynamical systems.

The dynamical system is described as follows: 
![equation](https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20%5Cfrac%7Bd%7D%7Bdt%7Dx%20%26%3D%20-%5Cmu_%7BC%7Dx%20%5Cleft%28%20ln%5Cfrac%7Bx%7D%7Bx_%7B%5Cinfty%7D%7D%20%5Cright%20%29%20-%20%5Cgamma%20xy%20-%20k_%7Bx%7Dxu%20%5C%5C%20%5Cfrac%7Bd%7D%7Bdt%7Dy%20%26%3D%20%5Cmu_%7BI%7D%28x%20-%20%5Cbeta%20x%5E%7B2%7D%29y%20-%20%5Cdelta%20y%20&plus;%20%5Calpha%20&plus;%20k_%7By%7Dyv%20%5Cend%7Balign*%7D)

The results of the simulation are displayed below: 
![alt text](https://github.com/miroslavgasparek/MPC_Cancer/blob/master/Cancer_Treatment_MPC.jpg)
