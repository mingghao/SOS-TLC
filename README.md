# SOS-TLC

This repository contains the official MATLAB implementation for our paper: **"Taylor-Lagrange Control for Sampled-Data Systems"**. 

It provides the codebase for running Sampled-Data TLC (SD-TLC) and Sum-of-Squares TLC (SOS-TLC) to ensure continuous-time safety for nonlinear systems with Zero-Order-Hold (ZOH) inputs.

## Dependencies

To run the code in this repository, you will need MATLAB along with the following optimization solvers and toolboxes:

*   **[MOSEK](https://www.mosek.com/)**: A high-performance solver for Quadratic Programs (QP) and Semidefinite Programs (SDP).
*   **[YALMIP](https://yalmip.github.io/)**: A MATLAB toolbox for modeling and solving optimization problems.
