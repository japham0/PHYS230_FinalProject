# PHYS 230 Final Project: Continuum Modeling for Fluid Phase-Separation

This is our submission for the final project for the Spring 2024 semester of
PHYS 230: Computational Modeling for the Biological and Interdisciplinary 
Sciences and Engineering. 


## Description

The goal of this project is to be able to model and visualize fluid phase-separation
using continuum modeling techniques. To do so, we solve the standard Cahn-Hillard equation 
using finite-differencing schemes. From there, we implement a vortex forcing term to our 
system to examine the effects of a type of fluid flow to our phase separation model. We 
also wrote a python script which calculates the average size of an oil droplet after the
phase separation, given an image of the system after Nt timesteps. 

## Usage

1. cahnhillardstandard.m models the standard cahn hillard equation phase separation over time
2. cahnhillardvortex.m models the cahn hillard equation with a vortex fluid flow added
3. convergence_test.m shows that the laplacian operator converges with order 2 accuracy
4. blobsize.ipynb prints the average size of a droplet given an image of the phase separation



