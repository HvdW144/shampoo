# Welcome to the SHAMPOO GitHub repository
SHAMPOO (**S**toc**HA**stic **M**onomer **P**r**O**cess**O**r) is a model simulating the evolution of individual dust particles in a stochastic fashion, faciliating a Monte-Carlo-like approach towards dust evolution in protoplanetary disks.

Developed by Mark Oosterloo (University of Groningen (UG), Vrije Universiteit Amsterdam (VUA)) in collaboration with Inga Kamp (UG) and Wim van Westerenen (VUA)

Please do consult the authors before publishing any results involving SHAMPOO to ensure the code has been set up appropriately.

# Introduction
Over the past years, many modelling efforts have explored the effects of dynamical, collisional and ice processing on dust in planet-forming disks. In this context, SHAMPOO provides a framework to explore the effects of these three processes combined on an individual dust particle in a disk that remains both physically and chemically static over the duration of the simulation. SHAMPOO is therefore aimed towards providing a complete description of these processes on the properties of dust particles in the planet-forming region on timescales shorter than 100 000 yr. 

# Getting started
This repository contains the data and tools to simulate and analyze the trajectories of monomers in similar fashion as in [1]. For the tools to run and analyze larger simulations containing many monomers akin to [2], we refer the user to the authors.

SHAMPOO is usually executed in a notebook fashion, but in general any python script should do. Any simulation peformed with SHAMPOO requires a background protoplanetary disk model generated with ProDiMo, along with an input file that has to end with ```shampooInput.csv```. 

The notebook ```shampooDemoIndividualMonomers.ipynb``` allows for the execution of SHAMPOO using the background model ```vFrag1``` in ```./BackgroundModels/vFrag1``` with the input file ```ShampooDemo_shampooInput.csv```, and is adjustable to allow for modification of the background disk model file and input file.

# Associated publications
[1] **Oosterloo, M.**, Kamp, I., van Westrenen, W., Dominik, C., **2023**. SHAMPOO: A stochastic model for tracking dust particles under the influence of non-local disk processes. Astronomy and Astrophysics 674, A124 

[2] **Oosterloo, M.**, Kamp, I., van Westrenen, W., **2024**. The effect of nonlocal disk processes on the volatile CHNOS budgets of planetesimal-forming material. Astronomy and Astrophysics 
