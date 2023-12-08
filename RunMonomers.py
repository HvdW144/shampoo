import shampoo
import pickle
import numpy as np
from time import process_time

#N = int(input("Give number of monomers to run: "))

def runH2O():
    """
    Testing the collisional equilibrium function. Do not use in paper.
    """
    N=150

    disk = shampoo.Disk(species=["H2O", "CO", "CO2", "CH4", "NH3", "H2S"], folder="../ShampooBackground", modelName="ProDiMo.out", 
                         t_index="{:03d}".format(5))
    savePath="./Simulations/WaterIceColEqTest/"


    for n in range(N):
        print("Running model ",n+1, "/", N)
        mod = shampoo.Model(disk=disk, parameter_folder="./InputFiles/H2OIceLineColEq", breakIce=False)
        mod.integrateMonomer()
        name = str(mod.seedStart)
        del mod.disk

        pickle.dump(mod, open(savePath+name+".pickle", "wb"))
        print("Saved model with seed: ", mod.seedStart)

def runLocal():
    """
    Systematic run for different grid points. Do not use in paper.
    """
    N=19

    disk = shampoo.Disk(species=["H2O", "CO", "CO2", "CH4", "NH3", "H2S", "SO2", "S"], folder="../ShampooBackground", modelName="ProDiMo.out", 
                         t_index="{:04d}".format(5))
    savePath="./Simulations/ModelComparison/"

    R = len(disk.rVals)
    Z = len(disk.zVals)

    
    toc = process_time()
    
    m=1
    
    for n in range(N):
        for r in range(R):
            for z in range(Z):
                
                r0 = disk.rVals[r]
                z0 = disk.zVals[z]*disk.rVals[r]
        
                if (r0<1)or(z0/r0>0.3)or(r%10>0)or(z%10>0):
                    print("-"*100)
                    print("Running monomer set ",n+1, "/", N, "-- set is {:.2f} % complete".format(m*100/(R*Z)))
                    print("Skipping r point ",r+1, "/", R, " -- r = {:.2f} AU".format(r0))
                    print("Skipping z point ",z+1, "/", Z, " -- z/r = {:.2f}".format(z0/r0))
                    print("-"*100)
                else:
                    print("-"*100)
                    print("Running monomer set ",n+1, "/", N, "-- set is {:.2f} % complete".format(m*100/(R*Z)))
                    print("Running r point ",r+1, "/", R, " -- r = {:.2f} AU".format(r0))
                    print("Running z point ",z+1, "/", Z, " -- z/r = {:.2f}".format(z0/r0))
                    print("-"*100)

                    mod = shampoo.Model(disk=disk, parameter_folder="./InputFiles/PopulationRunLocal", breakIce=False)

                    mod.integrateMonomer(r0=r0, z0=z0)

                    name = str(mod.seedStart)+"nr{:.0f}".format(r)+"nz{:.0f}".format(z)
                    del mod.disk
                    del mod.para
                    pickle.dump(mod, open(savePath+name+".pickle", "wb"))
                    print("Saved model with seed: ", mod.seedStart)
                m +=1

    tic = process_time()
    print("Finished ",N," monomer set(s) in {:.1f} CPUs".format(abs(tic-toc)))
       

def evaluateAtIceLines():
    """
    Makes a systematic study of the ice just behind the icelines. Studies for 4 populations:
    
    """
    
    abspath = "/net/zach/data/users/moosterloo/PhD/Papers/"
    #abspath = "../"
    
    diskfolder = abspath+"BackgroundModels/ShampooSciencePaper/vFrag1"
    modelfolder = abspath+"ShampooSciencePaper/Simulations/Section33BGComparison"
    
    disk = shampoo.Disk(species=["H2O", "CO", "CO2", "CH4", "CH3OH", "NH3", "H2S", "SO2", "OCS"], folder=diskfolder, modelName="ProDiMo.out", 
                         t_index="{:04d}".format(5))
    
    
    rLocs = [3,6,30,60]
    zLocs = [0,0,0,0] # in units of z/r!!
    C = len(rLocs)
    
    nx = [np.argmin(abs(r-disk.model.x[:,0])) for r in rLocs]
    nz = [np.argmin(abs(abs(z)-disk.model.z[-1,:]/disk.model.x[-1,:])) for z in zLocs]
    
    r0List = [disk.model.x[nx[c],0] for c in range(C)]
    z0List = [disk.model.z[-1,nz[c]]/disk.model.x[-1,nz[c]]*r0List[c] for c in range(C)]
    
    #print(r0List,z0List)
    
    N=100
    # MAKE SURE COLLISIONS ARE TRUE, MAY BE FALSE DUE TO MODEL TEST
    for n in range(N):
        for c in range(C):
            mod = shampoo.Model(disk=disk, parameter_folder="./InputFiles/Section33BGComparison")
            mod.integrateMonomer(r0=r0List[c], z0=z0List[c])
            name = str(mod.seedStart)
            del mod.disk

            pickle.dump(mod, open(modelfolder+"/"+str(rLocs[c])+"/"+name+".pickle", "wb"))
            print("Saved model with seed: ", mod.seedStart)

def runNonLocal():
    """
    Code to just run as many non local monomers as possible.
    
    Settings:
    - Initial positions between r in [0.5-100] AU, z/r in [-0.1,0.1]  
    """
    
    abspath = "/net/zach/data/users/moosterloo/PhD/Papers/"
    #abspath = "../"
    
    diskfolder = abspath+"BackgroundModels/ShampooSciencePaper/vFrag1"
    modelfolder = abspath+"ShampooSciencePaper/Simulations/NonLocalF3"
    
    disk = shampoo.Disk(species=["H2O", "CO", "CO2", "CH4", "CH3OH", "NH3", "H2S", "SO2", "OCS"], folder=diskfolder, modelName="ProDiMo.out", 
                         t_index="{:04d}".format(5))
    
    N=800
    
    for n in range(N):
        try:
            mod = shampoo.Model(disk=disk, parameter_folder="./InputFiles/Section3NonLocalF3")
            mod.integrateMonomer()
            name = str(mod.seedStart)
            del mod.disk

            pickle.dump(mod, open(modelfolder+"/"+name+".pickle", "wb"))
            print("Saved model with seed: ", mod.seedStart)
        except:
            print("A monomer broke SHAMPOO...")
            pass

import contextlib
import sys

with contextlib.redirect_stdout(None), contextlib.redirect_stderr(None):
    evaluateAtIceLines()
    
    
