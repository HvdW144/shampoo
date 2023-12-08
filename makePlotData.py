num = int(input("Which function to execute: "))

print("Loading libraries....")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show, contour, contourf, Circle, Rectangle, pie
from matplotlib import ticker, gridspec
import matplotlib.colors as colors
import prodimopy.read as pread
#import prodimopy.plot as pplot
from scipy.interpolate import griddata
from scipy.integrate import odeint
from scipy.stats import loguniform
from random import random, randint
import pickle
import os
from matplotlib.ticker import FormatStrFormatter
from time import process_time

import seaborn as sns

import shampoo
import pplot
print("Loaded libraries")



pathdict = {}

pathdict["disk"] = "/net/zach/data/users/moosterloo/PhD/Papers/BackgroundModels/ShampooSciencePaper/"
#pathdict["disk"] = "../BackgroundModels/ShampooSciencePaper/"
pathdict["simulation"] = "/net/zach/data/users/moosterloo/PhD/Papers/ShampooSciencePaper/Simulations/"
#pathdict["simulation"] = "./Simulations/"
pathdict["project"] = "/net/zach/data/users/moosterloo/PhD/Papers/ShampooSciencePaper/"
#pathdict["project"] = "./"

disk200 = shampoo.Disk(species=["H2O", "CO", "CO2", "CH4", "CH3OH", "NH3", "H2S", "SO2", "OCS"], folder=pathdict["disk"]+"vFrag1", modelName="ProDiMo.out", 
                         t_index="{:04d}".format(5))

# FUNCTIONS RELATED TO POPULATION STUDY

def countNum(analysis):
    
    R = analysis.disk.model.nx
    Z = analysis.disk.model.nz
    
    countGrid = np.zeros((R,Z))
    
    rArr = (analysis.monomerData["nx"]).values
    zArr = (analysis.monomerData["nz"]).values
    for r in range(R):
        intArrZ = zArr[rArr==r]
        for z in range(Z):
            print(r,z,end="\r")
            countGrid[r,z] = len(intArrZ[intArrZ==z])

    return countGrid

def countNumWeighted(analysis):
    R = analysis.disk.model.nx
    Z = analysis.disk.model.nz
    
    countGrid = np.zeros((R,Z))
    
    rArr = (analysis.monomerData["nx"]).values
    zArr = (analysis.monomerData["nz"]).values
    for r in range(R):
        for z in range(Z):
            print(r,z,end="\r")
            redData = analysis.monomerData[(analysis.monomerData["nx"]==r) & (analysis.monomerData["nz"]==z)] 
            countGrid[r,z] = np.sum(redData["weights"].values)

    return countGrid
    

def countUnique(analysis):
    
    R = analysis.disk.model.nx
    Z = analysis.disk.model.nz
    
    countGrid = np.zeros((R,Z))
    
    M = analysis.monoNum
     
    for r in range(R):
        for z in range(Z):
            print(r,z,end="\r")
            redData = analysis.monomerData[(analysis.monomerData["nx"]==r) & (analysis.monomerData["nz"]==z)]  
            countGrid[r,z] = redData["n"].nunique()
    
    return countGrid

def studyPopulation(analysis, disk200, loadPath=".", mode=None):

        
    if mode=="nUniq":
        tic = process_time()
        #print("Counted all data points in {:0f} CPUs".format(tic-toc))
        print("Counting unique monomer visits")
        grid2 = countUnique(analysis)
        #grid2[grid2<1] = 1
        pickle.dump(grid2, open(loadPath+"/"+analysis.name+"NonLocalCalculations_nMonUniq.analysis", "wb"))
        toc = process_time()
        print("Counted all unique monomers in {:0f} CPUs".format(toc-tic))
    elif mode=="nMon":
        toc = process_time()
        print("Counting data points")
        grid = countNum(analysis)   
        grid[grid<1] = 1
        pickle.dump(grid, open(loadPath+"/"+analysis.name+"NonLocalCalculations_nMonData.analysis", "wb"))
    else:
        tic = process_time()
        #print("Counted all data points in {:0f} CPUs".format(tic-toc))
        print("Counting weighted population positions")
        grid2 = countNumWeighted(analysis)
        #grid2[grid2<1] = 1
        pickle.dump(grid2, open(loadPath+"/"+analysis.name+"NonLocalCalculations_monomerDistributionWeighted.analysis", "wb"))
        toc = process_time()
        print("Calculated monomer distribution in {:0f} CPUs".format(toc-tic))
    
    
    
# CORRUPTION CHECK

def calcTimeFractionCorruption(analysis, disk200, loadPath="."):
       
    loadPath += "Simulations/NonLocalF3"

    fig = figure(figsize=(8,4))
    
    print("Calculating corruption fraction")
    grid = analysis.calcAvgQuant("corruption")
    pickle.dump((grid), open(loadPath+"/"+analysis.name+"NonLocalCalculations_corruption.analysis", "wb"))

# FIND MAX STOKES
def findMaxGrain(analysis, r, z, disk200):
    """
    Finds maximum grain size given radial and vertical position in AU.
    """
    
    rInd = np.argmin(abs(disk200.model.x[:,0]-r))
    zInd = np.argmin(abs(disk200.model.z[rInd,:]-z))
    
    grainData = (analysis.monomerData[(analysis.monomerData["nx"]==rInd)&(analysis.monomerData["nz"]==zInd)])["sa"].values
    if len(grainData)>0:
        saMax = np.max(grainData)
    else:
        saMax = 0 
    
    return saMax
    
    
def maxStokes(analysis, disk200, loadPath="."):
    """
    Calculates the Stokes number for the mean home aggregate size, for different heights above
    the disk midplane.
    """
    
    app = "Simulations/NonLocal1"

    gridList = pickle.load(open(loadPath+app+"/"+analysis.name+"NonLocalCalculations_AVGsAgg.analysis", "rb"))

    heights = [0]
    H = len(heights)
    
    fig = figure(figsize=(12,4))
    frameSA = fig.add_subplot(1,2,1)
    frameST = fig.add_subplot(1,2,2)
    
    model = shampoo.Model(disk=disk200, parameter_folder=loadPath+"/InputFiles/Section3NonLocal1")
    
    for h in range(H):
        rVals = disk200.model.x[:,0]*model.auTOm
        zVals = rVals*heights[h]
        
        S = len(rVals)
        dat_stop = np.zeros(S)
        dat_sagg = np.zeros(S)        
        
        
        for s in range(S):
            print(s,"/",S, end="\r")
            dat_sagg[s] = findMaxGrain(analysis, rVals[s]/model.auTOm, zVals[s]/model.auTOm, disk200)
            model.initGradients(size=dat_sagg[s])
            model.probeEnvironment(rVals[s]/model.auTOm, zVals[s]/model.auTOm, 0, inside_loop=False)
            dat_stop[s] = model.Stokes(0,rVals[s],zVals[s], size=dat_sagg[s], rhoMat=2094*float(model.paraDict["phi"]))
        frameSA.loglog(rVals/model.auTOm, dat_sagg, label="z/r = {:.2f}".format(heights[h]), linewidth=2,zorder=-1)
        frameST.loglog(rVals/model.auTOm, dat_stop, linewidth=2,zorder=-1)
        
        
            
            
  #  for d in range(D):
  #      color = ((D-d)/D,d/D,0)
  #      print(rVals_sagg[d,:])
  #      frame.loglog(rVals_sagg[d,:], rVals_stop[d,:],color=color, linewidth=2)
    
    for frame in [frameSA, frameST]:
        frame.set_xlabel("r (AU)", fontsize=16)
        frame.grid()
        frame.tick_params(labelsize=14)
        frame.set_xlim(0.01,1000)
    
    frameST.set_ylabel("Maximum Stokes number", fontsize=16)
    frameSA.set_ylabel("Maximum aggregate size", fontsize=16)
    frameSA.legend(fontsize=12, ncol=1)
    fig.tight_layout(pad=1)
    
    for fmt in ["svg","pdf","png"]:
        fig.savefig("./Figures/"+"32"+analysis.name+"StokesNumberSupport."+fmt, format=fmt, bbox_inches = 'tight', pad_inches = 0.1) 
        
    
    
# ORIGIN DIAGRAM    
    
def makeOriginHist(analysis):
    
    R = analysis.disk.model.nx
    
    originGrid = np.zeros((R,R))
    
    # select all t=0 points for all monomers.
    r0DataSelect = analysis.monomerData[(analysis.monomerData["t"]==0)] 
    
    for r in range(R):
        print(r+1,"/",R, end="\r")
        # select all monomers which have visited this radius nr.
        rDataSelect = analysis.monomerData[(analysis.monomerData["nx"]==r)]  
        
        # Select monomer Id's
        monomerNums = np.unique(rDataSelect["n"].values).astype(int)
        
        # Iterate over monomers.
        flag = False
        for num in monomerNums:
            r0 = (((r0DataSelect[(r0DataSelect["n"]==num)])["nx"].values).astype(int))
            if len(r0)==1:
                originGrid[r0[0],r] += np.sum((rDataSelect[rDataSelect["n"]==num])["weights"].values)
                flag=True
                
        if flag: #with the flag we make sure we don't get too many zero division errors (we skip empty slices)
            originGrid[:,r] /= np.sum(originGrid[:,r]) #normalize to 1 for each radial slice. 
       
    
    return originGrid

def makeOriginPlot(analysis, disk200, loadPath="."):
           
    loadPath += "Simulations/NonLocal1"   
    
    print("Calculating origin histogram")
    grid = makeOriginHist(analysis)
    pickle.dump(grid, open(loadPath+"/"+analysis.name+"NonLocalCalculations_origins.analysis", "wb"))
    
    print("Done!")
      
        
# AGGREGATE PROPERTIES

def calcAVGagg(analysis, disk200, loadPath=".", quant="all"):
    
    loadPath += "Simulations/NonLocal1"
    

    if quant in ["zm","both","all"]:
        tic = process_time()
        print("Calculating avg zm...")
        grid_zm = analysis.calcAvgQuant("zm")
        pickle.dump(grid_zm, open(loadPath+"/"+analysis.name+"NonLocalCalculations_AVGzMon.analysis", "wb"))
        toc = process_time()
        print("Done in {:.2f} CPUs".format(toc-tic))
    if quant in ["sa", "both","all"]:
        tic = process_time()
        print("Calculating avg sa...")
        grid_zm = analysis.calcAvgQuant("sa")
        pickle.dump(grid_zm, open(loadPath+"/"+analysis.name+"NonLocalCalculations_AVGsAgg.analysis", "wb"))
        toc = process_time()
        print("Done in {:.2f} CPUs".format(toc-tic))
    if quant in ["zm/sa", "all"]:
        tic = process_time()
        print("Calculating avg zm/sa...")
        ratio = analysis.calcAvgQuant("zm/sa")
        pickle.dump(ratio, open(loadPath+"/"+analysis.name+"NonLocalCalculations_AVGaggratio.analysis", "wb"))
        toc = process_time()
        print("Done in {:.2f} CPUs".format(toc-tic))
    
    print("Done!")
    
def checkGrainDistribution(analysis, disk, rPos, zPos, loadPath="."):
    
    loadPath += "Simulations/NonLocal1"
    
    r = np.argmin(abs(disk.model.x[:,0]-rPos))
    z = np.argmin(abs(disk.model.z[r,:]-zPos*rPos))
    
    sa_data = (analysis.monomerData[(analysis.monomerData["nx"]==r) & (analysis.monomerData["nz"]==z)])["sa"]
    zm_data = (analysis.monomerData[(analysis.monomerData["nx"]==r) & (analysis.monomerData["nz"]==z)])["zm"]
    weights = (analysis.monomerData[(analysis.monomerData["nx"]==r) & (analysis.monomerData["nz"]==z)])["weights"]

    fx = figure()
    aux = fx.add_subplot(1,1,1)
    
    data_sa = aux.hist(sa_data, bins=np.logspace(np.log10(5e-8),-3,100), weights=weights)
    data_zm = aux.hist(zm_data, bins=np.logspace(np.log10(5e-8),-3,100), weights=weights)
    
    fig = figure(figsize=(12,4))
    
    yVals_sa = data_sa[0]
    yVals_sa /= np.sum(yVals_sa)
    xVals_sa = ((data_sa[1])[1::]+(data_sa[1])[0:-1])/2
    yVals_zm = data_zm[0]
    yVals_zm /= np.sum(yVals_zm)
    xVals_zm = ((data_zm[1])[1::]+(data_zm[1])[0:-1])/2   
    
    color = next(aux._get_lines.prop_cycler)['color']
    
    frame_sa = fig.add_subplot(1,2,1)
    frame_zm = fig.add_subplot(1,2,2)
    
    frame_sa.loglog(xVals_sa, yVals_sa, linewidth=2)
    frame_zm.loglog(xVals_zm, yVals_zm, linewidth=2)
    
    frame_sa.set_xlabel("Aggregate size $s_a$", fontsize=16)
    frame_zm.set_xlabel("Monomer depth $z_m$", fontsize=16)
    
    # Also fetch the average calculated data to compare.
    grid_zm = pickle.load(open(loadPath+"/"+analysis.name+"NonLocalCalculations_AVGzMon.analysis", "rb"))
    grid_sa = pickle.load(open(loadPath+"/"+analysis.name+"NonLocalCalculations_AVGsAgg.analysis", "rb"))
    
    frame_sa.loglog([grid_sa[r,z]]*2, [1e-99,1], color=(0,0,0), linewidth=3, linestyle="--")
    frame_zm.loglog([grid_zm[r,z]]*2, [1e-99,1], color=(0,0,0), linewidth=3, linestyle="--")
    
    
    for frame in [frame_sa,frame_zm]:
        frame.set_xlim(1e-8, 1e-1)
        frame.set_ylim(1e-4,1e0)
        frame.tick_params(labelsize=12)
       
    for fmt in ["svg","pdf","png"]:
        fig.savefig("./Figures/"+"33"+analysis.name+"DistributionChecking."+fmt, format=fmt, bbox_inches = 'tight', pad_inches = 0.1) 
    
# CHECKING FOR CORRELATION BETWEEN R0 AND AGGREGATE SIZE     
    
def makeGrainSizeGrid(analysis, disk, r):
    # Check whether there is any relation with 
    
    R = disk.model.nx
    
    fx = figure()
    
    # Select the nearest index in the background model.
    rInd = np.argmin(abs(disk.model.x[:,0]-r))
    
    # Select data points at r.
    rDataSelect = analysis.monomerData[(analysis.monomerData["nx"]==rInd)]  
    
    # Look up all unique monomers by first limiting to the initial conditions.
    r0DataSelect = analysis.monomerData[(analysis.monomerData["t"]==0)] 
    # Select monomer Id's
    monomerNums = np.unique(rDataSelect["n"].values).astype(int)
    print("Number of monomers at this radius: ",len(monomerNums))
    
    
    flag1 = True
    M = len(monomerNums)
    m = 1
    
    for num in monomerNums:
        print(m,"/",M, end="\r")
        
        # Select the data for the histogram (note: this is at r=r)
        sa_data = (rDataSelect[(rDataSelect["n"].values==num)])["sa"].values
        weights = (rDataSelect[(rDataSelect["n"].values==num)])["weights"].values
        
        aux = fx.add_subplot(1,1,1)
        data_sa = aux.hist(sa_data, bins=np.logspace(np.log10(5e-8),-1,100), weights=weights)
    
        yVals_sa = data_sa[0]
        
        if flag1:
            xVals_sa = ((data_sa[1])[1::]+(data_sa[1])[0:-1])/2
            X = len(xVals_sa)
            originGrid = np.zeros((X,R))
            flag1 = False

        r0 = (((r0DataSelect[(r0DataSelect["n"]==num)])["nx"].values).astype(int))

        if len(r0)==1:
            originGrid[:,r0[0]] += yVals_sa
            
        m += 1
        
    #Lastly, we normalize every radial slice.
    for n in range(R):
        rSum = np.sum(originGrid[:,n])
        if rSum>0:
            originGrid[:,n] /= rSum
    
    
    
    return originGrid, xVals_sa
    

def calcAMAXR(disk):
    rhos = 2094 # aggregate density in kg/m3
    alphat = 1e-3 # turbulent alpha parameter

    sigma = 2*disk.model.sdg[:,0]*1e4/1e3
    cs = disk.model.soundspeed[:,0]*1e3
    ff = 1/2

    term1 = 2/(3*np.pi)
    term2 = sigma/(rhos*alphat)
    term3 = (1/cs)**2

    amaxr = ff*term1*term2*term3

    return amaxr # output in m
    
def makeGrainSizePlot(analysis, disk200, r=10, loadPath="."):
    """
    Check whether there is any correlation between r0 and the aggregate size at a given point.
    """
           
    loadPath += "Simulations/NonLocal1"   
    

    print("Calculating histograms")
    grid, xVals_sa = makeGrainSizeGrid(analysis, disk200, r)
    pickle.dump((grid, xVals_sa), open(loadPath+"/"+analysis.name+"NonLocalCalculations_grainsizedist"+str(r)+".analysis", "wb"))
    print("Done!")

    
def calcTimeFractionExposed(analysis, disk200, loadPath="."):
       
    loadPath += "Simulations/NonLocal1"

    print("Calculating exposure time fraction")
    grid = analysis.calcAvgQuant("exposed")
    pickle.dump((grid), open(loadPath+"/"+analysis.name+"NonLocalCalculations_exposed.analysis", "wb"))  
    print("Done!")
    
# MAKE NON-LOCAL ICE PLOTS
# LOCAL V.S. NON-LOCAL BAR PLOT

def defineMorvanData():
    # For now we use the maximum values
    
    data = {}
    data["H2O"] = 100
    data["CO"] = 23
    data["CO2"] = 30
    data["CH4"] = 1.5
    data["CH3OH"] = 6.3
    data["NH3"] = .71
    data["H2S"] = 1.5
    data["SO2"] = .23
    data["OCS"] = .40
    
    return data    


def makeBars(analysis, histframe, localData, nonlocalData, n, mode="absolute"):
    """
    Makes the bar diagram for a single panel.
    """
    
    labels = ["Abundance", "", "Species"]
    species = analysis.disk.iceList.copy()
    data = defineMorvanData()
    speciesList = ["H2O","CO","CO2","CH4","CH3OH","NH3","H2S","SO2","OCS"]
    
    if mode=="absolute":
        #print(species)
        species = [r'H$_2$O', 'CO', r"CO$_2$", r'CH$_4$', r'CH$_3$OH', r'NH$_3$', r'H$_2$S', r'SO$_2$', 'OCS']

        forbidden = []

        for popsp in forbidden:
            species.pop(species.index(popsp))

        M = 2
        I = len(species)

        dataVals = np.zeros((M*I, len(labels)), dtype="object")
        
        for m in range(M):
            for i in range(I):

                if species[i]=="R":
                    abundance = 1
                else:
                    if m==0:
                        abundance = localData[i]
                        labelname = "Local"
                    elif m==1:
                        abundance = nonlocalData[i]
                        labelname = "Nonlocal"

                ind = m*I+i# locate index in data array
                dataVals[ind,:] = [abundance, species[i], labelname]

        df = pd.DataFrame(dataVals, columns=labels)        


        barplotbank = sns.barplot(data=df, ax=histframe, x="", y="Abundance", hue="Species", 
                                  saturation=.8)
       
        
        if n==1:
            histframe.legend(bbox_to_anchor=(1.02, 1))
        else:
            histframe.legend([],[],frameon=False)


        histframe.set_yscale("log")
        histframe.set_ylim(1e-8, 1e0)
        #histframe.set_ylim(0, 1.1)
        histframe.grid(which='major', axis='y',zorder=-1)
        histframe.set_ylabel(r"Component mass ($m_m$)", fontsize=16)
        #histframe.set_xlabel(r"$v_{frag}$", fontsize=16)
        histframe.tick_params(labelsize=14)
        
    else:
        
        species = [r'H$_2$O', 'CO', r"CO$_2$", r'CH$_4$', r'CH$_3$OH', r'NH$_3$', r'H$_2$S', r'SO$_2$', 'OCS']

        forbidden = []

        for popsp in forbidden:
            species.pop(species.index(popsp))

        M = 3
        I = len(species)

        dataVals = np.zeros((M*I, len(labels)), dtype="object")

        for m in range(M):
            for i in range(I):

                if m==0:
                    abundance = localData[i]/localData[0]
                    labelname = "Local"
                elif m==1:
                    abundance = nonlocalData[i]/nonlocalData[0]
                    labelname = "Nonlocal"
                elif m==2:
                    abundance = data[speciesList[i]]/100
                    label = "Comets"

                ind = m*I+i# locate index in data array
                dataVals[ind,:] = [abundance, species[i], labelname]

        df = pd.DataFrame(dataVals, columns=labels)        


        barplotbank = sns.barplot(data=df, ax=histframe, x="", y="Abundance", hue="Species", 
                                  saturation=.8)
       
        
        if n==1:
            histframe.legend(bbox_to_anchor=(1.02, 1))
        else:
            histframe.legend([],[],frameon=False)


        histframe.set_yscale("log")
        histframe.set_ylim(1e-4, 1e1)
        #histframe.set_ylim(0, 1.1)
        histframe.grid(which='major', axis='y',zorder=-1)
        histframe.set_ylabel(r"Component mass", fontsize=16)
        #histframe.set_xlabel(r"$v_{frag}$", fontsize=16)
        histframe.tick_params(labelsize=14)


def fetchLocalIceComposition(iceData, n):
    """
    Converts the local monomer ice abundance data to a time-integated average (locations in AU).
    """
    
    I = len(analysis.disk.iceList)
    iceCompositions = [None]*I

    # --- iceData cheatsheet ---
    #iceData = np.zeros((len(tVals), maxNum, P, S))
    # Index 1: # of interpolated timesteps
    # Index 2: # of models per set
    # Index 3: # of plots (the 4 different locations)
    # Index 4: # of ices to compare

    for i in range(I):
        iceCompositions[i] = np.mean(iceData[:,:,n,i])
    
    return iceCompositions

def fetchNonLocalIceComposition(analysis, rPos,zPos):
    """
    Given radial and vertical position (both in AU), determine the local data sets from analysis instance. Also infers ice compositions.
    """
    
    I = len(analysis.disk.iceList)
    iceCompositions = [None]*I
    
    r = np.argmin(abs(analysis.disk.model.x[:,0]-rPos))
    z = np.argmin(abs(analysis.disk.model.z[r,:]-zPos*rPos))
       
    weights = (analysis.monomerData[(analysis.monomerData["nx"]==r) & (analysis.monomerData["nz"]==z)])["weights"]

    mMon = 4/3*np.pi*2094*(5e-8)**3
    
    for i in range(I):
        ice_data = (analysis.monomerData[(analysis.monomerData["nx"]==r) & (analysis.monomerData["nz"]==z)])[analysis.disk.iceList[i]]
        iceCompositions[i] = np.sum(ice_data*weights)/(np.sum(weights))
        iceCompositions[i] /= mMon
        
    return iceCompositions


def makeMasterBarPlot(analysis, mode="absolute", load=False):
    """
    For now, we make panels where we compare average local and average non-local ice composition. Positions are:
       
    rLocs = [3,6,30,60]
    zLocs = [0,0,0,0] in units of z/r!!

    """
    
    # Load the local data
    iceData = pickle.load(open(pathdict["simulation"]+"Section31BGComparisonJuly2023/plot.data", "rb"))

    rLocs = [3,6,30,60]
    zLocs = [0,0,0,0] # in units of z/r!!  
    N = len(zLocs)
    
    fig = figure(figsize=(12,6))
    frameList = [None]*N
    
    for n in range(N):
        if load:
            L_ice, NL_ice = pickle.load(open(pathdict["simulation"]+"Section31BGComparisonJuly2023/"+analysis.name+"histogramData"+str(n)+".analysis", "rb"))   
        else:
            L_ice = fetchLocalIceComposition(iceData, n)
            NL_ice = fetchNonLocalIceComposition(analysis, rLocs[n], zLocs[n]*rLocs[n])
            pickle.dump((L_ice,NL_ice), open(pathdict["simulation"]+"Section31BGComparisonJuly2023/"+analysis.name+"histogramData"+str(n)+".analysis", "wb"))  
            
        frameList[n] = fig.add_subplot(2,2,n+1) 
        makeBars(analysis, frameList[n], L_ice, NL_ice, n, mode)
        frameList[n].set_title("r = {:.2f} AU, z/r = {:.2f}".format(rLocs[n], zLocs[n]), fontsize=16)
    
    fig.tight_layout(pad=1.5)
    for fmt in ["svg","pdf","png"]:
        fig.savefig("./Figures/"+"34"+analysis.name+"NonlocalHistograms"+mode+"."+fmt, format=fmt, bbox_inches = 'tight', pad_inches = 0.1) 
      
# You can make these figures at different points in the disk (maybe again compare local v.s. non-local)


def calculateCompositionHists(analysis, rPos,zPos, quant):
    """
    Given radial position, determine the (radially) local data sets from analysis instance. Also infers ice compositions.
    """
    
    I = len(analysis.disk.iceList)
    iceCompositions = [None]*I
    
    r = np.argmin(abs(analysis.disk.model.x[:,0]-rPos))
    
    if zPos != None:
        z = np.argmin(abs(analysis.disk.model.z[r,:]-zPos*rPos))

        locData = (analysis.monomerData[(analysis.monomerData["nx"]==r) & (analysis.monomerData["nz"]<=z)])
    else:
        locData = (analysis.monomerData[(analysis.monomerData["nx"]==r)])

    mMon = 4/3*np.pi*2094*(5e-8)**3
    
    min_bin = np.log10(5e-8)
    max_bin = -1
    n_bins = 51
    bin_edges = np.logspace(min_bin, max_bin, n_bins)
    
    iceCompositions = np.zeros((n_bins-1, I))
    
    
    for n in range(len(bin_edges) - 1):
        start = bin_edges[n]
        end = bin_edges[n+1]
        red_df = locData[(locData[quant] >= start) & (locData[quant] < end)]
      
        for i in range(I):
            ice_data = red_df[analysis.disk.iceList[i]]
            weights = red_df["weights"]
            iceCompositions[n,i] = np.sum(ice_data*weights)/(np.sum(weights))
            iceCompositions[n,i] /= mMon
        
    binmeds = (bin_edges[1::]+bin_edges[0:-1])/2
        
    return iceCompositions, binmeds


def makeIceAggregateCorrelationPanels(analysis, pathdict):
    
    """
    Is there any correlation between ice composition/abundance and sa/zm?
    """   
        
    rPos = [2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,30,40,50,75,100]
    P = len(rPos)
    
    
    fig = figure(figsize=(12,10))
    frameList = [fig.add_subplot(P, 2, n+1) for n in range(2*P)]
    
    
    zPos = [0.1]*P # in units of z/r

    print("Doing calculations for aggregate correlation")
    for p in range(P):
        print(p+1,"/",P)
        iceCompositionsSA, binmeds = calculateCompositionHists(analysis, rPos[p], zPos[p], "sa")
        iceCompositionsZM, binmeds = calculateCompositionHists(analysis, rPos[p], zPos[p], "zm")    
        datstr = "r"+str(rPos[p])+"z"+str(zPos[p])
        pickle.dump((iceCompositionsSA, binmeds), open(pathdict["project"]+"Figures/"+datstr+analysis.name+"IceAggregateCorrelationSA.data", "wb"))
        pickle.dump((iceCompositionsZM, binmeds), open(pathdict["project"]+"Figures/"+datstr+analysis.name+"IceAggregateCorrelationZM.data", "wb"))
         
        I = len(analysis.disk.iceList)
        
        for i in range(I):
            frameList[2*p].loglog(binmeds, iceCompositionsSA[:,i], label=analysis.disk.iceList[i], linewidth=2)
            frameList[2*p+1].loglog(binmeds, iceCompositionsZM[:,i], label=analysis.disk.iceList[i], linewidth=2)
                    
        frameList[2*p].set_xlabel(r"$s_a$", fontsize=16)
        frameList[2*p+1].set_xlabel(r"$z_m$", fontsize=16)
                
    for frame in frameList:
        frame.set_ylabel(r"$m_{ice}/m_m$", fontsize=16)
        frame.set_ylim(1e-6,1e2)
        frame.set_xlim(5e-8,1e-3)
        frame.grid()  
        frame.tick_params(labelsize=14)
        
    frameList[1].legend(fontsize=12, ncol=3, loc="upper right")
        
    fig.tight_layout()
    for fmt in ["svg","pdf","png"]:
        fig.savefig("./Figures/"+"34"+analysis.name+"AggregateIceDistributions."+fmt, format=fmt, bbox_inches = 'tight', pad_inches = 0.1) 
     
    
    
# MAKE RADIAL ICE PROFILES
    
def fetchIces(analysis, pathdict, species):
    """
    Calculates mice/mmonomer as a function of r in the midplane for given species.
    """
    
    
    uTOkg = 1.660539067e-27 

    masses = {}
    dataName = np.loadtxt(pathdict["disk"]+"vFrag1/AdsorptionEnergies.in", dtype="str", comments="#", usecols=(0), encoding=None)
    dataNum = np.loadtxt(pathdict["disk"]+"vFrag1/AdsorptionEnergies.in", dtype="float", comments="#", usecols=(1,2), encoding=None)
        
    N = len(dataName)
    for n in range(N):
        masses["m"+dataName[n]] = dataNum[n,1]*uTOkg # in kg!   

    R = len(analysis.disk.model.x[:,0])
    abunProDiMo = np.zeros(R)
    abunSHAMPOO = np.zeros(R)         
    mMon = 4/3*np.pi*2094*(5e-8)**3
    for r in range(R):
        print(r+1,"/",R, end="\r")
        locData = (analysis.monomerData[(analysis.monomerData["nx"]==r) & (analysis.monomerData["nz"]==0)])
        iceMass = locData[species].values
        weights = locData["weights"].values
        abunSHAMPOO[r] = np.sum(iceMass*weights)/(np.sum(weights))/mMon
        abunProDiMo[r] = analysis.disk.expectedIce(analysis.disk.model.x[r,0], 0, species=species)
        
    return abunSHAMPOO, abunProDiMo

def showRadialIceBehaviour(analysis, pathdict):
        
    abuns = {}
    abunp = {}
    
    for species in ["H2O", "CO", "CO2", "CH4", "CH3OH", "NH3", "H2S", "SO2", "OCS"]:
        print("Fetching radial profile of: ", species)
        abuns[species], abunp[species] = fetchIces(analysis, pathdict, species)
        
           
    pickle.dump((abuns, abunp), open("./Figures"+"/RadialAbundances"+analysis.name+".data", "wb"))
  
# MAKE RATIO FIGS

def findweights(sp, A, B):
    wA = 0
    wB = 0
    # CAREFUL, IF WE GO TO H, WE NEED TO ACCOUNT FOR DOUBLE APPEARENCES AS CH3OH
    
    
    if A in sp:
        ind = sp.index(A)+1
        if ind<len(sp):
            if sp[ind].isnumeric():
                wA = int(sp[ind])
            else:
                wA = 1
        else:
            wA = 1
    if B in sp:
        ind = sp.index(B)+1
        if ind<len(sp):
            if sp[ind].isnumeric():
                wB = int(sp[ind])
            else:
                wB = 1
        else:
            wB = 1
    
    return wA, wB

def estimateRatio(analysis, A, B, pathdict):
    """
    Estimates the NUMBER ratio of element A and B.
    """
    
    I = len(analysis.disk.iceList)
    preA = [None]*I
    preB = [None]*I
    
    
    KB = 1.38064852e-23
    uTOkg = 1.660539067e-27 

    masses = {}
    dataName = np.loadtxt(pathdict["disk"]+"vFrag1/AdsorptionEnergies.in", dtype="str", comments="#", usecols=(0), encoding=None)
    dataNum = np.loadtxt(pathdict["disk"]+"vFrag1/AdsorptionEnergies.in", dtype="float", comments="#", usecols=(1,2), encoding=None)
        
    N = len(dataName)
    for n in range(N):
        masses["m"+dataName[n]] = dataNum[n,1] # in amu!   
    
    chnosM = [1.00797, 12.011, 14.0067, 15.9994, 32.06]
    chnos = ["H", "C", "N", "O", "S"]
    C = len(chnos)
    for c in range(C):
        masses["m"+chnos[c]] = chnosM[c]
    
    
    # First we fetch the appropriate prefactors
    for i in range(I):
        preA[i], preB[i] = findweights(analysis.disk.iceList[i], A, B)

    # Subsequently the average mass at each r for each species. convert to number.
    R = len(analysis.disk.model.x[:,0])
    massA = np.zeros(R)
    numA = np.zeros(R)
    massB = np.zeros(R)
    numB = np.zeros(R)   
    massApro = np.zeros(R)
    massBpro = np.zeros(R)
    numApro = np.zeros(R)
    numBpro = np.zeros(R)
    
    for r in range(R):
        print(r+1,"/",R, end="\r")
        locData = (analysis.monomerData[(analysis.monomerData["nx"]==r) & (analysis.monomerData["nz"]==0)])
        
        for i in range(I):
            # calculate mass of species from shampoo
            iceMass = locData[analysis.disk.iceList[i]]
            weights = locData["weights"]
            avgSPice = np.sum(iceMass*weights)/(np.sum(weights))
            
            # calculate mass of element A
            massA[r] += avgSPice*(preA[i]*masses["m"+A]/masses["m"+analysis.disk.iceList[i]])
            numA[r] += massA[r]/masses["m"+A]
            
            # "" element B
            massB[r] += avgSPice*(preB[i]*masses["m"+B]/masses["m"+analysis.disk.iceList[i]])
            numB[r] += massB[r]/masses["m"+B]
            
            # Calculate mass of species from prodimo
            prodimoNum = masses["m"+analysis.disk.iceList[i]]*10**((analysis.disk.data["iceAbun"+analysis.disk.iceList[i]])[r,0])
            massApro[r] += prodimoNum*(preA[i]*masses["m"+A]/masses["m"+analysis.disk.iceList[i]])
            numApro[r] += massApro[r]/masses["m"+A]
            
            # "" element B
            massBpro[r] += prodimoNum*(preB[i]*masses["m"+B]/masses["m"+analysis.disk.iceList[i]])
            numBpro[r] += massBpro[r]/masses["m"+B]
    
    
    massrat = {}
    numrat = {}
    
    massrat["shampoo"] = massA/massB
    numrat["shampoo"] = numA/numB
    massrat["prodimo"] = massApro/massBpro
    numrat["prodimo"] = numApro/numBpro
    
    return massrat, numrat
    


def makeRatioPlot(analysis, pathdict):
    
    R = len(analysis.disk.model.x[:,0])
    
    S = 3
    A = ["C", "N","S"]
    B = ["O"]*S
    
    for s in range(S):
        print(s+1,"/",S)
        massrat, numrat = estimateRatio(analysis, A[s], B[s], pathdict)
        pickle.dump((massrat, numrat), open("./Figures/"+A[s]+B[s]+"RatioFigures"+analysis.name+".data", "wb"))
                

def probeIceDistributions(analysis, pathdict, rPos, zPos=0):
    
    iceHists = {}
    iceBins = {}
    
    r = np.argmin(abs(rPos-analysis.disk.model.x[:,0]))
    
    locData = (analysis.monomerData[(analysis.monomerData["nx"]==r) & (analysis.monomerData["nz"]==0)])
    
    N = len(analysis.disk.iceList)
    binNo = 100
    
    mMon = 4/3*np.pi*2094*(5e-8)**3
    
    for n in range(N):
  
        iceMass = locData[analysis.disk.iceList[n]].values
    
        if np.any(np.isnan(iceMass)):
            print("nans found for ", analysis.disk.iceList[n])
    
        #if n==0:
            #print("Data points: ",len(iceMass))
            #print(iceMass)
        if len(iceMass)>0:
            minIce = max([np.min(iceMass)*0.9,1e-20*mMon])
            maxIce = max([np.max(iceMass)*1.1,minIce*1.1])
            
            iceMass[iceMass<minIce] = minIce*(1+1e-20)
            
            if np.any(iceMass>1e-20*mMon):
                binEdges = np.logspace(np.log10(minIce), np.log10(maxIce), binNo+1) # note that this is in kg

                hist, bins = np.histogram(iceMass, bins=binEdges)
                #if n==0:
                    #print(hist,bins)

                iceHists[analysis.disk.iceList[n]] = hist
                iceBins[analysis.disk.iceList[n]] = (bins[1::]-bins[0:-1])/(2*mMon)
            else:
                iceHists[analysis.disk.iceList[n]] = np.array([0])
                iceBins[analysis.disk.iceList[n]] = np.array([0])
        else:
            iceHists[analysis.disk.iceList[n]] = np.array([0])
            iceBins[analysis.disk.iceList[n]] = np.array([0])
    print(iceHists,iceBins)    
    pickle.dump((iceHists, iceBins), open(pathdict["project"]+"Figures/IceProbing"+analysis.name+"r"+str(r)+"z"+str(0)+".data", "wb"))
    
 

print("Loading data...")
analysis_name = "NonLocalF3_cleaned3"

print("File:",pathdict["simulation"]+"NonLocalF3"+"/"+analysis_name+".analysis")
analysis = pickle.load(open(pathdict["simulation"]+"NonLocalF3"+"/"+analysis_name+".analysis", "rb"))
analysis.name = analysis_name
print("Data loaded ({:.0f} monomers)".format(analysis.monoNum))
print("Number of data points: ", len((analysis.monomerData["n"]).values))


if num==0:
    studyPopulation(analysis, disk200, loadPath=pathdict["project"], mode="nUniq") #
elif num==15:
    studyPopulation(analysis, disk200, loadPath=pathdict["project"], mode="nMon") #
elif num==16:
    studyPopulation(analysis, disk200, loadPath=pathdict["project"], mode="weighted") #
elif num==1:
    calcTimeFractionCorruption(analysis, disk200, loadPath=pathdict["project"]) #
elif num==2:
    maxStokes(analysis, disk200, loadPath=pathdict["project"]) # 
elif num==3:
    makeOriginPlot(analysis, disk200, loadPath=pathdict["project"]) # 
elif num==4:
    calcAVGagg(analysis, disk200, loadPath=pathdict["project"], quant="zm") # 
elif num==5:
    calcAVGagg(analysis, disk200, loadPath=pathdict["project"], quant="sa") # 
elif num==6:
    calcAVGagg(analysis, disk200, loadPath=pathdict["project"], quant="zm/sa") # 
elif num==7:
    checkGrainDistribution(analysis, disk200, 10, 0.15, loadPath=pathdict["project"]) 
elif num==8:
    calcTimeFractionExposed(analysis, disk200, loadPath=pathdict["project"]) # 
elif num==9:
    makeGrainSizePlot(analysis, disk200, r=20, loadPath=pathdict["project"]) # 
elif num==10:
    makeMasterBarPlot(analysis, mode="absolute") # 
elif num==11:
    makeIceAggregateCorrelationPanels(analysis, pathdict) # 
elif num==12:
    showRadialIceBehaviour(analysis, pathdict) # 
elif num==13:
    makeRatioPlot(analysis, pathdict) # 
elif num==14:
    for r in range(150):
        print(r+1,"/",150)
        probeIceDistributions(analysis, pathdict, rPos=disk200.model.x[r,0])

print("Finished executing plotting scripts")
