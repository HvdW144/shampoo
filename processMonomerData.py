import shampoo 
import pickle

def doAnalysis():
    diskpath = "/net/zach/data/users/moosterloo/PhD/Papers/BackgroundModels/ShampooSciencePaper"
    homepath = "/net/zach/data/users/moosterloo/PhD/Papers/ShampooSciencePaper"
    
    disk200 = shampoo.Disk(species=["H2O", "CO", "CO2", "CH4", "CH3OH", "NH3", "H2S", "SO2", "OCS"], folder=diskpath+"/vFrag1", modelName="ProDiMo.out", t_index="{:04d}".format(5))
    
    loadPath= homepath+"/Simulations/NonLocalF3/"
    analysis = shampoo.Analysis(disk=disk200)
    analysis.loadModels(read=True, monoNum=64000, cleaningIndex=3, removeCorrupted="selective", loadPath=loadPath)
    ##### SHOULD SELECTIVELY REMOVE 
    print("Writing file...")
    pickle.dump(analysis, open(loadPath+"NonLocalF3_cleaned3.analysis", "wb"))
    print("Finished file preparation!")
    
doAnalysis()
