from rdkit import Chem
from openbabel import pybel
import glob
import pandas as pd
import os
import multiprocessing as mp
# from wrapt_timeout_decorator import *
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')
rdBase.DisableLog('rdApp.warning')

ob_log_handler = pybel.ob.OBMessageHandler()
ob_log_handler.SetOutputLevel(0)

carbon = Chem.MolFromSmarts("[#6]")



def is_transition_metal(at):
    n = at.GetAtomicNum()
    return (n>=22 and n<=29) or (n>=40 and n<=47) or (n>=72 and n<=79)


def write_output(filenames,outputname="out.sdf"):
    w = Chem.SDWriter(outputname)
    for filename in filenames:
        try:
            if os.stat(filename).st_size / (1024 * 1024)<2.0:
                mol = next(pybel.readfile("cif",str(filename))) 
                molblock=mol.write("mol")    #write out as a molfile string and ship that now into rdkit
                m = Chem.MolFromMolBlock(molblock,removeHs=False,sanitize=True)
                m.SetProp('COD',os.path.basename(filename).split(".")[0])
                w.write(m)
        except Exception:
            pass
    w.close()

def select_molecule(filename):

    blacklist=["/Users/peter/Downloads/cod/cif/2/31/17/2311717.cif","/Users/peter/Downloads/cod/cif/2/10/46/2104629.cif","/Users/peter/Downloads/cod/cif/2/10/59/2105953.cif"]
    if(filename not in blacklist):
        try:
            mol = next(pybel.readfile("cif",str(filename))) 
            molblock=mol.write("mol")    #write out as a molfile string and ship that now into rdkit
            m = Chem.MolFromMolBlock(molblock,removeHs=False,sanitize=True)
            
            if m is not None and len(m.GetSubstructMatches(carbon))>0 and m.GetNumAtoms()> 6  :
                
                if(True not in [is_transition_metal(atom) for atom in m.GetAtoms()]):
                    return filename
        except Exception:
            return None
        return None
    return None


if __name__ == '__main__':
    files=glob.glob('/Volumes/WD3Tb-2/peter/cod/cif/**/*.cif', recursive=True)
    n=0

    pool = mp.Pool(mp.cpu_count())
    results = pool.map(select_molecule, files)
    
    validresults=[el for el in results if el is not None]
    codids=[int(os.path.basename(filename).split(".")[0]) for filename in validresults]
    write_output(validresults)

    # df = pd.read_table('/Users/peter/Downloads/COD_2020jun13.txt',  header=0)

    # dwr=df["COD Number"]

    # common=list(set(dwr) & set(codids))

    # print("common molecules")
    # print(len(common))

    # import numpy as np
    # print("in mine, not in datawarrior")
    # intersect1=np.setdiff1d(codids,dwr)
    # np.savetxt("out_intersect1.csv",intersect1.astype(int),delimiter=",",fmt='%i')
    # print(len(intersect1))

    # print("in datawarrior, not in mine")
    # intersect2=np.setdiff1d(dwr,codids)
    # np.savetxt("out_intersect2.csv",intersect2.astype(int),delimiter=",",fmt='%i')
    # print(len(intersect2))