# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 14:03:27 2023

@author: Bangchen Yin

A script which automatically generate input folder for multiple examples
"""
import os
import argparse
import ase
from ase.io import read
def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument("PARAM", type=str,default="param_enhanced_sampling.json", help="The parameters of the generator")
    parser.add_argument(
        "MACHINE", type=str,default="machine.json", help="The settings of the machine running the generator"
    )
    args = parser.parse_args()

 
    import ruamel
    from monty.serialization import dumpfn, loadfn

  
    jobs=3
    n_batch=10000

    current_directory = os.getcwd()
    subdirectories = [name for name in os.listdir(current_directory) if os.path.isdir(os.path.join(current_directory, name))]
    
    #print(subdirectories)
    for sub in subdirectories:
        if os.path.exists(sub+"/POSCAR"):
            jdata = loadfn(args.PARAM)
            mdata = loadfn(args.MACHINE)
    
            atoms=read(sub+"/POSCAR")
            f=open(sub+"/POSCAR","r")
            lines=f.readlines()
            name=''.join(lines[5].strip().split())
            os.system("mkdir "+name)
            os.system("cp "+sub+"/POTCAR "+name+"/")
            os.system("mkdir "+name+"/Data")
            os.system("mkdir "+name+"/Data/devi")
            for i in range(3):
                os.system("mkdir "+name+"/Data/devi/0"+str(i))
                os.system("cp  "+sub+"/POSCAR "+name+"/Data/devi/0"+str(i))
            os.system("cp -r "+sub+"/deepmd "+name+"/Data/")
            os.system("cp INCAR "+sub)
            os.system("cp dimer.py "+name)
            os.system("cp run.py "+name)
            os.system("cp arginfo.py "+name)
            os.system("cp INCAR "+name)
            os.system("cp input.lammps "+name)
            jdata['fp_pp_files']=['POTCAR']
            jdata["fp_incar"]="INCAR"
            jdata['type_map']=lines[5].strip().split()
            jdata['mass_map']='auto'
            jdata['default_training_param']['model']['type_map']=lines[5].strip().split()
            print(lines[5].strip())
            user_input = input("sel_number for each element, Input separated by commas: ")
            num_list=[]
            for fn in  user_input.split(','):
                num_list.append(int(fn))
            jdata['default_training_param']['model']['descriptor']["sel"]=num_list 
            jdata['default_training_param']["training"]["stop_batch"]=n_batch
           # absolute_path = os.path.abspath(__file__)
           #TODO: These can be modified based on your situation
            command="CUDA_VISIBLE_DEVICES=1 srun python "+current_directory+"/"+name+"/dimer.py "
           
            para=" True  "+str(len(atoms)-1)+" 0.05 "+  str(len(atoms)-1)+"  2"
            mdata["model_devi"][0]["command"]=command+para
            dumpfn(jdata,name+"/par.json")
            dumpfn(mdata,name+"/machine.json")                
        else:
            print("POSCAR not found!")



if __name__ == "__main__":
    _main()

