# -*- coding: utf-8 -*-
"""
This code could generate muliple dimer search to get saddle points and corresponding barriers for the given initial configuration.

Created on Thu Jul  6 14:03:27 2023

@author: Bangchen Yin
"""

import os
from ase.dimer import DimerControl, MinModeAtoms, MinModeTranslate
from ase.io import read
from deepmd.calculator import DP
import argparse

def main(i):
 
    #TODO: modify the path of initial configuration and path of the trained model
    atoms=read(args.ini_conf_path)
    atoms.calc = DP(model="args.model_path")
    E0=atoms.get_potential_energy()
    #TODO: modify the dimer parameters
    with DimerControl(initial_eigenmode_method='displacement',maximum_translation=0.1,
                      displacement_method='gauss', displacement_center=20,dimer_separation=0.001,displacement_radius=3,logfile=None, #TODO: make it adjustable
                      ) as d_control:
        d_atoms = MinModeAtoms(atoms, d_control)
       
        d_atoms.displace()

        with MinModeTranslate(d_atoms, trajectory='dimer_method.traj',
                              logfile='log') as dim_rlx:
            dim_rlx.run(fmax=0.02,steps=200)
            f=open("barrier_sad.txt","a")
            os.system("mkdir saddle")
            d_atoms.write("saddle/"+str(i)+'.vasp')
            atoms=read("saddle/"+str(i)+'.vasp')
            atoms.calc=DP(model="args.model_path")
            f.write(str(atoms.get_potential_energy()-E0))
            f.write("  ")
            f.write("\n")
            return            atoms.get_potential_energy()-E0


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("ini_conf_path", type=str,default=".", help="path of initial configuration, better be a poscar")
  parser.add_argument("model_path", type=str,default=".", help="path of the trained model")
  args = parser.parse_args()
  for i in range(10):
     main(i,args)

