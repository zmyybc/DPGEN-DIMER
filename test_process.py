# -*- coding: utf-8 -*-
"""
This code could generate muliple search trajectory for the given initial configuration, and give the predicted energy.

Created on Thu Jul  6 14:03:27 2023

@author: Bangchen Yin
"""

import os
from ase.dimer import DimerControl, MinModeAtoms, MinModeTranslate
from ase.io import read
from ase.io.trajectory import Trajectory
from deepmd.calculator import DP
import numpy as np
import argparse
def write_lamptrj(atoms,forces,step,name):
    f_conf=open("conf.lmp","r")
    lines=f_conf.readlines()
    for line in lines:
        if "atom types" in line:
            n_atoms=int(line.split()[0])
    fname="./"+name+"/"+str(step)+".lammpstrj"
   
    posname="./"+name+"/"+str(step)+".poscar"
    atoms.write(posname)
    f=open(fname,"w")
    f.write("ITEM: TIMESTEP\n")
    f.write(str(step)+'\n')
    f.write("ITEM: NUMBER OF ATOMS\n")
    f.write(str(len(atoms))+'\n')
    f.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
    o_line="0.00 "+str(atoms.cell[0][0])+" 0.00\n"
    t_line="0.00 "+str(atoms.cell[1][1])+" 0.00\n"
    th_line="0.00 "+str(atoms.cell[2][2])+" 0.00\n"
    f.write(o_line)
    f.write(t_line)
    f.write(th_line)
    print(atoms.cell)
    f.write("ITEM: ATOMS id type x y z fx fy fz\n")
    count=1
    pre=""
    for atom in atoms:
        xyz=atom.position
        f.write(str(atom.index+1)+" ")
        if(atom.symbol==pre)or(pre==""):
          f.write(str(count)+" ")
        else:
            count+=1
            f.write(str(count)+" ")
        pre=atom.symbol
        for ii in xyz:
            f.write(str(ii)+" ")
        for idf in forces[atom.index]:
            f.write(str(idf)+" ")
        f.write("\n")
def main(i):
 
    #TODO: modify the path of initial configuration and path of the trained model
    atoms=read(args.ini_conf_path)
    atoms.calc = DP(model="args.model_path")
    #TODO: modify the dimer parameters
    with DimerControl(initial_eigenmode_method='displacement',maximum_translation=0.1,
                      displacement_method='gauss', displacement_center=20,dimer_separation=0.001,displacement_radius=3,logfile=None, #TODO: make it adjustable
                      ) as d_control:
        d_atoms = MinModeAtoms(atoms, d_control)
       
        d_atoms.displace()

        with MinModeTranslate(d_atoms, trajectory='dimer_method.traj',
                              logfile='log') as dim_rlx:
            
            
            dim_rlx.run(fmax=0.02,steps=200)
            d_atoms.write('POSCAR')
            IS=read("POSCAR")
            FS=read("POSCAR")
            IS.calc=DP(model="args.model_path")
            FS.calc=DP(model="args.model_path")
            IS.positions+=d_atoms.eigenmodes[0]
            FS.positions-=d_atoms.eigenmodes[0]
            dyn1=BFGS(IS,trajectory="IS.traj")
            dyn2=BFGS(FS,trajectory="FS.traj")
            dyn1.run(fmax=0.02,steps=200)
            dyn2.run(fmax=0.02,steps=200)
            IS.write("IS.vasp")
            FS.write("FS.vasp")
            traj1 = Trajectory('IS.traj')
            traj2=Trajectory('FS.traj')
            traj=[]
            for atoms in traj1:
                traj.append(atoms)
            for atoms in traj2:
                traj.append(atoms)
            if not os.path.exists(name):
              os.system("mkdir "+name)
            count=0
            for atoms in traj:
                forces= np.zeros([len(atoms),3])
                write_lamptrj(atoms,forces,count,name)
                count+=1
            coor=np.array([coor])
            box=np.array([box])
            atoms=traj[-1]
            atoms.calc=DP(model="args.model_path")
            return            atoms.get_potential_energy()


if __name__ == "__main__":
  f=open("ener.txt","w")
  parser = argparse.ArgumentParser()
  parser.add_argument("ini_conf_path", type=str,default=".", help="path of initial configuration, better be a poscar")
  parser.add_argument("model_path", type=str,default=".", help="path of the trained model")
  args = parser.parse_args()
  for i in range(10):

     print("SEARCH "+str(i+1),main(i,args),file=f)

