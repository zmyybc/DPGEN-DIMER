# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 14:03:27 2023

@author: Bangchen Yin
"""
import os
from ase.build import bulk
from ase import Atom, Atoms
from ase.dimer import DimerControl, MinModeAtoms, MinModeTranslate
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.io import read,write
from ase.io.trajectory import Trajectory
from deepmd.infer import calc_model_devi
from deepmd.calculator import DP
import numpy as np
import faulthandler
from ase.calculators.vasp.vasp import Vasp
def write_lamptrj(atoms,forces,step,name):
    f_conf=open("conf.lmp","r")
    lines=f_conf.readlines()
    for line in lines:
        if "atom types" in line:
            n_atoms=int(line.split()[0])
    fname="./"+name+"/"+str(step)+".lammpstrj"
   
    posname="./"+name+"/"+str(step)+".poscar"
    #f1=open(posname,"w")
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
def parse_parameters():
    f=open("input.lammps","r")
    lines=f.readlines()
    parameters={}
    for line in lines:
       if not line.startswith("#"):
        start=line.find(" ")
        if "mass" in line:
            if "masses" in parameters.keys():
                parameters["masses"].append(line[start+1:].strip())
            else:
                parameters["masses"]=[line[start+1:].strip()]
        elif "pair_coeff" in line:
             parameters["pair_coeff"]=[line[start+1:].strip()]
        else:
             parameters[line[:start]]=line[start+1:].strip()
    return parameters
def read_conf():
   
    atoms=read("conf.lmp",format="lammps-data",style="atomic")
    
    return atoms
def main(i):
    name="traj"+str(i)
    parameters=parse_parameters()
    #TODO: modify the path of initial configuration
    atoms=read("/home/yinbc/test/dpgen/dpgen/graph/iter.000009/01.model_devi/confs/000.0000.poscar")
    atoms.calc = DP(model="../graph.001.pb")
    with DimerControl(initial_eigenmode_method='displacement',maximum_translation=0.1,
                      displacement_method='gauss', displacement_center=20,dimer_separation=0.001,displacement_radius=3,logfile=None, #TODO: make it adjustable
                      ) as d_control:
        d_atoms = MinModeAtoms(atoms, d_control)
       
        d_atoms.displace()

        with MinModeTranslate(d_atoms, trajectory='dimer_method.traj',
                              logfile='log') as dim_rlx:
            
            
            dim_rlx.run(fmax=0.02,steps=200)
            traj = Trajectory('dimer_method.traj')
            if not os.path.exists(name):
              os.system("mkdir "+name)
            count=0
            box=[]
            coor=[]
            atype=[]
            k=0
            pre=""
            for atom in atoms:
               if atom.symbol != pre:
                 pre=atom.symbol
                 k+=1
               atype.append(k-1)
            for atoms in traj:
              
                if box ==[]:
                    box= np.array(atoms.cell).reshape([1,-1])[0]
                else:
                   box=np.append(box,np.array(atoms.cell).reshape([1,-1])[0])
                 
                if coor==[]:
                    coor=np.array(atoms.positions).reshape([1,-1])[0]
                else:
                    coor=np.append(coor,np.array(atoms.positions).reshape([1,-1])[0])
              
                forces= np.zeros([21,3])
                write_lamptrj(atoms,forces,count,name)
                 
                 
                count+=1
            coor=np.array([coor])
            box=np.array([box])
            atoms=traj[-1]
            atoms.calc=DP(model="../graph.001.pb")
            return            atoms.get_potential_energy()


if __name__ == "__main__":
  f=open("ener.txt","w")
  for i in range(10):

     print("SEARCH "+str(i+1),main(i),file=f)

