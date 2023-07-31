# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 14:03:27 2023

@author: ASUS
"""
import os
import argparse
from ase.build import bulk
from ase import Atom, Atoms
from ase.dimer import DimerControl, MinModeAtoms, MinModeTranslate
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.io import read,write
from ase.io.trajectory import Trajectory
from deepmd.infer import calc_model_devi
#from ase.calculators.lammpsrun import LAMMPS
from deepmd.calculator import DP
import numpy as np
import faulthandler
def write_lamptrj(atoms,forces,step):
    f_conf=open("conf.lmp","r")
    lines=f_conf.readlines()
    for line in lines:
        if "atom types" in line:
            n_atoms=int(line.split()[0])
    fname="./traj/"+str(step)+".lammpstrj"
   
    posname="./traj/"+str(step)+".poscar"
    #f1=open(posname,"w")
   # atoms.write(posname)
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
    #f=open("conf.lmp","r")
    atoms=read("conf.lmp",format="lammps-data",style="atomic")
    
    return atoms
def main(args):
    parameters=parse_parameters()
    #parameters["pair_style"]="morse 4.5"
    #atoms=read_conf()
   # atoms=read("water.cif")
#    parameters = {'pair_style': 'eam/alloy',
 #             'pair_coeff': ['* * NiAlH_jea.eam.alloy H Ni']}

  #  files = ['NiAlH_jea.eam.alloy']

   # Ni = bulk('Ni', cubic=True)
    #H = Atom('H', position=Ni.cell.diagonal()/2)
    atoms=read("../confs/"+"000.0000.poscar")
   # atoms = read_conf()
    if args.fix:

      mask=[atom.index<args.index for atom in atoms]
      print(mask)
    #atoms.set_constraint(FixAtoms(mask=mask))
 #   mask=[atom.tag<(len(atoms)-1) for atom in atoms]
      atoms.set_constraint(FixAtoms(mask=mask))
   # atoms.write("POSCAR")
   # lammps = LAMMPS(parameters=parameters)
    atoms.calc = DP(model="../graph.001.pb")
    with DimerControl(initial_eigenmode_method='displacement',maximum_translation=args.trans,
                      displacement_method='gauss', displacement_center=args.center,dimer_separation=0.0005,displacement_radius=args.radius,logfile=None, #TODO: make it adjustable
                      ) as d_control:
        d_atoms = MinModeAtoms(atoms, d_control)
        # Displace the atoms
        d_atoms.displace()
        # Converge to a saddle point
     #   f1=open("model_devi_tot.out","w")
       # f1.write("#       step         max_devi_v         min_devi_v         avg_devi_v         max_devi_f         min_devi_f         avg_devi_f\n")
        with MinModeTranslate(d_atoms, trajectory='dimer_method.traj',
                              logfile='log') as dim_rlx:
            #print(1)
            dim_rlx.run(fmax=0.02,steps=200)
            traj = Trajectory('dimer_method.traj')
            if not os.path.exists("traj"):
              os.system("mkdir traj")
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
               # atoms.calc=DP(model="../graph.000.pb")
                #atoms.write("0.traj")
                if box ==[]:
                    box= np.array(atoms.cell).reshape([1,-1])[0]
                else:
                   box=np.append(box,np.array(atoms.cell).reshape([1,-1])[0])
                 #print(box)
                 #print(np.array(atoms.positions))
                if coor==[]:
                    coor=np.array(atoms.positions).reshape([1,-1])[0]
                else:
                    coor=np.append(coor,np.array(atoms.positions).reshape([1,-1])[0])
              #  if count%2==0:
              #   print(atoms)
                #import numpy as np
                forces= np.zeros([len(atoms),3])#atoms.get_forces()
                write_lamptrj(atoms,forces,count)
                 
                 
                count+=1
          #  if (200-count>1):
           #         for k in range(5):
            #            atoms=traj[-1]
             #           atoms.calc=DP(model="../graph.000.pb")
              #          forces= atoms.get_forces()
               #         write_lamptrj(atoms,forces,count+k)
           # f1.close()
            #os.system("cp model_devi_tot.out model_devi.out")
            #print(len(coor))
            #print(len(box))
           # print(atype)
            coor=np.array([coor])
            box=np.array([box])
            from deepmd.infer import DeepPot as DPot
            graphs = [DPot("../graph.000.pb"), DPot("../graph.001.pb"), DPot("../graph.002.pb"), DPot("../graph.003.pb")]
            model_devi = calc_model_devi(coor, box, atype, graphs,fname="model_devi.out",frequency=1)
            print(1)  # view(atoms)
    # Analyze atoms



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fix", type=bool,default=False, help="if fix")
    parser.add_argument("index", type=int,default=-1, help="tag<index will be fixed")
    parser.add_argument("trans",type=float,default=0.1,help="max trans")
    parser.add_argument("center",type=int,default=0,help="displace center")
    parser.add_argument("radius",type=float,default=2,help="displace radius")
#    args = parser.parse_args()

    args = parser.parse_args
    main(args)

