# DPGEN-DIMER
A modification of  [DPGEN](https://github.com/deepmodeling/dpgen) code, we add  the [dimer method](http://theory.cm.utexas.edu/henkelman/pubs/henkelman99_7010.pdf) as a method of model deviation, which could make our MLFF model predict the saddle points with high precision and few training data. 

## Requirements:
python>=3.8.0

dpgen

## Usage:
You just need to put  our `run.py` and  `arginfo.py` into your work directory, then prepare the param.json file and  machine.json file required by dpgen. You also need a file named `input.lammps`, no mattter what content is, we didn't use it but the code needs it.
To use the dimer method, you need to set the `model_devi_engine` to `dimer` in the machine.json file, and set the `command` to your command to run the python code `dimer.py`. For instance:
```
CUDA_VISIBLE_DEVICES=1 srun python /home/yinbc/test/dpgen/dpgen/Li-TiO2//dimer.py True  24  0.05  24  2
```
You can also use other model devi methods of dpgen, such as enhanced sampling. And you can also try to mix them.

## Tools:
`auto.py` in the example folder: A script which could automatically generate input files for multiple examples.
To use this, you need to put your initial vasp MD caculation folders together, and use [dpdata](https://github.com/deepmodeling/dpdata) to generate the initial data, for instance:
```
import dpdata
import numpy as np
vaspmd_dir = '.'
data=dpdata.LabeledSystem('OUTCAR')
data.to_deepmd_raw('deepmd')
data.to_deepmd_npy('deepmd')
```
Besides the calculation folders, put run.py, arginfo.py, dimer.py, Templates of INCAR, param.json, machine.json, input.lammps  together also. Then you can run `python auto.py`. Then you may adjust some of the parameters for each specific example.

`test_process.py`:This code could generate muliple search trajectory for the given initial configuration, and give the predicted energy. The energy are in the `ener.txt`, and we would get each trajectory folder. 

`test_saddle.py`: This code could generate muliple dimer search to get saddle points and corresponding barriers for the given initial configuration. The barriers are in the `barrier_sad.txt`. 

## PS:
In our code, you just need one POTCAR, not one POTCAR for each element. 