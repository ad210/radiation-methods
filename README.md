# radiation methods
This repository is a collection of several programs used to simulate radiation by electrons in various fields

## field-solver
An integration code to produce electron trajectories for various fields


## Clara2
Clara2 is a code borrowed from [another project](https://github.com/ComputationalRadiationPhysics/clara2).

After generating a trajectory with the field-solver, run
```
./executable
./process_data
```
to generate the spectrum of your trajectory 

This spectrum can be displayed with by calling 
```
../tools/plotRadiation ./my_spectrum_all_000.dat --dataExtend [min_angle] max
```


## splice-script
A collection of scripts to study the time evolution of radiation 


