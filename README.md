# radiation-methods
This repository is a collection of several programs used to simulate radiation by electrons in various fields

## field-solver
An integration code to produce electron trajectories for various fields. 

`ode_solver.py` is the generic integration code with no fields while `undulator.py` and `chicane.py` are copies with the specific fields already defined. 

There are two funtions which define the fields for all values of position and time:

> def E_fun(r,t,params):
>     x, y, z = r
>     q, m, c = params
>     Ex = 0.
>     Ey = 0.
>     Ez = 0.
>     return np.array([Ex,Ey,Ez])

The initial values defined in the main method allow the user to change the electron energy and initial direction. It is important to note that Clara2 assumes by default that the particle will be traveling along the x-axis. This can be changed in ```~/clara2-dev/src/all_directions.cpp``` if desired. The params vector is useful for passing variables about the fields to the field function.

The integrator runs for a default of 2000 time steps. This can be changed for finer resolution if the fields change too rapidly. Because Clara2 implements an interpelation scheme, it is unlikely you will need more than 2000 time steps for to resolve the radiation. 

`odeint` takes an optional argument `hmax=None` which defines the maximum time step size. This is important for fields that only change suddenly in small regions. To save computation time, odeint will take larger time steps if there is no change in a single time step. This can cause it to step over regions of field if ```hmax``` is not set. 

The code pads the start of the trajectory by ten times the length of the integrator. This zero padding is required by Clara2 because of the fourier analysis used in calculating the far field radiation. A factor of ten was found to be a good balance between resolution and computation speed. If the spectrum looks sharp or fractured, try a longer padding length.

Finally, Clara2 takes trajectories in .csv files of the form:
```
 x   y   z   Bx  By  Bz  t
--- --- --- --- --- --- ---
 0   0   0   0   0   0   0
...
 0   0   0   0   0   0   0
```

The position (x,y,z) is in meters, the velocity is in (m/s) per c, and time is in seconds.

`tools.py` is a collection of functions which are used in main for various purposes.

## Clara2
Clara2 is a code borrowed from [another project](https://github.com/ComputationalRadiationPhysics/clara2) and licensed under GPLv3+. See that page for a list of dependencies and licence and reference information. The version located in this repository has been edited slightly for this file structure. 

The simulation parameters can be edited in `settings.hpp`

> omega_max : The maximum frequency of the spectrum
> theta_max : The simulation looks from -theta -> +theta
> N_spectrum: Frequency resolution of the spectrum
> N_theta   : Angular resolution of the spectrum
> N_trace   : Number of particles (due to an OB1 error subtract 1 from your actual particle number)

Be sure to call the makefile any time `settings.hpp` is edited.

After generating a trajectory with the field-solver, run

> ./executable
> ./process_data

to generate the spectrum of your trajectory 

This spectrum can then be displayed with by calling 

> ../tools/plotRadiation ./my_spectrum_all_000.dat --dataExtend 0 [max_omega] [min_angle] [max_angle]

The flag --dataExtend allows you to label the axes with the proper angle and frequency

## splice-script
A collection of scripts to study the time evolution of radiation 


