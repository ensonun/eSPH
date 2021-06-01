# eSPH
eSPH is a simple, lightweight 2D SPH code written in MATLAB. The method is based on a high-order, low-dissipation Riemann solver SPH architecture.

It is part of Enson's final year project for a MEng degree at Imperial College London. Please refer to my [thesis]() for more information and references.

## Dependencies

The code uses the functions in this repository and MATLAB built-in functions only. After downloading, remember to put **all** member functions in the same directory as the ```eSPH.m```.

To enable parallel computing, please ensure that the [Parallel Computing Toolbox](https://uk.mathworks.com/products/parallel-computing.html) is installed in your local MATLAB.

The current version of the code is tested with MATLAB 2019a. Please report any conflicts with newer MATLAB versions.

## I/O
### Inputs

The code is run by calling the function ```eSPH("$FNAME.mat")```.

The input ```.mat``` file contains the followings (must be in exact names):

1. A  Nx11 double array ```fluid```:

| Entry             | Parameter                             |
| :---------------- | ------------------------------------- |
| ```fluid(:,1)```  | x-coordinate                          |
| ```fluid(:,2)```  | y-coordinate                          |
| ```fluid(:,3)```  | density                               |
| ```fluid(:,4)```  | mass (constant throughout simulation) |
| ```fluid(:,5)```  | pressure                              |
| ```fluid(:,6)```  | x-velocity                            |
| ```fluid(:,7)```  | y- velocity                           |
| ```fluid(:,8)```  | *(Reserved)*                          |
| ```fluid(:,9)```  | reference density of fluid, rho_0     |
| ```fluid(:,10)``` | artificial speed of sound, c_0        |
| ```fluid(:,11)``` | kinematic viscosity, nu               |

2. A  Nx4 double array ```wall```:

| Entry             | Parameter                         |
| :---------------- | --------------------------------- |
| ```wall(:,1)```  | x-coordinate                      |
| ```wall(:,2)```  | y-coordinate                      |
| ```wall(:,3)``` | x-velocity                        |
| ```wall(:,4)```  | y-velocity                       |

3. A function handle ```f``` for the body force

   ```f = @(x,y,t) [f_x; f_y]```

   *Note*: x, y are Nx1 arrays, so the output of ```f``` has to be Nx2. You can do this naively by ```f = @(x,y,t) [f_x; f_y] + 0*x'```.

4. A 1x8 double array ```settings```:

| Entry             | Parameter                                                    | Options                                                      |
| :---------------- | ------------------------------------------------------------ | :----------------------------------------------------------- |
| ```settings(1)``` | Kernel function                                              | 3 - Cubic b-spline<br>4 - Quartic b-spline<br>5 - 5th order Wendland |
| ```settings(2)``` | Kernel support radius, kh                                    | /                                                            |
| ```settings(3)``` | Gamma in equation of state (usually 7 for water, 1.4 for air) | /                                                            |
| ```settings(4)``` | Reconstruction scheme                                        | 1 - MUSCL piecewise linear <br>2 - MUSCL piecewise parabolic (not TVD)<br>3 - 3rd order WENO-Z<br>otherwise - piecewise constant |
| ```settings(5)``` | Riemann solver                                               | 0 - Classical SPH w/o dissipation <br>1 - Roe solver         |
| ```settings(6)``` | Switch for second derivative                                 | 1 - on<br>0 - off                                            |
| ```settings(7)``` | Simulation end time (non-dimensionalised by t_ref)           | /                                                            |
| ```settings(8)``` | CFL number (must be <=1)                                     | /                                                            |

5. A double ```dt_save```, indicating the time between which each output file is written
6. A string ```dir_name```, specifying the name of the directory where the output files are stored in

### Outputs

The code outputs 

1. a series of ```SPHout_$NSTEP.mat``` files according to the inputs
2. a ```.txt``` file containing the settings and runtime

## Contribution

All kinds of discussions and developments are more than welcomed. However, since the project is finished, please do not expect a significant update to the code.

## Licence

This source code is licensed under the MIT license, which can be found in the [LICENSE](LICENSE) file.
