# eSPH
eSPH is a simple, lightweight 2D Riemann solver based high-order SPH code written in MATLAB.

It is part of Enson's final year project for a MEng degree at Imperial College London. Please refer to my [thesis]() for more information and references.

## Dependencies

The code uses the functions in this repository and MATLAB built-in functions only. After downloading, remember to put **all** member functions in the same directory as the ```mainSPH.m```.

To enable parallel computing, please ensure that the [Parallel Computing Toolbox](https://uk.mathworks.com/products/parallel-computing.html) is installed in your local MATLAB.

The current version of the code is tested with MATLAB 2019a. Please report any conflicts with newer MATLAB versions.

## I/O
### Inputs

The main code ```mainSPH.m``` takes in a ```.mat``` file which default name is ```input.mat```.

The input file contains the following variables:

- A  double array ```fluid```

  ```fluid(:,1)``` = x-coordinate

  ```fluid(:,2)``` = y-coordinate

  ```fluid(:,3)``` = density

  ```fluid(:,4)``` = mass (constant throughout simulation)

  ```fluid(:,5)``` = pressure

  ```fluid(:,6)``` = x velocity

  ```fluid(:,7)``` = y velocity

  ```fluid(:,8)``` = *(To be added)*

  ```fluid(:,9)``` = reference density of fluid, rho_0

  ```fluid(:,10)``` = artificial speed of sound, c_0

  ```fluid(:,11)``` = kinematic viscosity, nu

- A double array ```wall```

  ```wall(:,1)``` = x-coordinate

  ```wall(:,2)``` = y-coordinate

  ```wall(:,3)``` = x velocity

  ```wall(:,4)``` = y velocity

- A function handle ```f```

  ```f = @(x,y,t) [f_x; f_y]```

  *Note*: x, y are N x 1 arrays, so the output of ```f``` has to be N x 2

- A double array ```settings```

```settings(1:6)``` are parameters related to```timeDer.m```:

| Entry             | Parameter                                                    | Options                                                      |
| :---------------- | ------------------------------------------------------------ | :----------------------------------------------------------- |
| ```settings(1)``` | Kernel function type                                         | 2 - Wendland 5th order<br>3 - cubic spline                   |
| ```settings(2)``` | Smoothing length, h                                          | /                                                            |
| ```settings(3)``` | gamma in equation of state (usually 7 for water, 1.4 for air) | /                                                            |
| ```settings(4)``` | Reconstruction scheme                                        | 0 - piecewise constant <br>1 - MUSCL piecewise linear <br>2 - MUSCL piecewise parabolic <br>5 - WENO (Zhang, 2019) |
| ```settings(5)``` | Riemann solver type                                          | 0 - classical SPH w/o dissipation <br>1 - Roe solver         |
| ```settings(6)``` | Riemann solver limiter parameter                             | /                                                            |

```settings(7:9)``` are time integration related:

```settings(7)``` = time integration order (=2)
```settings(8)``` = simulation end time (non-dimensionalised by sqrt(g/H))
```settings(9)``` = max dt

- (Optional) A string ```file_name```

- (Optional) A string array ```catch_err``` *(To be added)*

### Outputs

The code outputs a series of ```.mat``` files according to the settings. The post-processing routine ```mat2vtu.m``` can be used to convert the ```.mat``` files to ```.vtu``` files readable to [ParaView](https://www.paraview.org/) or other visualisation software.

## Contribution

All kinds of discussions and developments are more than welcomed. However, since the project is finished, please do not expect a significant update to the code.

## Licence

This source code is licensed under the MIT license, which can be found in the [LICENSE](LICENSE) file.
