
# 1D slab neutron transport solver
````
  - discrete ordinates method (SN)
  - Lumped linear discontinuous spatial discretization scheme
  - one-group
  - steady state
  - fixed source, eigenvalue, or critical size calculation
  - source iteration
  - non-uniform cell spacing
  - vacuum, reflective, or specified incoming flux BC
  - nonlinear source acceleration using Quasi Diffusion (QD)
````
# How to build
simply change to src directory and type make (assuming you have a standard compiler).
```
cd ./src
make
```
This creates the exectuable `project-2.exe`

# How to run
Switch to a test directory and use:
````
../src/project-2.exe input-filename.inp
````
Example input files are included in each directory along with their corresponding output (excluding Test B).
