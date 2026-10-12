
# 1D slab neutron transport solver
````
  - discrete ordinates method (SN)
  - linear discontinuous galerkin spatial discretization scheme
  - one-group
  - steady state
  - fixed source, eigenvalue, or critical size calculation
  - source iteration
  - non-uniform cell spacing
  - vacuum, reflective, or specified incoming flux BC
````
<img width="300" height="300" alt="psibar" src="https://github.com/user-attachments/assets/c75a110d-c06b-4729-ba82-0eec0fce56b7" />

# How to build
simply change to src directory and type make (assuming you have a standard compiler).
```
cd ./src
make
```
This creates the exectuable `project-1.exe`

# How to run
Switch to a test directory and use:
````
../src/project-1.exe input-filename.inp
````
Example input files are included in each directory along with their corresponding output (excluding Test B).
