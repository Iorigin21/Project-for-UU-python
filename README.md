# Project-for-UU-python:  Numerical solution of a 1D PN junction using Python
Task: Give the size of the PN diode, doping concentration, Apply voltage.Solve the field equations via Newton-Raphson method. Get electron, hole distribution. Draw energy band, IV curve
(1) cell1d.py,device.py used to build the PN diode(Like SDE in TCAD sentaurus)
(2) parameters.py contains the material parameters used in the simulation.
(3) models.py contain the equation need to solve
(4) gen_solve_matrix.py used to create and solve the Jacobian matrix. And NetwonSolver need to call this module
(5) NetwonSolver.py is the core part in this project. Which used to solve the Equations of mathematical physics
(6) Main.py is the program you need to run
