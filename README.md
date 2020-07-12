# Woofer Julia Simulation

## Overview
This repository contains the simulation of a QP based MPC controller which runs inside of MuJoCo. The current implementation is a fixed gait sequence with linearized rigid body dynamics model without leg dynamics. We use OSQP and Parametron in order to solve the QP at each time step. The Raibert heuristic is used in order to plan footsteps.

## Installation for Simulation
1. Acquire a license for MuJoCo at http://mujoco.org/. You can get a free trial of the professional license for a month, or with a student account, a free year.

2. Save the license ```mjkey.txt``` somewhere and set the environment variable ```MUJOCO_KEY_PATH``` to that location. One way to set the environment variable is through your bash profile. On a mac this is done by adding the line
```
export MUJOCO_KEY_PATH=[YOUR PATH]/mjkey.txt
```
to your ~/.bash_profile.

This can be done by executing:
```
echo 'export MUJOCO_KEY_PATH=[YOUR PATH]/mjkey.txt' >> ~/.bash_profile
```

3. Install Julia by visiting https://julialang.org/downloads/

4. Clone this repository.

5. Enter the julia REPL through the command line.

6. Enter the Package Manager by typing the `]` key and activate the environment with
```shell
pkg> activate .
```

7. Install MuJoCo.jl while still in the Package Manager with  
```shell
pkg> add https://github.com/klowrey/MuJoCo.jl
pkg> build MuJoCo
```

8. Then instantiate the environment, which should download all other dependencies.
```shell
pkg> instantiate
```
## Run Simulation
1. When inside the WooferQP repo, `cd src/Simulator`
2. Enter the Julia REPL inside of the folder by typing `julia` (may need to alias path to julia).
2. Run
```julia
include("WooferSim.jl")
Main.WooferSim.simulate()
```
3. The MuJoCo simulator should then pop up in a new window with various interactive options. Press space to start the simulation.
- Click and drag with the left mouse button to orbit the camera, and with the right mouse button to pan the camera.
- To perturb the robot, double click on the body you want to perturb, then hold Control and click and drag with the mouse. Using the left mouse button will apply a rotational torque while the right button will apply a translational force.
- Press space to play or pause the simulation.
