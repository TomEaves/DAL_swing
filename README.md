# DAL_swing
Computation of minimal disturbances for desynchronisation in the swing equation (a second-order Kuramoto model of power grids) using the Direct-Adjoint-Looping method (DAL).

## power_grid_DAL.m
Usage: [y0,gains,resids] = power_grid_DAL(E0,T,G,alpha,Ps,K,th0s,y0)

Searches for optimally growing initial conditions over a time T of energy E0 in the vicinity of initial condition y0, in the graph G with parameters alpha (damping coefficient) and K (coupling strength). The distribution of power in the grid is Ps. 

The steady state (synchronised) phases at each node are inputted as th0s.
Outputs a new initial condition y0 along with a list of gains (kinetic energy at time T) and residuals for each initial condition considered during the iteration.

## autoDAL.m
Uses power_grid_DAL.m to 'automatically' implement the algorithm described in [...]. Reduces the initial condition energy E0 systematically until convergence to a minimal disturbance is acheived.

## G1.m
A simple script for use with ode45 (in autoDAL.m) to simulate the swing equation. Currently set up for a particular simple 4-node graph.
