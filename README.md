## Computational Methods

This repository will contain my first homework for Computation Methods in fall 2024.

### Files
1. `Computational_Methods_HW1.pdf`: My solutions
1. `HW1_handout.pdf`: The assignment
1. `Q2_intergration.m`: Comparison of integration techniques
    - `integrate_utility.m`: Performs integration
1. `Q3_optimization.m`: Compares the performance of optimization methods
1. `Q4_EE_PO_allocations.m`: Solves for Pareto optimal allocations in an endowment economy
    - `PO_allocation.m`: Does the solving
        - `SP_objective.m`: Objective function of the social planner
        - `SP_constraints.m`: Social planner resource constraint
1. `Q5_EE_prices.m`: Computes equilibrium prices in an endowment economy
    - `MC_conditions.m`: Market clearing conditions for the endowment economy
1. `Q6_value_function_interation.m`: Value Function Iteration for Neoclassicial Growth Model with Distortionary Taxation
    - `SP_ss_conditions.m`: Social planner steady state conditions
    - `ss_conditions.m`: Decentralized steady state conditions
    - `initial_guesses.m`: Initial guesses for value function iteration
    - `compute_VF_exp.m`: Computes household expected value function
    - `VF_objective.m`: Household value function
    - `VF_constraints.m`: Constraints for the household problem
    - `utility.m`: Household utility function
    - `mVF_exp_interp.m`: Interpolates household expected value function
1. `print_matrix.m`: prints matrices with some # of decimal points