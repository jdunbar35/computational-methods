clear variables; close all; clc

%% 6.1 Social Planner Steady State

% Set parameters
% [beta, gamma, alpha, delta, psi, tau, z]
parameters = [0.97, 0.2, 0.33, 0.1, 0.05, 0.25, 0];

% Initial guess
% [c, g, l, k, i]
v0 = [1, 1, 1, 1, 1];

% Solve for steady state
% [c, g, l, k, i]
options = optimoptions('fsolve', 'Display', 'off');
conds = @(v) SP_ss_conditions(v, parameters);
v_SP_ss = fsolve(conds, v0, options);

%% 6.2 Decentralized Steady State

% Set parameters
% [beta, gamma, alpha, delta, psi, tau, z]
parameters = [0.97, 0.2, 0.33, 0.1, 0.05, 0.25, 0];

% Initial guess
% [c, g, l, k, i]
v0 = [1, 1, 1, 1, 1];

v_ss = compute_ss(v0, parameters);

%% 6.3 Value Function Iteration with a Fixed Grid

% Set parameters
% [beta, gamma, alpha, delta, psi]
parameters = [0.97, 0.2, 0.33, 0.1, 0.05];

vZ = [-0.0673; -0.0336; 0; 0.0336; 0.0673];
mZ_trans = [0.9727, 0.0273, 0.0000, 0.0000, 0.0000;
               0.0041, 0.9806, 0.0153, 0.0000, 0.0000;
               0.0000, 0.0082, 0.9837, 0.0082, 0.0000;
               0.0000, 0.0000, 0.0153, 0.9806, 0.0041;
               0.0000, 0.0000, 0.0000, 0.0273, 0.9727];

