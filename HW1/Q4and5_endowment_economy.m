% Description
% Jack Dunbar
% Due: October 31, 2024

%% Pareto Optimal Allocations

clear variables; close all; clc

% Parameters
    m = 2;      % # of rows/goods
    n = 3;      % # of columns/people
    
    % alpha - the weights people place on goods
    vAlpha = linspace(1, 1, m)';

    % omega - utility function parameters
    % Preallocate the matrix
    mOmega = zeros(m, n);

    % Define constants for variation
    a = -0.1;  % Varies by good
    b = 0;  % Varies by agent
    c = -1;  % Base value

    % Loop through each agent and good to fill in omegas
    for j = 1:n
        for i = 1:m
            mOmega(i, j) = a * i + b * j + c;
        end
    end

    % e - endowments
    mE = ones(m, n);

    % lambda - Pareto weights
    vLambda = linspace(1, 2, n)';
    vLambda = vLambda / sum(vLambda); % Normalize to 1

% Optimization
    vX0 = ones(m*n, 1);

    obj = @(vX) -SP_objective(vX, vAlpha, mOmega, vLambda);
    lb = zeros(m*n, 1);
    cons = @(vX) SP_constraints(vX, mE);

    options = optimoptions('fmincon', 'Display', 'off');

    vX_opt = fmincon(obj, vX0, [], [], [], [], lb, [], cons, options);

    mX_opt = reshape(vX_opt, m, n);

%% Equilibrium Prices

close all; clc

% Solving for prices
vP0 = ones(m-1, 1)*2;

conds = @(vP) MC_conditions(vP, vAlpha, mOmega, mE);
vP_opt = fsolve(conds, vP0);
