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

    % Define constants for variation
    a = -0.1;  % Varies by good
    b = 0;  % Varies by agent
    c = -1;  % Base value

    vGoods = (1:m)';
    vAgents = (1:n)';

    % mOmega(i, j) = a * i + b * j + c;
    mOmega= a * vGoods + b * vAgents' + c;

    % e - endowments
    mE = ones(m, n);

    % lambda - Pareto weights
    vLambda = linspace(1, 2, n)';
    vLambda = vLambda / sum(vLambda); % Normalize to 1

% Optimization
    vX0 = ones(m*n, 1);

    obj = @(vX) -SP_objective(vX, vAlpha, mOmega, vLambda);
    cons = @(vX) SP_constraints(vX, mE);

    options = optimoptions('fmincon', 'Display', 'off');

    vX_opt = fmincon(obj, vX0, [], [], [], [], [], [], cons, options);

    mX_opt = reshape(vX_opt, m, n);

%% Equilibrium Prices

close all; clc

% Solving for prices
vP0 = ones(m-1, 1)*2;

conds = @(vP) MC_conditions(vP, vAlpha, mOmega, mE);
vP_opt = fsolve(conds, vP0);
