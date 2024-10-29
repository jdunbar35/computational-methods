%% Q5: Compute equilibrium prices in an endowment economy
% Jack Dunbar
% Due: October 31, 2024

clear variables; close all; clc

m = 3;      % # of rows/goods
n = 3;      % # of columns/people

% alpha - good-specific weights
vAlpha = linspace(2, 1, m)';

% omega - agent-good-specific utility function parameters

a = -0.1;      % Varies by good
b = -0.1;      % Varies by agent
c = -2;     % Base value

vGoods = (1:m)';
vAgents = (1:n)';

mOmega= a * vGoods + b * vAgents' + c; % mOmega(i, j) = a * i + b * j + c;

% e - endowments
mE = ones(m, n);

% Solving for prices
vP0 = ones(m-1, 1);

conds = @(vP) MC_conditions(vP, vAlpha, mOmega, mE);
% Levenberg-Marquardt works on non-square systems
options = optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');
vP_opt = fsolve(conds, vP0, options);

vP_opt = [1; vP_opt];