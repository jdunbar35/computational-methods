%% Q5: Compute equilibrium prices in an endowment economy
% Jack Dunbar
% October 31, 2024

clear variables; close all; clc

m = 10;      % # of rows/goods
n = 10;      % # of columns/people

% alpha - good-specific weights
vAlpha = linspace(1, 1, m)';

% omega - agent-good-specific utility function parameters

a = -0.1;   % Varies by good
b = -0;     % Varies by agent
c = -2;     % Base value

vGoods = (1:m)';
vAgents = (1:n)';

mOmega = a * vGoods + b * vAgents' + c; % mOmega(i, j) = a * i + b * j + c;

% e - endowments
mE = ones(m, n) .* (1:n); %

% Solving for prices (keep one fixed)
vP0 = ones(m-1, 1);

conds = @(vP) MC_conditions(vP, vAlpha, mOmega, mE);
% Levenberg-Marquardt works on non-square systems
options = optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');

vP_opt = fsolve(conds, vP0, options);
vP_opt = [1, reshape(vP_opt, 1, m-1)]; % Add back the price kept fixed

fprintf("m = n = 10, heterogeneity\n")
print_matrix(vP_opt, 2)