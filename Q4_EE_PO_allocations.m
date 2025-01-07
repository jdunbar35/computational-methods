%% Q4: Solves for Pareto optimal allocations in an endowment economy
% Jack Dunbar
% October 31, 2024

clear variables; close all; clc

%% 1. m = n = 3, no heterogeneity

m = 3;      % # of rows/goods
n = 3;      % # of columns/people

% alpha - good-specific weights
vAlpha = linspace(1, 1, m)';

% omega - agent-good-specific utility function parameters

a = 0;      % Varies by good
b = 0;      % Varies by agent
c = -2;     % Base value

vGoods = (1:m)';
vAgents = (1:n)';

mOmega= a * vGoods + b * vAgents' + c; % mOmega(i, j) = a * i + b * j + c;

% e - endowments
mE = ones(m, n);

% lambda - Pareto weights
vLambda = linspace(1, 2, n)';
vLambda = vLambda / sum(vLambda); % Normalize to 1

mX_PO_1 = PO_allocation(mE, vAlpha, mOmega, vLambda);
%fprintf("1. m = n = 3, no heterogeneity\n")
%print_matrix(mX_PO_1, 2)

%% 2. m = n = 10, heterogeneity

m = 10;      % # of rows/goods
n = 10;      % # of columns/people

% alpha - good-specific weights
vAlpha = linspace(1, 1, m)';

% omega - agent-good-specific utility function parameters

a = -0.1;   % Varies by good
b = 0;      % Varies by agent
c = -2;     % Base value

vGoods = (1:m)';
vAgents = (1:n)';

mOmega= a * vGoods + b * vAgents' + c;

% e - endowments
mE = ones(m, n);

% lambda - Pareto weights
vLambda = linspace(1, 2, n)';
vLambda = vLambda / sum(vLambda); % Normalize to 1

mX_PO_3 = PO_allocation(mE, vAlpha, mOmega, vLambda);
%fprintf('\nm = n = 10, heterogeneity\n')
%print_matrix(mX_PO_3, 2)
