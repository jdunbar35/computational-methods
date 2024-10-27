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

% Solve for steady state
% [c, g, l, k, i]
options = optimoptions('fsolve', 'Display', 'off');
conds = @(v) ss_conditions(v, parameters);
v_ss = fsolve(conds, v0, options);

k_ss = v_ss(4);
i_ss = v_ss(5);

%% 6.3 Value Function Iteration with a Fixed Grid

% Set parameters
% [beta, gamma, alpha, delta, psi]
parameters = [0.97, 0.2, 0.33, 0.1, 0.05];

% Productivity (z) process
vZ = [-0.0673; -0.0336; 0; 0.0336; 0.0673];
mZ_trans = [0.9727, 0.0273, 0.0000, 0.0000, 0.0000;
            0.0041, 0.9806, 0.0153, 0.0000, 0.0000;
            0.0000, 0.0082, 0.9837, 0.0082, 0.0000;
            0.0000, 0.0000, 0.0153, 0.9806, 0.0041;
            0.0000, 0.0000, 0.0000, 0.0273, 0.9727];

% Labor tax rate (tau) process
vTau = [0.2; 0.25; 0.3];
vTau_trans = [0.9, 0.1, 0;
              0.05, 0.9, 0.05;
              0, 0.1, 0.9];

% Capital grid
nK_grid = 250;
wK_grid = 0.3;
vK_grid = linspace(k_ss * (1-wK_grid), k_ss * (1+wK_grid), nK_grid)';

% Investment yesterday (i-) grid
nIm_grid = 50;
wIm_grid = 0.5;
vIm_grid = linspace(i_ss * (1-wIm_grid), i_ss * (1+wIm_grid), nIm_grid)';

% Required matrices and vectors

mOutput           = zeros(nGridCapital,nGridProductivity);
mValue    = zeros(nGridCapital,nGridProductivity);
mValue_new = zeros(nGridCapital,nGridProductivity);
mPolicy  = zeros(nGridCapital,nGridProductivity);
m = zeros(nGridCapital,nGridProductivity);


% Value function iteration

maxDifference = 10;
tolerance = 10e-6;
iteration = 0;

while (maxDifference>tolerance)  
    
    % Compute continuation
    expectedValueFunction = mValueFunction*mTransition';
    
    for nProductivity = 1:nGridProductivity
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for nCapital = 1:nGridCapital
                        
            valueHighSoFar = -1000.0;
            capitalChoice  = vGridCapital(1);
            
            for nCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital
                
                consumption = mOutput(nCapital,nProductivity)-vGridCapital(nCapitalNextPeriod);
                valueProvisional = (1-bbeta)*log(consumption)+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);
            
                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = vGridCapital(nCapitalNextPeriod);
                    gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                    break; % We break when we have achieved the max
                end    
                  
            end
            
            mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar;
            mPolicyFunction(nCapital,nProductivity) = capitalChoice;
            
        end
        
    end
    
    maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
    mValueFunction = mValueFunctionNew;
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
    end
           
end