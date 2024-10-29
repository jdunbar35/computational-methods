clear variables; close all; clc

%% 6.1 Social Planner Steady State

% Set parameters
% [beta, gamma, alpha, delta, psi, z, tau]
parameters_ss = [0.97, 0.2, 0.33, 0.1, 0.05, 0, 0.25];

% Initial guess
% [c, g, l, k, i]
v0 = [1, 1, 1, 1, 1];

% Solve for steady state
% [c, g, l, k, i]
options = optimoptions('fsolve', 'Display', 'off');
conds = @(v) SP_ss_conditions(v, parameters_ss);
v_SP_ss = fsolve(conds, v0, options);

%% 6.2 Decentralized Steady State

% Set parameters
% [beta, gamma, alpha, delta, psi, z, tau]
parameters_ss = [0.97, 0.2, 0.33, 0.1, 0.05, 0, 0.25];

% Initial guess
% [c, g, l, k, i]
v0 = [1, 1, 1, 1, 1];

% Solve for steady state
% [c, g, l, k, i]
options = optimoptions('fsolve', 'Display', 'off');
conds = @(v) ss_conditions(v, parameters_ss);
v_ss = fsolve(conds, v0, options);

k_ss = v_ss(4);
i_ss = v_ss(5);

%% 6.3 Value Function Iteration with a Fixed Grid

% Set parameters
% [beta, gamma, alpha, delta, psi]
parameters = [0.97, 0.2, 0.33, 0.1, 0.05];

% Productivity (z) process
nZ = 5;
vZ = [-0.0673; -0.0336; 0; 0.0336; 0.0673];
mZ_trans = [0.9727, 0.0273, 0.0000, 0.0000, 0.0000;
            0.0041, 0.9806, 0.0153, 0.0000, 0.0000;
            0.0000, 0.0082, 0.9837, 0.0082, 0.0000;
            0.0000, 0.0000, 0.0153, 0.9806, 0.0041;
            0.0000, 0.0000, 0.0000, 0.0273, 0.9727];

% Labor tax rate (tau) process
nTau = 3;
vTau = [0.2; 0.25; 0.3];
mTau_trans = [0.9, 0.1, 0;
              0.05, 0.9, 0.05;
              0, 0.1, 0.9];

% Capital grid
nK = 250;
vK = linspace(k_ss * 0.7, k_ss * 1.3, nK)';

% Investment yesterday (i-) grid
nIm = 50;
vIm = linspace(i_ss * 0.5, i_ss * 1.5, nIm)';

% Initialize policy functions
mPol_c = zeros(nZ, nTau, nK, nIm);
mPol_l = zeros(nZ, nTau, nK, nIm);
mPol_kp = zeros(nZ, nTau, nK, nIm);
mPol_i = zeros(nZ, nTau, nK, nIm);

% Transition
% [Z, Tau, Z', Tau']
mOnes = ones(nZ, nTau, nZ, nTau);
mTau_trans = mOnes .* reshape(mTau_trans, 1, nTau, 1, nTau);
mZ_trans = mOnes .* reshape(mZ_trans, nZ, 1, nZ, 1);
mTrans = mTau_trans .* mZ_trans;

% Initial guess for wage, interest rate
v0 = v_ss; % Guess is deterministic steady state
[mW, mR] = initial_w_r(v0, vZ, vTau, vK, vIm, parameters);

% Also take a sick guess at the value function...
mVF_curr = zeros(nZ, nTau, nK, nIm); % initial_VF(VF0, vZ, vTau, parameters)
mVF_new = zeros(nZ, nTau, nK, nIm);

% Wage, interest rate iteration
max_diff_w_r = 10;
tolerance_w_r = 10e-6;
iteration_w_r = 0;

% Loop over wage, interest rate...

% Value function iteration
max_diff_VF = 10;
tolerance_VF = 10e-6;
iteration_VF = 0;

% Do value function loop...

% Taking your value function as given...

mVF_exp = compute_VF_exp(mVF_curr, mTrans);

for iZ = 1:nZ
    z = vZ(iZ);
    for iTau = 1:nTau
        tau = vTau(iTau);
        for iK = 1:nK
            k = vK(iK);
            for iIm = 1:nIm
                im = vIm(iIm);

                mVF_exp_curr = squeeze(mVF_exp(iZ, iTau, :, :));

                end_states = [k, im];

                w = mW(iZ, iTau, iK, iIm);
                r = mR(iZ, iTau, iK, iIm);
                g = 1; % Doesn't affect HH choices. Just here for transparency
                given = [w, r, g];
    
                % [c, l, k', i]
                % Initial guess = deterministic steady state...
                v0 = ones(4, 1);

                obj = @(v) -VF_objective(v, mVF_exp_curr, vK, vIm, parameters, given);
                cons = @(v) VF_constraints(v, parameters, end_states, given);

                options = optimoptions('fmincon', 'Display', 'off');

                % [c, l, k', i]
                v_opt = fmincon(obj, v0, [], [], [], [], lb, [], cons, options);

                c = v_opt(1);
                l = v_opt(2);
                kp = v_opt(3);
                i = v_opt(4);
                g = (1 - tau) * w * l;

                beta = parameters(1);
                gamma = parameters(2);

                flow_value = utility(c, g, l, gamma);
                cont_value = beta * mVF_exp_interp(mVF_exp_curr, kp, i, vK, vIm);

                mVF_new(iZ, iTau, iK, iIm) = flow_value + cont_value;
            end
        end
    end
end