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

c_ss = v_ss(1);
g_ss = v_ss(2);
l_ss = v_ss(3);
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

% Full state transition
% [Z, Tau, Z', Tau']
mOnes = ones(nZ, nTau, nZ, nTau);
mTau_trans = mOnes .* reshape(mTau_trans, 1, nTau, 1, nTau);
mZ_trans = mOnes .* reshape(mZ_trans, nZ, 1, nZ, 1);
mTrans = mTau_trans .* mZ_trans;

% Capital grid
nK = 5;
vK = linspace(k_ss * 0.7, k_ss * 1.3, nK)';

% Investment yesterday (i-) grid
nIm = 3;
vIm = linspace(i_ss * 0.5, i_ss * 1.5, nIm)';

% Take filthy initial guesses...
v0 = v_ss;
[mPol_c, mPol_l, mPol_kp, mPol_i, mW_curr, mR_curr, mVF_curr] = initial_guesses(v0, vZ, vTau, vK, vIm, parameters);

% Initialize new stuff (values don't matter)
mW_new = zeros(nZ, nTau, nK, nIm);
mW_tmp = zeros(nZ, nTau, nK, nIm);
mR_new = zeros(nZ, nTau, nK, nIm);
mR_tmp = zeros(nZ, nTau, nK, nIm);
mVF_new = zeros(nZ, nTau, nK, nIm);

% Wage, interest rate iteration
max_diff_w_r = 10;
tol_w_r = 10e-2;

% Value function iteration
tol_VF = 10e-2;

% Loop
tic;
iter_w_r = 0;
while (max_diff_w_r > tol_w_r)
    iter_w_r = iter_w_r + 1;

    max_diff_VF = 10;
    iter_VF = 0;

    while (max_diff_VF > tol_VF)
        iter_VF = iter_VF + 1;

        mVF_exp = compute_VF_exp(mVF_curr, mTrans); % FIX THIS
        
        for iZ = 1:nZ
            z = vZ(iZ);
            for iTau = 1:nTau
                tau = vTau(iTau);
                for iK = 1:nK
                     k = vK(iK);
                    for iIm = 1:nIm
                        im = vIm(iIm);
        
                        mVF_exp_curr = squeeze(mVF_exp(iZ, iTau, :, :));
        
                        states = [z, tau, k, im];
        
                        w = mW_curr(iZ, iTau, iK, iIm);
                        r = mR_curr(iZ, iTau, iK, iIm);
                        g = g_ss; % Doesn't affect HH choices. Just here for transparency
                        given = [w, r, g];
            
                        % [c, l, k', i]
                        % Initial guess = the policy functions from last iteration
                        v0 = [mPol_c(iZ, iTau, iK, iIm), mPol_l(iZ, iTau, iK, iIm), ...
                              mPol_kp(iZ, iTau, iK, iIm), mPol_i(iZ, iTau, iK, iIm)];
        
                        obj = @(v) -VF_objective(v, mVF_exp_curr, vK, vIm, parameters, given);
                        cons = @(v) VF_constraints(v, parameters, states, given);
        
                        % Add lb and ub
                        lb = [0, 0, vK(1), vIm(1)];
                        ub = [1000, 1000, vK(nK), vIm(nIm)];
        
                        options = optimoptions('fmincon', 'Display', 'off');
                        
                        % [c, l, k', i]
                        % DO THIS WAY BETTER
                        v_opt = fmincon(obj, v0, [], [], [], [], lb, ub, cons, options);
        
                        c = v_opt(1);
                        l = v_opt(2);
                        kp = v_opt(3);
                        i = v_opt(4);
                        
                        % Parameters
                        beta = parameters(1);
                        gamma = parameters(2);
                        alpha = parameters(3);
        
                        % Calculate new wages...
                        w = (1-alpha) * exp(z) * (k / l)^(alpha);
                        r = alpha * exp(z) * (k / l)^(alpha-1);
        
                        g = (1 - tau) * w * l;
        
                        flow_value = utility(c, g, l, gamma);
                        cont_value = beta * mVF_exp_interp(mVF_exp_curr, kp, i, vK, vIm);
        
                        mVF_new(iZ, iTau, iK, iIm) = flow_value + cont_value;
                        mPol_c(iZ, iTau, iK, iIm) = c;
                        mPol_l(iZ, iTau, iK, iIm) = l;
                        mPol_kp(iZ, iTau, iK, iIm) = kp;
                        mPol_i(iZ, iTau, iK, iIm) = i;

                        % Store w, r just in case...
                        mW_tmp(iZ, iTau, iK, iIm) = w;
                        mR_tmp(iZ, iTau, iK, iIm) = r;
                    end
                end
            end
        end
        
        diff_VF = mVF_new - mVF_curr;
        max_diff_VF = max(abs(diff_VF(:)));
        
        if (max_diff_VF < tol_VF)
            break;
        end
        
        mVF_curr = mVF_new;
    end

    % Update... DONT DO THIS SO HARSHLY
    mW_new = mW_tmp;
    mR_new = mR_tmp;

    diff_w = mW_new - mW_curr;
    diff_r = mR_new - mR_curr;
    
    max_diff_w_r = max(max(abs(diff_w(:))), max(abs(diff_r(:))));
    
    if (max_diff_w_r < tol_w_r)
        break;
    end
    
    mW_curr = mW_new;
    mR_curr = mR_new;
end
toc;