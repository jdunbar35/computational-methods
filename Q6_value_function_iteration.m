%% Q6: Value Function Iteration for Neoclassicial Growth Model with Distortionary Taxation
    % 6.3 takes ~40 seconds to run on a 5x5 capital-lagged investment grid with 0.1 tolerances
% Jack Dunbar
% October 31, 2024

clear variables; close all; clc

%% 6.1 Social Planner Steady State with z = 0, tau = 0.25

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

l_SP_ss = v_SP_ss(3);
k_SP_ss = v_SP_ss(4);

MPL_SP_ss = (1-0.33) * (k_SP_ss/l_SP_ss)^(-0.33);
MPK_SP_ss = 0.33 * (k_SP_ss/l_SP_ss)^(0.33-1);

%print_matrix(v_SP_ss', 2);
%print_matrix([MPL_SP_ss; MPK_SP_ss], 2);

%% 6.2 Decentralized Steady State with z = 0, tau = 0.25

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

w_ss = (1-0.33) * (k_ss/l_ss)^(-0.33);
r_ss = 0.33 * (k_ss/l_ss)^(0.33-1);

%fprintf("\n")
%print_matrix(v_ss', 2);
%print_matrix([w_ss; r_ss], 2);

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

% Capital (k) grid
nK = 5;
vK = linspace(k_ss * 0.7, k_ss * 1.3, nK)';

% Investment yesterday (i-) grid
nIm = 5;
vIm = linspace(i_ss * 0.5, i_ss * 1.5, nIm)';

% Take initial guesses based on deterministic steady states
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
tol_w_r = 1e-1;

% Value function iteration
tol_VF = 1e-1;

% Loop over wages and interest rates
tic;
iter_w_r = 0;
while (max_diff_w_r > tol_w_r)
    iter_w_r = iter_w_r + 1;

    max_diff_VF = 10;
    iter_VF = 0;

    while (max_diff_VF > tol_VF)
        iter_VF = iter_VF + 1;

        % Compute the expectation of the value function
        mVF_exp = compute_VF_exp(mVF_curr, mTrans);
        
        % Loop over states and update the value function at each
        % EFFICIENCY: Some vectorization could make this much faster
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

                        lb = [0, 0, vK(1), vIm(1)];
                        ub = [1000, 1000, vK(nK), vIm(nIm)];
        
                        options = optimoptions('fmincon', 'Display', 'off');
                        
                        % [c, l, k', i]
                        % EFFICIENCY: The code spends most of its time in here
                        v_opt = fmincon(obj, v0, [], [], [], [], lb, ub, cons, options);
        
                        c = v_opt(1);
                        l = v_opt(2);
                        kp = v_opt(3);
                        i = v_opt(4);

                        beta = parameters(1);
                        gamma = parameters(2);
                        alpha = parameters(3);
        
                        % Calculate new wages
                        w = (1-alpha) * exp(z) * (k / l)^(alpha);
                        r = alpha * exp(z) * (k / l)^(alpha-1);
        
                        g = (1 - tau) * w * l;
        
                        flow_value = utility(c, g, l, gamma);
                        % EFFICIENCY: Shouldn't be interpolating every time like this
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
        
        % Check for convergence
        diff_VF = mVF_new - mVF_curr;
        max_diff_VF = max(abs(diff_VF(:)));
        
        if (max_diff_VF < tol_VF)
            break;
        end
        
        mVF_curr = mVF_new;
    end

    % Update wage and interest rates
    % EFFICIENCY: Harsh updating may not be ideal
    mW_new = mW_tmp;
    mR_new = mR_tmp;

    % Check for convergence
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


%% 6.3: Plot Value Function, Capital Policy, and Wages

close all
lw = 1.5;

% Value function
fig = figure;
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile;
    plot(vK, squeeze(mVF_new(1, 2, :, 3)), 'LineWidth', lw); hold on;
    plot(vK, squeeze(mVF_new(2, 2, :, 3)), 'LineWidth', lw); 
    plot(vK, squeeze(mVF_new(3, 2, :, 3)), 'LineWidth', lw); 
    plot(vK, squeeze(mVF_new(4, 2, :, 3)), 'LineWidth', lw); 
    plot(vK, squeeze(mVF_new(5, 2, :, 3)), 'LineWidth', lw); 
    title("Value Function v (i^{-} and \tau fixed)"); xlabel("k");
    legend({sprintf('z = %.2f', vZ(1)), sprintf('z = %.2f', vZ(2)), ...
        sprintf('z = %.2f', vZ(3)), sprintf('z = %.2f', vZ(4)), ...
        sprintf('z = %.2f', vZ(5))}, 'Location', 'southeast');

nexttile;
    plot(vIm, squeeze(mVF_new(1, 2, 3, :)), 'LineWidth', lw); hold on;
    plot(vIm, squeeze(mVF_new(2, 2, 3, :)), 'LineWidth', lw); 
    plot(vIm, squeeze(mVF_new(3, 2, 3, :)), 'LineWidth', lw); 
    plot(vIm, squeeze(mVF_new(4, 2, 3, :)), 'LineWidth', lw); 
    plot(vIm, squeeze(mVF_new(5, 2, 3, :)), 'LineWidth', lw); 
    title("Value Function v (k and \tau fixed)"); xlabel("i^{-}");
    legend({sprintf('z = %.2f', vZ(1)), sprintf('z = %.2f', vZ(2)), ...
        sprintf('z = %.2f', vZ(3)), sprintf('z = %.2f', vZ(4)), ...
        sprintf('z = %.2f', vZ(5))}, 'Location', 'southeast');

nexttile;
    plot(vK, squeeze(mVF_new(3, 1, :, 3)), 'LineWidth', lw); hold on;
    plot(vK, squeeze(mVF_new(3, 2, :, 3)), 'LineWidth', lw); 
    plot(vK, squeeze(mVF_new(3, 3, :, 3)), 'LineWidth', lw); 
    title("Value Function v (i^{-} and z fixed)"); xlabel("k");
    legend({sprintf('\\tau = %.2f', vTau(1)), sprintf('\\tau = %.2f', vTau(2)), ...
        sprintf('\\tau = %.2f', vTau(3))}, 'Location', 'southeast');

nexttile;
    plot(vIm, squeeze(mVF_new(3, 1, 3, :)), 'LineWidth', lw); hold on;
    plot(vIm, squeeze(mVF_new(3, 2, 3, :)), 'LineWidth', lw); 
    plot(vIm, squeeze(mVF_new(3, 3, 3, :)), 'LineWidth', lw);
    title("Value Function v (k and z fixed)"); xlabel("i^{-}");
    legend({sprintf('\\tau = %.2f', vTau(1)), sprintf('\\tau = %.2f', vTau(2)), ...
        sprintf('\\tau = %.2f', vTau(3))}, 'Location', 'southeast');

fig.Position = [100, 100, 900, 500];
saveas(fig, "figures/value.png")

% Capital policy
fig = figure;
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile;
    plot(vK, squeeze(mPol_kp(1, 2, :, 3)), 'LineWidth', lw); hold on;
    plot(vK, squeeze(mPol_kp(2, 2, :, 3)), 'LineWidth', lw); 
    plot(vK, squeeze(mPol_kp(3, 2, :, 3)), 'LineWidth', lw); 
    plot(vK, squeeze(mPol_kp(4, 2, :, 3)), 'LineWidth', lw); 
    plot(vK, squeeze(mPol_kp(5, 2, :, 3)), 'LineWidth', lw); 
    title("k^{\prime} policy function (i^{-} and \tau fixed)"); xlabel("k");
    legend({sprintf('z = %.2f', vZ(1)), sprintf('z = %.2f', vZ(2)), ...
        sprintf('z = %.2f', vZ(3)), sprintf('z = %.2f', vZ(4)), ...
        sprintf('z = %.2f', vZ(5))}, 'Location', 'southeast');

nexttile;
    plot(vIm, squeeze(mPol_kp(1, 2, 3, :)), 'LineWidth', lw); hold on;
    plot(vIm, squeeze(mPol_kp(2, 2, 3, :)), 'LineWidth', lw); 
    plot(vIm, squeeze(mPol_kp(3, 2, 3, :)), 'LineWidth', lw); 
    plot(vIm, squeeze(mPol_kp(4, 2, 3, :)), 'LineWidth', lw); 
    plot(vIm, squeeze(mPol_kp(5, 2, 3, :)), 'LineWidth', lw); 
    title("k^{\prime} policy function (k and \tau fixed)"); xlabel("i^{-}");
    legend({sprintf('z = %.2f', vZ(1)), sprintf('z = %.2f', vZ(2)), ...
        sprintf('z = %.2f', vZ(3)), sprintf('z = %.2f', vZ(4)), ...
        sprintf('z = %.2f', vZ(5))}, 'Location', 'southeast');

fig.Position = [100, 100, 900, 300];
saveas(fig, "figures/policy_kp.png")

% Wages
fig = figure;
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile;
    plot(vK, squeeze(mW_new(1, 2, :, 3)), 'LineWidth', lw); hold on;
    plot(vK, squeeze(mW_new(2, 2, :, 3)), 'LineWidth', lw); 
    plot(vK, squeeze(mW_new(3, 2, :, 3)), 'LineWidth', lw); 
    plot(vK, squeeze(mW_new(4, 2, :, 3)), 'LineWidth', lw); 
    plot(vK, squeeze(mW_new(5, 2, :, 3)), 'LineWidth', lw); 
    title("Wages w (i^{-} and \tau fixed)"); xlabel("k");
    legend({sprintf('z = %.2f', vZ(1)), sprintf('z = %.2f', vZ(2)), ...
        sprintf('z = %.2f', vZ(3)), sprintf('z = %.2f', vZ(4)), ...
        sprintf('z = %.2f', vZ(5))}, 'Location', 'southeast');

nexttile;
    plot(vIm, squeeze(mW_new(1, 2, 3, :)), 'LineWidth', lw); hold on;
    plot(vIm, squeeze(mW_new(2, 2, 3, :)), 'LineWidth', lw); 
    plot(vIm, squeeze(mW_new(3, 2, 3, :)), 'LineWidth', lw); 
    plot(vIm, squeeze(mW_new(4, 2, 3, :)), 'LineWidth', lw); 
    plot(vIm, squeeze(mW_new(5, 2, 3, :)), 'LineWidth', lw); 
    title("Wages w (k and \tau fixed)"); xlabel("i^{-}");
    legend({sprintf('z = %.2f', vZ(1)), sprintf('z = %.2f', vZ(2)), ...
        sprintf('z = %.2f', vZ(3)), sprintf('z = %.2f', vZ(4)), ...
        sprintf('z = %.2f', vZ(5))}, 'Location', 'southeast');

fig.Position = [100, 100, 900, 300];
saveas(fig, "figures/wages.png")