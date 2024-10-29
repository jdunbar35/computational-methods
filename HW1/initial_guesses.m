function [mC, mL, mK, mI, mW, mR, mVF] = initial_guesses(v0, vZ, vTau, vK, vIm, parameters)
    % Dimensions
    nZ = length(vZ);
    nTau = length(vTau);
    nK = length(vK);
    nIm = length(vIm);

    % Parameters
    beta = parameters(1);
    gamma = parameters(2);
    alpha = parameters(3);

    mC_ss = zeros(nZ, nTau);
    mG_ss = zeros(nZ, nTau);
    mL_ss = zeros(nZ, nTau);
    mK_ss = zeros(nZ, nTau);
    mI_ss = zeros(nZ, nTau);

    mVF_ss = zeros(nZ, nTau);
    
    for iZ = 1:nZ
        z = vZ(iZ);
        for iTau = 1:nTau
            tau = vTau(iTau);
            parameters_ss = [parameters, z, tau];

            % Solve for steady state
            % [c, g, l, k, i]
            options = optimoptions('fsolve', 'Display', 'off');
            conds = @(v) ss_conditions(v, parameters_ss);
            v_ss = fsolve(conds, v0, options);
            
            mC_ss(iZ, iTau) = v_ss(1);
            mG_ss(iZ, iTau) = v_ss(2);
            mL_ss(iZ, iTau) = v_ss(3);
            mK_ss(iZ, iTau) = v_ss(4);
            mI_ss(iZ, iTau) = v_ss(5);

            
            flow_value = utility(v_ss(1), v_ss(2), v_ss(3), gamma);
            mVF_ss(iZ, iTau) = flow_value / (1-beta);
        end
    end
    
    % Meshgrid
    mOnes = ones(nZ, nTau, nK, nIm);
    mC = mOnes .* reshape(mC_ss, nZ, nTau, 1, 1);
    mL = mOnes .* reshape(mL_ss, nZ, nTau, 1, 1);
    mK = mOnes .* reshape(mK_ss, nZ, nTau, 1, 1);
    mI = mOnes .* reshape(mI_ss, nZ, nTau, 1, 1);
    mVF = mOnes .* reshape(mVF_ss, nZ, nTau, 1, 1);

    % Do wages, interest rate

    % Meshgrid of relevant variables
    mZ = mOnes .* reshape(vZ, nZ, 1, 1, 1);
    mK_full = mOnes .* reshape(vK, 1, 1, nK, 1);

    % Firm's FOCs
    mW = (1-alpha) * exp(mZ) .* (mK_full ./ mL) .^ (alpha);
    mR = alpha * exp(mZ) .* (mK_full ./ mL) .^ (alpha-1);
end
