function [mW, mR] = initial_w_r(v0, vZ, vTau, vK, vIm, parameters)
    % Dimensions
    nZ = length(vZ);
    nTau = length(vTau);
    nK = length(vK);
    nIm = length(vIm);

    mL_ss = zeros(nZ, nTau);
    
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
            
            mL_ss(iZ, iTau) = v_ss(3);
        end
    end

    % Meshgrid of relevant variables
    % Probably some of this is unecessary, but the cost is low and 4D matrices hurt my little brain
    mOnes = ones(nZ, nTau, nK, nIm);
    mZ = mOnes .* reshape(vZ, nZ, 1, 1, 1);
    mK = mOnes .* reshape(vK, 1, 1, nK, 1);
    mL_ss = mOnes .* reshape(mL_ss, nZ, nTau, 1, 1);

    % Firm's FOCs
    alpha = parameters(3);
    mW = (1-alpha) * exp(mZ) .* (mK ./ mL_ss) .^ (alpha);
    mR = alpha * exp(mZ) .* (mK ./ mL_ss) .^ (alpha-1);
end
