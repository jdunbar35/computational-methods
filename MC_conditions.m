%% Q5: Computes excess demand at a set of prices
% Jack Dunbar
% October 31, 2024

function conds = MC_conditions(vP, vAlpha, mOmega, mE)
    n = size(mOmega, 2);

    vP = [1; vP];

    % Demand as a function of the Lagrange multiplier mu
    mX_demand = @(vMu) ((vP * vMu') ./ vAlpha) .^ (1 ./ mOmega);

    % Find the Mu that statisfies the market clearing constraint
    vMu0 = ones(n, 1);
    BC_conds = @(vMu) vP' * (mE - mX_demand(vMu));
    options = optimoptions('fsolve', 'Display', 'off');
    vMu_opt = fsolve(BC_conds, vMu0, options);

    % Compute the allocation at that mu
    mX_opt = mX_demand(vMu_opt);

    % Return excess demand
    conds = sum(mE - mX_opt, 2);
end

