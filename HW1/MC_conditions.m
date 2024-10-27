function conds = MC_conditions(vP, vAlpha, mOmega, mE)
    n = size(mOmega, 2);

    vP = [1; vP];

    % figure out demand as a function of prices and m
    vMu0 = ones(n, 1);

    mX_func = @(vMu) ((vP * vMu') ./ vAlpha) .^ (1 ./ mOmega);
    %mX_func(vMu0)
    BC_conds = @(vMu) vP' * (mE - mX_func(vMu));
    
    vMu_opt = fsolve(BC_conds, vMu0);
    %vMu_opt
    mX_opt = mX_func(vMu_opt);
    mX_opt

    conds = sum(mE - mX_opt, 2);
end

