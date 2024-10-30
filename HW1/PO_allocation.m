%% Q4: Solves for Pareto optimal allocations in an endowment economy
% Jack Dunbar
% October 31, 2024

function mX_opt = PO_allocation(mE, vAlpha, mOmega, vLambda)
    [m, n] = size(mE);
    
    % Find the allocation that solves the social planner's problem
    vX0 = ones(m*n, 1);
    obj = @(vX) -SP_objective(vX, vAlpha, mOmega, vLambda);
    cons = @(vX) SP_constraints(vX, mE);
    options = optimoptions('fmincon', 'Display', 'off');
    
    vX_opt = fmincon(obj, vX0, [], [], [], [], [], [], cons, options);
    
    mX_opt = reshape(vX_opt, m, n);
end

