%% Q4: Computes slack in the social planner's resource constrain
% Jack Dunbar
% Due: October 31, 2024

function [ineq, eq] = SP_constraints(vX, mE)
    ineq = [];

    [m , n] = size(mE);
    mX = reshape(vX, m, n);
    
    eq = sum(mE - mX, 2);
end

