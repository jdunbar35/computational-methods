function [ineq, eq] = SP_constraints(vX, mE)
    ineq = [];

    [m , n] = size(mE);
    mX = reshape(vX, m, n);
    eq = sum(mE - mX, 2);
end

