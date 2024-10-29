%% Q4: Computes social planner's objective function
% Jack Dunbar
% Due: October 31, 2024

function obj_val = SP_objective(vX, vAlpha, mOmega, vLambda)
    [m , n] = size(mOmega);
    mX = reshape(vX, m, n);

    % Fast but inscrutable vectorization
    vUtility = @(mX, vAlpha, mOmega) (mX.^(1+mOmega) ./ (1+mOmega))' * vAlpha;
    obj_val = vLambda' * vUtility(mX, vAlpha, mOmega);
end
