function obj_val = SP_objective(vX, vAlpha, mOmega, vLambda)
    [m , n] = size(mOmega);
    mX = reshape(vX, m, n);
    vUtility = @(mX, vAlpha, mOmega) (mX.^(1+mOmega) ./ (1+mOmega))' * vAlpha;
    obj_val = vLambda' * vUtility(mX, vAlpha, mOmega);
end
