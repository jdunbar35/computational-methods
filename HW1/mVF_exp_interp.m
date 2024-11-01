%% Q6: Interpolates household expected value function
% Jack Dunbar
% October 31, 2024

function value = mVF_exp_interp(mVF_exp, kp, i, vK, vIm)
    k_lower = find(vK <= kp, 1, 'last');  % closest lower index for capital
    k_upper = find(vK >= kp, 1, 'first'); % closest upper index for capital
    
    i_lower = find(vIm <= i, 1, 'last');  % closest lower index for investment
    i_upper = find(vIm >= i, 1, 'first'); % closest upper index for investment

    % Extract the values for interpolation
    Q11 = mVF_exp(k_lower, i_lower);
    Q12 = mVF_exp(k_lower, i_upper);
    Q21 = mVF_exp(k_upper, i_lower);
    Q22 = mVF_exp(k_upper, i_upper);
    
    % Bilinear interpolation weights
    denom = (vK(k_upper) - vK(k_lower)) * (vIm(i_upper) - vIm(i_lower));
    if denom == 0
        value = mVF_exp(k_lower, i_lower); % Return nearest value if denom is zero
    else
        value = (1 / denom) * ...
            (Q11 * (vK(k_upper) - kp) * (vIm(i_upper) - i) + ...
             Q21 * (kp - vK(k_lower)) * (vIm(i_upper) - i) + ...
             Q12 * (vK(k_upper) - kp) * (i - vIm(i_lower)) + ...
             Q22 * (kp - vK(k_lower)) * (i - vIm(i_lower)));
    end
end

