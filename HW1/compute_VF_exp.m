function mVF_exp = compute_VF_exp(mVF_curr, mTrans)
    [nZ, nTau, nK, nIm] = size(mVF_curr);
    
    mVF_exp = zeros(nZ, nTau, nK, nIm);

    for iIm = 1:nIm
        for iK = 1:nK
            VF_states = mVF_curr(:, :, iK, iIm);
            for iTau = 1:nTau
                for iZ = 1:nZ
                    mTrans_state = mTrans(iZ, iTau, :, :);
                    mVF_exp(iZ, iTau, iK, iIm) = sum(VF_states .* mTrans_state, 'all');
                end
            end
        end
    end
end

