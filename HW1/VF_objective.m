function obj_val = VF_objective(v, mVF_exp, vK, vIm, parameters, given)
    c = v(1);
    l = v(2);
    kp = v(3);
    i = v(4);

    %v

    beta = parameters(1);
    gamma = parameters(2);

    g = given(3);

    flow_val = utility(c, g, l, gamma);
    cont_val = beta * mVF_exp_interp(mVF_exp, kp, i, vK, vIm);
    obj_val = flow_val + cont_val;
end

