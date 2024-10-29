function [ineq, eq] = VF_constraints(v, parameters, end_states, given)
    ineq = [];

    c = v(1);
    l = v(2);
    kp = v(3);
    i = v(4);

    delta = parameters(4);
    psi = parameters(5);

    k = end_states(1);
    im = end_states(2);

    w = given(1);
    r = given(2);

    eq(1) = (1-tau)*w*l + r*k - c - i;
    eq(2) = (1-delta)*k + (1 - psi*(i/im - 1)^2)*i - kp;
end

