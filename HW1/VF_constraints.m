%% Q6: Constraints for the household value function
% Jack Dunbar
% October 31, 2024

function [ineq, eq] = VF_constraints(v, parameters, states, given)
    ineq = [];

    c = v(1);
    l = v(2);
    kp = v(3);
    i = v(4);

    delta = parameters(4);
    psi = parameters(5);

    tau = states(2);
    k = states(3);
    im = states(4);

    w = given(1);
    r = given(2);

    eq(1) = (1-tau)*w*l + r*k - c - i;
    eq(2) = (1-delta)*k + (1 - psi*(i/im - 1)^2)*i - kp;
end

