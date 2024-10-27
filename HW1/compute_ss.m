function v_ss = compute_ss(v0, parameters)
    % Name parameters
    beta = parameters(1);
    gamma = parameters(2);
    alpha = parameters(3);
    delta = parameters(4);
    %psi = parameters(5);
    tau = parameters(6);
    z = parameters(7);

    % Solve for steady state
    % [c, g, l, k, i]
    options = optimoptions('fsolve', 'Display', 'off');
    conds = @(v) ss_conditions(v, parameters);
    v_ss = fsolve(conds, v0, options);
end