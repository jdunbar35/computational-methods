%% Q6: Returns decentralized steady state conditions
% Jack Dunbar
% October 31, 2024

% Function: Create steady state optimality conditions for solving
function ss_conds = ss_conditions(v, parameters) 
    % Name inputs
    c = v(1);
    g = v(2);
    l = v(3);
    k = v(4);
    i = v(5);

    % Name parameters
    beta = parameters(1);
    %gamma = parameters(2);
    alpha = parameters(3);
    delta = parameters(4);
    %psi = parameters(5);
    z = parameters(6);
    tau = parameters(7);

    % Define some helper variables
    U_c = 1 / c;
    U_l = -l;
    %f = exp(z) * k^alpha * l^(1-alpha);
    f_l = (1-alpha) * exp(z) * (k/l)^alpha;
    f_k = alpha * exp(z) * (k/l)^(alpha-1);
    h = i;
    h_i = 1;
    h_im = 0;

    % Equation time
    ss_conds(1) = U_l + (1-tau) * f_l * U_c;
    ss_conds(2) = (beta * f_k / (1-beta*(1-delta))) * (beta * h_im + h_i) - 1;
    ss_conds(3) = (1-tau) * f_l * l + f_k * k - c - i;
    ss_conds(4) = (1-delta)*k + h - k;
    ss_conds(5) = tau * f_l * l - g;
end
