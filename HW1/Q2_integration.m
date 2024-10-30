%% Q2: Comparison of integration techniques on a lifetime utility function
% Jack Dunbar
% October 31, 2024

clear variables; close all; clc
rng(123123);

% Setup
    % Parameters
    rho = 0.04;
    lambda = 0.02;
    
    % Grid
    lb = 0;
    ub = 100;
    n_bins = 10000000;
    vT = linspace(lb, ub, n_bins+1);
    bin_width = (ub - lb) / n_bins;
    
    % Function
    f = @(t) -exp(exp(-lambda*t) - 1 - rho*t);

% Midpoint quadrature
    int_mid = 0;
    
    for bin_num = 1:n_bins
        midpoint = lb + (bin_num - 0.5) * bin_width;
        int_mid = int_mid + f(midpoint);
    end
    
    int_mid = (int_mid / n_bins) * (ub-lb);

% Trapezoid quadrature
    int_trap = trapz(vT, f(vT));

% Simpson's quadrature (1/3 rule)
    int_simp = f(lb) + f(ub);
    
    for num_point = 2:n_bins
        point_value = lb + num_point*bin_width;
        if mod(num_point, 2) == 0
            int_simp = int_simp + 4 * f(point_value);
        elseif mod(num_point, 2) == 1
            int_simp = int_simp + 2 * f(point_value);
        end
    end
    
    int_simp = int_simp * (bin_width/3);

% Monte Carlo
    vT_random = lb + (ub - lb) * rand(n_bins, 1);
    f_random = f(vT_random);
    int_monte = (ub-lb) * mean(f_random);