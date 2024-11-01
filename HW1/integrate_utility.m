%% Integrates lifetime utility function for a # of bins/draws
    % Compares midpoint, trapezoid, Simpson's, and Monte Carlo
    % Returns value and time for each
% Jack Dunbar
% October 31, 2024

function mIntegrals = integrate_utility(n_bins)
    % Parameters
    rho = 0.04;
    lambda = 0.02;
    
    % Grid
    lb = 0;
    ub = 100;
    vT = linspace(lb, ub, n_bins+1);
    bin_width = (ub - lb) / n_bins;
    
    % Lifetime utility function
    f = @(t) -exp(exp(-lambda*t) - 1 - rho*t);

    % Initialize output
    % [midpoint; trapizoid; Simpson's; Monte Carlo]
    mIntegrals = zeros(4, 2);

    % Vectorized midpoint quadrature
    tic;
    mIntegrals(1, 1) = sum(f(lb + (1:n_bins - 0.5) * bin_width)) / n_bins * (ub-lb);
    mIntegrals(1, 2) = toc;
    
    % Trapezoid quadrature -- built into Matlab
    tic;
    mIntegrals(2, 1) = trapz(vT, f(vT));
    mIntegrals(2,2) = toc;
    
    % Simpson's quadrature (1/3 rule)
    tic;
    int_simp = f(lb) + f(ub);
    for num_point = 2:n_bins
        point_value = lb + num_point*bin_width;
        if mod(num_point, 2) == 0
            int_simp = int_simp + 4 * f(point_value);
        elseif mod(num_point, 2) == 1
            int_simp = int_simp + 2 * f(point_value);
        end
    end
    
    mIntegrals(3, 1) = int_simp * (bin_width/3);
    mIntegrals(3, 2) = toc;
    
    % Monte Carlo
    tic;
    vT_random = lb + (ub - lb) * rand(n_bins, 1);
    f_random = f(vT_random);

    mIntegrals(4, 1) = (ub-lb) * mean(f_random);
    mIntegrals(4, 2) = toc;
end