%% Q2: Comparison of integration techniques on a lifetime utility function
% Jack Dunbar
% October 31, 2024

clear variables; close all; clc
rng(123123);

% Try 3 different # bins/draws
for n_bins = [1000, 100000, 10000000]
    if n_bins == 1000
        mIntegral = integrate_Rosenbrock(n_bins);
    else
        mIntegral = [mIntegral, integrate_utility(n_bins)];
    end
end

%fprintf("Lifetime utility function integrals, time to compute\n")
%print_matrix(mIntegral, 4);