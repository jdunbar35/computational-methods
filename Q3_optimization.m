%% Q3: Compares the performance of optimization methods on the Rosenbrock function 
% Jack Dunbar
% October 31, 2024

%% 1. Setup

clear variables; close all; clc

% Every method gets the same parameters
tol = 1e-6;
max_iter = 10000;
v0 = [0; 0];

% Rosenbrock function
f = @(v) 100 * (v(2) - v(1)^2)^2 + (1 - v(1))^2;
grad_f = @(v) [-400*v(1)*(v(2)-v(1)^2) - 2*(1-v(1)); ...
               200*(v(2)-v(1)^2)];
hess_f = @(v) [-400*(v(2)-v(1)^2) - 800*v(1)^2 + 2, -400*v(1); ...
               -400*v(1), 200];
f_with_grad = @(v) deal(f(v), grad_f(v));

% Times
times = zeros(1, 4);

%% 2. Newton-Raphson

tic;

iter = 0;
v_curr = v0;

while iter < max_iter
    iter = iter + 1;

    grad_curr = grad_f(v_curr);
    hess_curr = hess_f(v_curr);

    % Update based on second-order Taylor approximation
    v_new = v_curr - (hess_curr \ grad_curr);

    if norm(grad_curr) < tol
        break;
    end

    v_curr = v_new;
end

v_NR = v_new;

times(1) = toc;
%fprintf("Newton-Raphson does %d iterations\n", iter)

%% 3. BFGS - built into Matlab

tic;
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off');
[v_bfgs, val_bfgs] = fminunc(f, v0, options);
times(2) = toc;

%% 4.  Steepest descent

tic;

iter = 0;
v_curr = v0;

while iter < max_iter
    iter = iter + 1;
    
    grad_curr = grad_f(v_curr);

    % Direction
    d = -grad_curr / norm(grad_curr);

    % Code for line search
    %g = @(alpha) f(v_curr + alpha*d);
    %alpha = fminunc(g, 0.1);

    % Step size
    alpha = 0.999 ^ iter; % Faster than line search

    % Update in direction of steepest descent
    v_new = v_curr + alpha * d;

    if norm(grad_curr) < tol
        break;
    end

    v_curr = v_new;
end

v_steep = v_new;

times(3) = toc;
%fprintf("Steepest Descent does %d iterations\n", iter)

%% 5. Conjugate descent

tic;

% Initialize first iteration
iter = 1;
grad_prev = grad_f(v0);
d = -grad_prev / norm(grad_prev); % Direction
v_curr = v0 + 0.1 * d;

while iter < max_iter
    iter = iter + 1;
    
    grad_curr = grad_f(v_curr);
    
    alpha = 0.999 ^ iter; % Runs faster than line search

    % Update direction
    beta = (grad_curr'*grad_curr) / (grad_prev'*grad_prev);
    d = -grad_curr + beta*d;
    d = d / norm(d);
    
    % Code for line search
    %g = @(alpha) f(v_curr + alpha*d);
    %options = optimoptions('fminunc', 'Display', 'off');
    %alpha = fminunc(g, 0.1, options);

    % Update guess
    v_new = v_curr + alpha * d;

    if norm(grad_curr) < tol
        break;
    end

    grad_prev = grad_curr;
    v_curr = v_new;
end

v_conj = v_new;

times(4) = toc;
%fprintf("Conjugate Descent does %d iterations\n", iter)
%fprintf("Times for methods are the following\n")
%print_matrix(times, 4)