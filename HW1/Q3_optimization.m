% Description...
% Jack Dunbar
% Due: October 31, 2024

% Setup
    clear variables; close all; clc
    
    tol = 1e-6;
    
    % f
    f = @(v) 100 * (v(2) - v(1)^2)^2 + (1 - v(1))^2;
    grad_f = @(v) [-400*v(1)*(v(2)-v(1)^2) - 2*(1-v(1)); ...
                   200*(v(2)-v(1)^2)];
    hess_f = @(v) [-400*(v(2)-v(1)^2) - 800*v(1)^2 + 2, -400*v(1); ...
                   -400*v(1), 200];
    f_with_grad = @(v) deal(f(v), grad_f(v));
    v0 = [0; 0];

% Newton-Raphson
    max_iter = 1000;
    iter = 0;
    v_curr = v0;
    
    while iter < max_iter
        iter = iter + 1;
    
        grad_curr = grad_f(v_curr);
    
        hess_curr = hess_f(v_curr);
    
        v_new = v_curr - (hess_curr\grad_curr);
    
        if norm(v_new - v_curr) < tol
            break;
        end
    
        v_curr = v_new;
    end
    
    v_NR = v_new;

% BFGS
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off');
    [v_bfgs, val_bfgs] = fminunc(f, v0, options);

% Steepest descent
    max_iter = 100000;
    iter = 0;
    v_curr = v0;
    
    while iter < max_iter
        iter = iter + 1;
        
        grad_curr = grad_f(v_curr);
    
        d = -grad_curr / norm(grad_curr);
    
        %g = @(alpha) f(v_curr + alpha*d);
        %alpha0 = 0.01;
        %alpha = fminunc(g, alpha);
        alpha = 1 / iter; % Faster than line search for the same quality
    
        v_new = v_curr + alpha * d;
    
        if norm(grad_curr) < tol
            break;
        end
    
        v_curr = v_new;
    end
    
    v_steep = v_new;

% Conjugate descent
    max_iter = 1000;
    iter = 0;
    v0 = [0; 0];
    v_curr = v0;
    
    alpha = 0.01;
    
    while iter < max_iter
        iter = iter + 1;
        
        grad_curr = grad_f(v_curr);
        
        %alpha = 0.99^iter; % Runs faster than line search
    
        if iter == 1
            d = -grad_curr;
            g = @(alpha) f(v_curr + alpha*d);
            alpha = fminunc(g, alpha);
            v_new = v_curr + alpha * d;
        else
            beta = (grad_curr'*grad_curr) / (grad_prev'*grad_prev);
            vd = -grad_curr + beta*d;
            g = @(alpha) f(v_curr + alpha*d);
            alpha = fminunc(g, alpha);
            v_new = v_curr + alpha * vd;
        end
    
    
        if norm(grad_curr) < tol
            break;
        end
    
        grad_prev = grad_curr;
        v_curr = v_new;
    end
    
    v_conj2 = v_new;