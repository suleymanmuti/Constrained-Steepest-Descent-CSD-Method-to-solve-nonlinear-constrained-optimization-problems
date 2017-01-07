clc; clear all; clc;

% define the problem
[f_obj, g_b, g_e, h_b, h_e, x, e1, e2, x_init, R_init, gamma] = define_problem();

% Determine numbers of each type of variable.
n = length(x); % number variables in objective function
p = length(h_b); % number of equality constraints
m = length(g_b); % number of inequality constraints

% Create extra variables in the Lagrange function
% Populate needed number of Lagrange multipliers
v = sym('v',[p 1]); % resulting from equality constraints.
u = sym('u',[m 1]); % resulting from inequality constraints.
% Populate needed number of slack variables 
s = sym('s',[m 1]); % resulting from inequality constraints.

% Variables in QP subproblem.
d = sym('d',[n 1]);

% Merge the variables of the Langrange function in a vector for
% differentiation.
kx = vertcat(x,v,u,s);
kd = vertcat(d,v,u,s); % contains all variables.
kl = length(kd); % number of all variables.


syms a;

dp = x_init.';
R = R_init;


norm_of_grad = e2 +1; 


% Gradients of the cost function, equality constraints, and inequality
% constraints.
grad_f = gradient(f_obj, x);

for i_iter = 1:p 
    grad_h_b(:, i_iter) = gradient(h_b, x);
end

for i_iter = 1:m 
    grad_g_b(:, i_iter) = gradient(g_b(i_iter), x).';
end
 


k = 1; % design point number.  k = 0 is initialdesign point.

% enter into a loop if initial design does not satisfy the condition of
% length of the gradient being smaller than epsilon.
while norm_of_grad > e2 && k < 5
    % condition was true, we need a new design point.
   
    
    % Values of the objective function, constraints, and gradients at the
    % current design point.
    f_value = subs(f_obj, x, dp);
    c = subs(grad_f, x, dp);
    
    h_value_max = 0;
    for i_iter = 1:p
        h_value(i_iter) = subs(h_b(i_iter), x, dp);
        h_value_max = max(h_value(i_iter));
    end
    
    g_value_max = 0;
    for i_iter = 1:m
        g_value(i_iter) = subs(g_b(i_iter), x, dp);
        g_value_max = max(g_value(i_iter));
    end
    
    V = max(0, max(h_value_max, g_value_max));
    
    
    if p > 0
        grad_h_b_value = subs(grad_h_b, x, dp);
    end
    
    if m > 0
        grad_g_b_value = subs(grad_g_b, x, dp);
    end
    
    
    

    
    % Defineth QP subproblem.
    
    f_lin = c.'*d + 0.5*d.'*d;

    h_e_lin = h_e;
    h_lin = h_b;
    if p > 0
        for i_iter = 1:m
            h_lin(i_iter) = h_value(i_iter) + grad_h_b_value(:,i_iter).'*d;
        end
    end
    
    g_e_lin = g_e;
    g_lin = g_b;
    if m > 0
        for i_iter = 1:m
            g_lin(i_iter) = g_value(i_iter) + grad_g_b_value(:,i_iter).'*d;
        end
    end
   
    [d_opt_sorted, f_lin_value] = opt_hw1(f_lin, n, h_lin, h_e_lin, g_lin, g_e_lin);
    
    d_num = d_opt_sorted(1:n);
    v_num = d_opt_sorted(n+1:p);
    u_num = d_opt_sorted(n+p+1:n+p+m);
    s_num = d_opt_sorted(n+p+m+1:n+p+m+m);
    
    
    norm_of_grad = norm(double(d_num));
    
    r = 0;
    if p > 0
        r = r + sum(abs(v_num));
    end
    if m > 0
        r = r + sum(u_num);
    end
    
    R = max(R, r);
   
    
    % Inexact line search
    
    mu = 0.5;
    phi = f_value + R*V;
	beta = gamma*norm_of_grad.^2;
    
    j = 0;    
    cond = false;
    while cond == false && j < 10
        
        t = (mu).^(j);
        %
        dp_next = dp.' + t*d_num.';
        f_value_next =  subs(f_obj, x, dp_next.');
        
        h_value_max = 0;
        
        
        f_value = subs(f_obj, x, dp_next.');
        
        for i_iter = 1:p
            h_value(i_iter) = subs(h_b(i_iter), x, dp_next.');
            h_value_max = max(h_value);
        end
        
        
        g_value_max = 0;
        for i_iter = 1:m
            g_value(i_iter) = subs(g_b(i_iter), x, dp_next.');
            g_value_max = max(g_value);
        end
        
        V_next = max(0, max(h_value_max, g_value_max));
        
        phi_next = f_value_next + R*V_next;
        
        if double(phi_next) < double(phi - t*beta)
            cond = true;
            a = t;
            dp = dp_next.';
            break;
        end
        
        j = j+1;
        
    end
    
    k = k + 1;
end

summary(f_obj, h_b, h_e, g_b, g_e, dp, f_value, n, p, m);

fprintf('Results:\n\n');

fprintf('\nOptimum point:\n');
for i_iter = 1:n
    fprintf('%f*\t', double(dp(i_iter)));
end
fprintf('\n');

fprintf('\nObjective function''s value at the optimum point:\nf* = %f\n',double(f_value));
fprintf('\n\nSee "summary_and_results.txt" for details...\n\n');