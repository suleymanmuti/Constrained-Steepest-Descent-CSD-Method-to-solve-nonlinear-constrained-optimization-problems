function [opt_sorted, f_value] = opt_hw1(f, n, h_b, h_e, g_b, g_e)

% Define symbolic variables.
d = sym('d',[n 1]); % variables of the objective function, namely x and y.


% Determine numbers of each type of variable.
n = length(d); % number variables in objective function
p = length(h_b); % number of equality constraints
m = length(g_b);% number of inequality constraints

% Create extra variables in the Lagrange function
% Populate needed number of Lagrange multipliers
v = sym('v',[p 1]); % resulting from equality constraints.
u = sym('u',[m 1]); % resulting from inequality constraints.
% Populate needed number of slack variables 
s = sym('s',[m 1]); % resulting from inequality constraints.

% Calculate Lagrange function.
% Since different number of equality and inequality constraints can be
% submittted, calculated intwo steps.

% Step #1: Sum of equality portion in Lagrange equation.
l_e = sym(zeros(p,1)); 
for i_iter = 1:p
    l_e(i_iter) = v(i_iter)*( h_b(i_iter)-h_e(i_iter) );
end 

% Step #1: Sum of inequality portion in Lagrange equation.
l_i= sym(zeros(m,1));
for j_iter = 1:m
    l_i(j_iter) = u(j_iter)*( g_b(j_iter)-g_e(j_iter) + s(j_iter)^2 );
end 

% Lagrange function
l = f + sum(l_e)+ sum(l_i);

% Merge the variables of the Langrange function in a vector for
% differentiation.
k = vertcat(d,v,u,s); % contains all variables.
kl = length(k); % number of all variables.

% Symbolic differentiation and creation of Karush Kuhn Tucker equations.
sym kkt;
for i_iter = 1:kl
   kkt(i_iter) = (diff(l,k(i_iter))==0);
end

% Solve KKT equations
sol = solve(kkt);

% Extract values for candidate optimum points and all other variables.
for i_iter = 1:length(k)
    sol_val(:,i_iter) = subs(sol.(char(k(i_iter))));
end

[sr,sc] = size(sol_val); % get row and column numbers.

% Extract each variable's numeric value, only if that type of variable
% exists.
x_numeric = sol_val(1:sr,1:n); % x-type of variables.

v_lb = n+1; % lower bound for v in solutions matrix.
v_ub = n+p; % upper bound for v in solutions matrix.
if p~= 0
    v_numeric = sol_val(1:sr,n+1:n+p);  % v-type of variables.
end

u_lb = n+p+1; % lower bound for u in solutions matrix.
u_ub = n+p+m; % upper bound for u in solutions matrix.
if m ~= 0
    u_numeric = sol_val(1:sr,u_lb:u_ub);  % u-type of variables.
    s_numeric = sol_val(1:sr,u_lb+m:u_ub+m);  % s-type of variables.
end


% solution needs to be real number, clean numerical value vectors 
% from any violation.
for i_iter = 1: sr
    if m~= 0
        if ~isreal(u_numeric(i_iter,:)) || ~isreal(s_numeric(i_iter,:))
            x_numeric(i_iter,:)= 0;
            if p~= 0
                v_numeric(i_iter,:)= 0;
            end
            if m~= 0
                u_numeric(i_iter,:)= 0;
                s_numeric(i_iter,:)= 0;
            end
        end
    end
end

% if exist, u-type of variables need to be greater than or equal to zero.
% clean numerical value vectors if there is a violation.
for i_iter = 1: sr
    if m~= 0
        if any(u_numeric(i_iter,:)<0) || any(s_numeric(i_iter,:)<0)
            x_numeric(i_iter,:)= 0;
            if p~= 0
                v_numeric(i_iter,:)= 0;
            end
            if m~= 0
                u_numeric(i_iter,:)= 0;
                s_numeric(i_iter,:)= 0;
            end
        end
    end
end

% Sort numeric values, merge with existing variables into a vector.
% i.e., sorted candidate optimum.
opt_sorted = x_numeric;
if p~= 0
    opt_sorted = [opt_sorted, v_numeric];
end
if m~= 0
    opt_sorted = [opt_sorted, u_numeric, s_numeric];
end

% Clear any all zero rows and columns.
opt_sorted(all(opt_sorted==0,2),:)=[];
opt_sorted = opt_sorted';


% Calculate Hessian for the remaining candidate points.
h = hessian(f, d);
heig = eig(h);
[x_op, y_op] = size(opt_sorted); % get row and column numbers.
for i_iter = 1: y_op
    % Send in variables of objective function, 
    % and eigen values of the Hessian matrix to check for optimality.
    for j_iter = 1 : n
     x_candidate(j_iter) = opt_sorted(j_iter, i_iter);
    end
    
    copt = check_opt(heig, x_candidate', n); % return 0 if candidate violates Hessian condition.
    
    % if violates, delete that point from sorted candidate points.
    if copt == 0
        opt_sorted(:,i_iter)= []; % remove false candidates.
    end
end

% Calculate the value of the objective function at the optimum point
f_value = subs(f, k, opt_sorted);
% Calculate the calue of the Lagrange function at the optimum point
lagrange_value = subs(l, k, opt_sorted);

% Arrenge variables to substitude them into the Lagrange function.
del= k;
del(1,:)=[]; % remove variables of the objective function.
del(1,:)=[]; % leave behind only the extra variables.
op_del = opt_sorted;
op_del(1,:)=[];  % remove variables of the objective function.
op_del(1,:)=[]; % leave behind only the extra variables.


end