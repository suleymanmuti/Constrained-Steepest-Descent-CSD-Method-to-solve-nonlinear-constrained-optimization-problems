function [f, g_b, g_e, h_b, h_e, x, e1, e2, x_init, R_init, gamma] = define_problem()

% Example 13.4
n = 2; % Enter the number of variable in the objective function here.
x = sym('x',[n 1]); % variables of objective function
f = x(1).^2 + x(2).^2 - 3*x(1)*x(2);  


% Inequality constraints
g_b = [(1./6)*x(1).^2 + (1./6)*x(2).^2 - 1, -1*x(1), -1*x(2)];
g_e = [0, 0, 0];

% Equality constraints
h_b = [];
h_e = [];

% initial design
x_init = [1, 1];


R_init = 10;
gamma = 0.5;

e1 = 0.001;
e2 = 0.001;

end