function [return_opt] = check_opt(heig, x_candidate, n)

x = sym('x',[n 1]); % variables of objective function

return_opt = 1; % set to 1 assuming candidate point is optimum. 
for i_iter = 1: length(heig)
    check = subs(heig(i_iter),x, x_candidate);
    if check < 0
        return_opt = 0; % set to 0 if an eigenvalue is zero.
    end
end

end