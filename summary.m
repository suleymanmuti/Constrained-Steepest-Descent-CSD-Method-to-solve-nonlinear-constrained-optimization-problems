function summary(f, h_b, h_e, g_b, g_e, dp, f_value, n, p, m)
% Store the summary and results of the problem in a text file.
fileID = fopen('summary_and_results.txt','wt');
fprintf(fileID,'*** Welcome to Constrained Steepest Descent (CSD) Method Solver ***\n\nAuthor: Suleyman Muti');
fprintf(fileID,'\n\n\nHere is the summary and results of the defined optimization problem:\n');
fprintf(fileID,'\n\nObjective function: f = %s\n',char(f));
fprintf(fileID,'\nSubject to:\n');
if p~= 0
    for i_iter = 1:p
        fprintf(fileID, '%s = %s\n',char(h_b(i_iter)), num2str(h_e(i_iter)));
    end
end
if m~= 0
    for j_iter = 1:m
        fprintf(fileID, '%s <= %s\n',char(g_b(j_iter)), num2str(g_e(j_iter)));
    end
end

fprintf(fileID,'\n\n\nOptimum point:\n');
for i_iter = 1:n
    fprintf(fileID, '%f*\t', double(dp(i_iter)));
end
fprintf(fileID, '\n');

fprintf(fileID, '\n\nObjective function''s value at the optimum point:\nf* = %f',double(f_value));
fclose(fileID);

end