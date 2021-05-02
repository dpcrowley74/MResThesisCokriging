% Setting up for loop
diff = zeros(100,1);
vardiff = zeros(100,1);
maxiter = zeros(100,1);
% parallel for loop if a single for loop is required use for
parfor i = 1:100
    % for single fidelity
    [model, benefit, iter] = minRastd(i,24,80,8,20,0.1);
    % saving results
    [m, miter] = max(benefit);
    if iter == 0
        % if the stopping condition is not triggered
        diff(i) = 21;
    else
        diff(i) = iter - miter;
    end
    % variable fidelity
    [modelvar, varbenefit, viter, maxiter(i)] = minvarRastd(i,24,80,8,2.2,0.1,0.01);
    [mv, mviter] = max(varbenefit);
    if viter == 0
        vardiff(i) = maxiter(i)+1;
    else
        vardiff(i) = viter - mviter;
    end
end
save experimentTwo8d diff vardiff maxiter