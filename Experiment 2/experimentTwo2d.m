% setting up for the loop
diff = zeros(100,1);
vardiff = zeros(100,1);
maxiter = zeros(100,1);
% change to for for single loop
parfor i = 1:100
    [model, benefit, iter] = minRastd(i,6,20,2,20,0.05);
    [m, miter] = max(benefit);
    if iter == 0
        diff(i) = 21;
    else
        diff(i) = iter - miter;
    end
    [modelvar, varbenefit, viter, maxiter(i)] = minvarRastd(i,6,20,2,1.1,0.05,0.005);
    [mv, mviter] = max(varbenefit);
    if viter == 0
        vardiff(i) = maxiter(i)+1;
    else
        vardiff(i) = viter - mviter;
    end
end
save experimentTwo2d diff vardiff maxiter
