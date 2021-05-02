% Initialise the loop
s = zeros(100,1);
ystarvar = ones(100,1);
xstarvar = ones(100,2);
costvar = zeros(100,1);
cost = ones(100,1)*0.55;
ystar = ones(100,1);
xstar = ones(100,2);
highfi = zeros(100,1);

% For a parallel for loop run with prafor otherwise replace with a for loop
parfor i = 1:100
    % Variable fidelity optimisation method
    [modelvar, f] = minvarRastd(i,6,20,2,0.55,0.005,0.05);
    % Results
    x = minModelPred(modelvar,2,-1,1);
    y = rastrigin(x,10000);
    if y < min(modelvar.Ye)
        ystarvar(i) = y;
        xstarvar(i,:) = x;
    else
        [ystarvar(i),indx] = min(modelvar.Ye);
        xstarvar(i,:) = modelvar.Xe(indx,:);
    end
    highfi(i) = sum(1-f);
    costvar(i) = sum((1-f)*0.05 + 0.005);
    
    
    %single Fidlelity optimisation method
    model = minRastd(i,6,20,2,10);
    %results
    x1 = minModelPred(model,2,-1,1);
    y1 = rastrigin(x1,10000);
    if y1 < min(model.Ye)
        ystar(i) = y1;
        xstar(i,:) = x1;
    else
        [ystar(i),indx] = min(model.Ye);
        xstar(i,:) = model.Xe(indx,:);
    end
    if ystarvar(i) < ystar(i)
        s(i) = 1;
 
   end
end
% Save the results
save experiment12d s ystarvar xstarvar costvar ystar xstar cost model modelvar highfi
