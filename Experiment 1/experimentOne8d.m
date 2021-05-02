s = zeros(100,1);
ystarvar = ones(100,1);
xstarvar = ones(100,8);
costvar = zeros(100,1);
cost = ones(100,1)*2.2;
ystar = ones(100,1);
xstar = ones(100,8);
highfi = zeros(100,1);


parfor i = 1:100
    [modelvar, f] = minvarRastd(i,24,80,8,2.2,0.01,0.1);
    x = minModelPred(modelvar,8,-1,1);
    y = rastrigin(x,10000);
    if y < min(modelvar.Ye)
        ystarvar(i) = y;
        xstarvar(i,:) = x;
    else
        [ystarvar(i),indx] = min(modelvar.Ye);
        xstarvar(i,:) = modelvar.Xe(indx,:);
    end
    costvar(i) = sum((1-f)*0.1 + 0.01);
    highfi(i) = sum(1-f);
    
    
    
    model = minRastd(i,24,80,8,20);
    x1 = minModelPred(model,8,-1,1);
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
save experimentOne8d s ystarvar xstarvar costvar ystar xstar cost model modelvar highfi
