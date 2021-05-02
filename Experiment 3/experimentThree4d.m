ystar1 = ones(100,1);
xstar1 = ones(100,4);
ystar2 = ones(100,1);
xstar2 = ones(100,4);
ystar4 = ones(100,1);
xstar4 = ones(100,4);




parfor i = 1:100
    [model1,f1] = minvarparRastd(i, 12, 40, 4, 30, 1, 0.05, 0.005);
    x = minModelPred(model1,4,-1,1);
    y = rastrigin(x,10000);
    if y < min(model1.Ye)
        ystar1(i) = y;
        xstar1(i,:) = x;
    else
        [ystar1(i),indx] = min(model1.Ye);
        xstar1(i,:) = model1.Xe(indx,:);
    end
    [model2,f2] = minvarparRastd(i, 12, 40, 4, 30, 2, 0.05, 0.005);
    x = minModelPred(model2,4,-1,1);
    y = rastrigin(x,10000);
    if y < min(model2.Ye)
        ystar2(i) = y;
        xstar2(i,:) = x;
    else
        [ystar2(i),indx] = min(model2.Ye);
        xstar2(i,:) = model2.Xe(indx,:);
    end
    [model4,f4] = minvarparRastd(i, 12, 40, 4, 30, 4, 0.05, 0.005);
    x = minModelPred(model4,4,-1,1);
    y = rastrigin(x,10000);
    if y < min(model4.Ye)
        ystar4(i) = y;
        xstar4(i,:) = x;
    else
        [ystar4(i),indx] = min(model4.Ye);
        xstar4(i,:) = model4.Xe(indx,:);
    end
end
save experimentThree4d ystar1 xstar1 ystar2 xstar2 ystar4 xstar4