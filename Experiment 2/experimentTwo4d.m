% setting up for loop, greater granular detail was stored for the four
% dimnesional case
stop = zeros(100,1);
stopvar = zeros(100,1);
peak = zeros(100,1);
peakvar = zeros(100,1);
maxiter = zeros(100,1);
maximum = zeros(100,1);
maxvar = zeros(100,1);
parfor i = 1:100
    [model, benefit, stop(i)] = minRastd(i,12,40,4,20,0.1);
    [maximum(i), peak(i)] = max(benefit);
    if stop(i) == 0
        stop(i) = 21;
    end
    [modelvar, varbenefit, stopvar(i), maxiter(i)] = minvarRastd(i,12,40,2,2.2,0.1,0.01);
    [maxvar(i), peakvar(i)] = max(varbenefit);
    if stopvar(i) == 0
        stopvar(i) = maxiter(i)+1;
    end
end
save experimentTwo4d stop stopvar maxiter peak peakvar maximum maxvar