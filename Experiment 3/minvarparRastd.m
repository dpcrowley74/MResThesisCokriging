function[model,f, i] = minvarparRastd(seed, numhf, numlf, dim, n, parallel, hfidelity, lfidelity)
% this is a script that returns a model that approximates the Total
% Varaince function and determins the minimum point via bayesian
% optimisation over one level over on level of fidelity.
%
% Input Variables:
%   seed - Is the seed for the random number generator to insure consistancy
%   numhf - Is the number of high fidelity points
%   numlf - Is the number of low fidleity points
%   dim - Is the number of dimensions the problem is based on
%   n - Is the number of steps for the bayesian optimisation
%   parallel - number of points found before the model is retrained
%
% Output Varaibles:
%
%   model - is a data structure used in the pred function, one can
%   determine the optimum by either optimising the pred function or
%   by using min(model.Ye) and location model.Xe(model.Ye == min(model.Ye))
rng(seed);
% As in this iteration we will only be considering points from the high
% fidelity
% f is zero as we are only considering the expected improvement for the
% high fidelity function
% cheap points
Xc = initialdesign(numlf,dim,-1,1);
% Expensive Points using the slection function from Surrogates Toolbox
% Version 3.0 by Felipe A. C. Viana
Xe = srgtsDOESubSample(Xc,numhf);


%evaluate the expensive function
Ye = rastrigin(Xe,10000,1);
% evaluate the cheap function
Yc = rastrigin(Xc,1000,1);
% fit an initial Cokriging model
model = Cokriging(Xe,Ye,Xc,Yc);
% set paramters for the maxVarExpImp
lb = ones(dim,1)*-1;
ub = ones(dim,1)*1;
lf = numlf - numhf;
i = 1;
h = 1;
% main loop runs a set number of times.
while h < n && i < 2*n
    % maximise the expected Improvement and find x associated with the maximum
    % expected improvement
    j = 1;
    k = 1;
    [x, fi, Imp] = maxVarExpImp(model,1,lb,ub, hfidelity, lfidelity);
    lfid = 0;
    while k < parallel
        %parallel loop for finding parallel points
        if fi(j) == 0
            %storing the found points in the correct location
            model.Ye = [model.Ye; pred(x(j,:),model,0)];
            model.Xe = [model.Xe; x(j,:)];
            model.Yc = [model.Yc; pred(x(j,:),model,1)];
            model.X = [model.X; x(j,:);];
            k = k+1;
            % breaking out of the loop
            if k == parallel, break; end
            lfid = lfidelity;
        else
            %again for storing in the correct location
            model.Yc = [model.Yc(1:lf); pred(x(j,:),model,1);model.Yc(lf+1:end)];
            model.X = [model.X(1:lf,:); x(j,:);model.X(lf+1:end,:)];
            Imp(j) = Imp(j) -lfid + lfidelity;
            % keeping track of the low fidelity points
            lf = lf+1;
            lfid = lfid + lfidelity;
        end
        j = j+1;
        % adjust the model to the new points
        model = buildmodel(model.Yc,model.X,model.T,model.p,model.Ye, model.Xe, model.r, model.Te, model.pd);
        [xtemp, fid, I] = maxVarExpImp(model,1,lb,ub, hfidelity, lfid);
        % keeping track of the points and associate fidelity's
        Imp = [Imp;I];
        x = [x;xtemp];
        fi = [fi;fid];
    end
    % parallel stopping condition taking the average
    if sum(Imp) >= 0
        f(i) = sum(fi);
        model = buildmodel(Yc,Xc,model.T,model.p,Ye, Xe, model.r, model.Te, model.pd);
        break;
    end
    y = rastrigin(x(fi == 0,:),10000,1);
    h = h + sum(fi == 0);
    % evaluate the low fidelity
    yc = rastrigin(x,1000,1);
    % store these values
    Xe = [Xe;x(fi == 0,:)];
    Ye = [Ye;y];
    Yc = [Yc;yc];
    Xc = [Xc;x];
    
    % record number of low fidelity
    f(i) = sum(fi);
    % Fit a new model
    model = Cokriging(Xe,Ye,Xc,Yc);
    i = i+1;
end
end