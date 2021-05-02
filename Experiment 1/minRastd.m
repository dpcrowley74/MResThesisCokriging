function[model] = minRastd(seed, numhf, numlf, dim, n)
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
lb = -ones(dim,1);
ub = ones(dim,1);
% main loop runs a set number of times.
for i = 1:n
    % maximise the expected Improvement and find x associated with the maximum
    % expected improvement
    x = maxVarExpImp(model, 0, lb, ub);
    % evaluate its high fidelity point
    y = rastrigin(x,10000,1);
    % evaluate the low fidelity
    yc = rastrigin(x,1000,1);
    % store these values
    Xe = [Xe;x]; Xc = [Xc;x];
    Ye = [Ye;y]; Yc = [Yc;yc];
    % fit a new model
    model = Cokriging(Xe,Ye,Xc,Yc);
end
end
