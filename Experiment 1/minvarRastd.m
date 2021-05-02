function[model,f] = minvarRastd(seed, numhf, numlf, dim, totalcost, lfidcost, hfidcost)
% this is a script that returns a model that approximates the Rastrigin function and determins the minimum point via bayesian
% optimisation over one level over on level of fidelity.
%
% Input Variables:
%   seed - Is the seed for the random number generator to insure consistancy
%   numhf - Is the number of high fidelity points
%   numlf - Is the number of low fidleity points
%   dim - Is the number of dimensions the problem is based on
%   totalcost - the total cost of highfidelity computations
%   hfidcost - the cost of a single high fidelity evaluation
%
% Output Varaibles:
%
%   model - is a data structure used in the pred function, one can
%   determine the optimum by either optimising the pred function or
%   by using min(model.Ye) and location model.Xe(model.Ye == min(model.Ye))
%   f - is the number of high fidelity an low fidelity

% Setting default inputs for totalcost and hfidcost
if nargin < 7, hfidcost = 0.01; end
if nargin < 6, lfidcost = 0.1*hfidcost; end
if nargin < 5, totalcost = 0.33; end

rng(seed);
% As in this iteration we will only be considering points from the high
% fidelity

% Generate cheap points, given an upper bound and lower bound, change this
% for the different function
Xc = initialdesign(numlf,dim,-1,1);
% Expensive Points using the slection function from Surrogates Toolbox
% Version 3.0 by Felipe A. C. Viana
Xe = srgtsDOESubSample(Xc,numhf);


% Evaluate the expensive function change this to the function you wish to
% optimise
Ye = rastrigin(Xe,10000,1);

% evaluate the cheap function change this to the function that is the low
% fidelity data
Yc = rastrigin(Xc,1000,1);

% fit an initial Cokriging model
model = Cokriging(Xe,Ye,Xc,Yc);
% set paramters for the maxVarExpImp
lb = -ones(dim,1);
ub = ones(dim,1);
cost = 0;
i = 1;
nc = numlf - numhf;
n = numlf;
% main loop runs a set number of times.
while cost < totalcost
    % maximise the expected Improvement and find x associated with the maximum
    % expected improvement
    [x, fi] = maxVarExpImp(model, 1, lb, ub,hfidcost,lfidcost);
    % assuming the low fidelity cost is 1/10th the high fidelity cost
    % evaluate the low fidelity, change this to the desired low fidelity
    % function
    
    
    if cost + lfidcost >= totalcost, break; end
    yc = rastrigin(x,1000,1);
    cost = cost + lfidcost;
    % evaluate its high fidelity point if fi indicates high fidleity
    if fi == 0
        if sum(ismember(x,Xc,"rows")) == 0
            Yc = [Yc;yc];
            Xc = [Xc;x];
        else
            cost = cost - lfidcost;
        end
        if cost + hfidcost >= totalcost, break; end
        cost = cost + hfidcost;
        %Insuring that the cost of computation does not exceed that of the
        %highfidelity evaluation
        
        % Change this function to the desired function
        y = rastrigin(x,10000,1);
        Xe = [Xe;x];
        Ye = [Ye;y];
    else
        % store these values
        Yc = [Yc;yc];
        Xc = [Xc;x];
    end
    % fit a new model
    [model] = Cokriging(Xe,Ye,Xc,Yc);
    f(i) = fi;
    i = i+1;
    n = n+1;
    nc = nc+1;
end
[model] = Cokriging(Xe,Ye,Xc,Yc);
end

