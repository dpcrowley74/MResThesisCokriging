function[model,benefit, iter, i] = minvarRastd(seed, numhf, numlf, dim, totalcost, hfidcost, lfidcost)
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
if nargin < 6, hfidcost = 0.01; end
if nargin < 7, lfidcost = 0.1*hfidcost; end
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
iter = 0;
benefit = min(Ye);
% main loop runs a set number of times.
while cost < totalcost
    % maximise the expected Improvement and find x associated with the maximum
    % expected improvement
    [x, fi, Imp] = maxVarExpImp(model, 1, lb, ub, hfidcost, lfidcost);
    if Imp >= 0 && iter == 0
        iter = i+1;
    end
    % evaluate its high fidelity point if fi indicates high fidleity
    if fi == 0
        cost = cost + hfidcost;
        %Insuring that the cost of computation does not exceed that of the
        %highfidelity evaluation
        if cost > totalcost, break; end
        % Change this function to the desired function
        y = rastrigin(x,10000,1);
        if y < min(Ye)
            benefit = [benefit, benefit(i) + (min(Ye) - y) - (hfidcost - lfidcost)];
        else
            benefit = [benefit, benefit(i) - (hfidcost - lfidcost)];
        end
        Xe = [Xe;x];
        Ye = [Ye;y];
        f(i) = fi;
    end
    cost = cost + lfidcost;
    if cost > totalcost, break; end
    f(i) = fi;
    % evaluate the low fidelity, change this to the desired low fidelity
    % function
    yc = rastrigin(x,1000,1);
    if fi == 1
        benefit = [benefit, benefit(i) - lfidcost];
    end
    % store these values
    if ismember(x,Xc,'Rows') == 0
        Yc = [Yc;yc];
        Xc = [Xc;x];
    else
        cost = cost - lfidcost;
    end
    % fit a new model
    f(i) = fi;
    model = Cokriging(Xe,Ye,Xc,Yc);
    i = i+1;
end
end

