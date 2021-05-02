function [yhat] = pred(x, model, type)
% Calculates the kriging prediction at x and estimates the mean square
% error
%
% Inputs:
%       x - 1xk vector of design variables
%       model - a data structure called model
%       type - determine whether Kriging or Cokriging model is fitted
%
% Output:
%       yhat- the expected value given by the model
%       s - the variance given by the model
% extract variables for data structure
if nargin < 3, type = 0; end
X = model.X;
Y = model.Y;
T = 10.^model.T;
Uc = model.Uc;
Sc = model.Sc;
pc = 2.^model.p;
mu = model.mu;
Yc = model.Yc;
n = size(X,1);
c = ones(n,1);
one = ones(n,1);
if type == 0
    Ye = model.Ye;
    Te = 10.^model.Te;
    r = model.r;
    Sd = model.Sd;
    pd = 2.^model.pd;
    ne = size(Ye,1);
    nc = n - ne;
    U = model.U;
end


%fill psi vector by model type
if type ==0
    for i = 1:nc
        c(i) = r*Sc*exp(-sum(T.*abs(X(i,:) -x).^pc));
    end
    for i = nc+1:n
        c(i) = (r^2)*Sc*exp(-sum(T.*abs(X(i,:) -x).^pc))+Sd*exp(-sum(Te.*abs(X(i,:)-x).^pd));
    end
    yhat = mu + c'*(U\(U'\(Y -mu*one)));
else
    for i = 1:n
        c(i) = exp(-sum(T.*abs(X(i,:) -x).^pc));
    end
    yhat = mu + c'*(Uc\(Uc'\(Yc -mu*one)));
end
end
