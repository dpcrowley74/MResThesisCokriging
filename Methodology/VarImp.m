function Imp = VarImp(x, model, fid, hificost, lowficost)
% Calculates the expected Improvement at x
%
% Inputs:
%       x - 1xdim vector of design variables
% Output:
%       Imp- the expected Improvement at the point x(1:dim) with
%       fidelity x(dim+1)
% extract variables for data structure
X = model.X;
n = size(X,1);
Y = model.Y;
T = 10.^model.T;
Uc = model.Uc;
Sc = model.Sc;
pc = 2.^model.p;
mu = model.mu;
r = 1;
if model.type == "Cokriging"
    Ye = model.Ye;
    Te = 10.^model.Te;
    r = model.r;
    Sd = model.Sd;
    pd = 2.^model.pd;
    ne = size(Ye,1);
    nc = n - ne;
    U = model.U;
end
c = zeros(n,1);
one = ones(n,1);

%fill psi vector
if model.type == "Cokriging"
    for i = 1:nc
        c(i) = r*Sc*exp(-sum(T.*abs(X(i,:) -x).^pc));
    end
    for i = nc+1:n
        for i = nc+1:n
            c(i) = (r^2)*Sc*exp(-sum(T.*abs(X(i,:)-x).^pc))+Sd*exp(-sum(Te.*abs(X(i,:)-x).^pd));
        end
    end
    yhat = mu + c'*(U\(U'\(Y -mu*one)));
    if fid == 0
        S = (r^2*Sc + Sd - c'*(U\(U'\c)));
        cost = hificost;
    end
end
if fid == 1 || model.type == "Kriging"
    for i = 1:n
        c(i) = exp(-sum(T.*abs(X(i,:) -x).^pc));
    end
    if model.type == "Kriging"
        yhat = mu + c'*(Uc\(Uc'\(Y -mu*one)));
    end
    S = (r^2*Sc*(1 - c'*(Uc\(Uc'\c))));
    cost = lowficost;
end
if S <= 0.01
    Imp = cost;
else
    SSqr = sqrt(S);
    ymin = min(Ye);
    c = normcdf((ymin-yhat)./SSqr,0,1);
    p = normpdf((ymin-yhat)./SSqr,0,1);
    Imp = ((ymin-yhat).*c + SSqr.*p) - cost;
end
end