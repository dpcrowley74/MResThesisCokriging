function[x] = minModelPred(model,nvars,lb, ub)
% function for calculate the maximum variable expected improvement, with
% input the model and the parameters for the ga function of interest.
% outputs the location of the parametes and the fidelity of the parameters.
% Take Information out of the model
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
lower = ones(nvars,1)*lb;
upper = ones(nvars,1)*ub;
x = particleswarm(@pred, nvars, lower, upper);





    function yhat = pred(x)
        % Calculates the expected Improvement at x
        %
        % Inputs:
        %       x - 1xdim vector of design variables
        % Output:
        %       Imp- the expected Improvement at the point x(1:dim) with
        %       fidelity x(dim+1)
        % extract variables for data structure
        c = ones(n,1);
        one = ones(n,1);
        
        %fill psi vector
        if model.type == "Cokriging"
            for i = 1:nc
                c(i) = r*Sc*exp(-sum(T.*abs(X(i,:) -x).^pc));
            end
            for i = nc+1:n
                c(i) = ((r^2)*Sc+Sd)*exp(-sum(Te.*abs(X(i,:)-x).^pd));
            end
            yhat = mu + c'*(U\(U'\(Y -mu*one)));
        else
            for i = 1:n
                c(i) = exp(-sum(T.*abs(X(i,:) -x).^pc));
            end
            yhat = mu + c'*(Uc\(Uc'\(Y -mu*one)));
        end
    end
end

