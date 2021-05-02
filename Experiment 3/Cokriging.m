function [model] = Cokriging(Xe, Ye, Xc, Yc)
% Function for Building the Cokriging model.
%
% Inputs:
%       Xe - the predictor for the expensive function
%       Ye - the response for the expensive function
%       Xc - the predictor for the cheap function
%       Yc - the response for the cheap function
%       m - the number of attempts at finding maximum likelihood or a
%       positive definite matrix
%
% Outputs:
%       model - a model used for the predicition from the Cokriging
%       regressiosn


% Initialise loop
[n,dim] =size(Xe);
n2 = size(Xc,1);


%Estimate the cheap model parameters
lb = -3*ones(dim,1);
ub = 2*ones(dim,1);
lbd = [lb;0];
ubd = [ub;1-eps];
flag = 1;
i = 0;
xk0 = zeros(dim,1)';
xc0 = [zeros(dim,1);0.9]';
while flag ~= 0
    if i > 3
        lb = lb - 1;
        ub = ub + 1;
    end
    params = simulannealbnd(@likelihoodc, xk0 , lb, ub);
    Tc = params;
    pc = 1;
    % Once the kriging model is built we determine how many of the high
    % fidelity values have a matching low fidelity value
    if dim == 1
        [Xce, ie] = intersect(Xc,Xe);
    else
        [Xce, ie] = intersect(Xc,Xe,'rows');
    end
    Yce = Yc(ie);
    % We then sort the code so that the high-fidleity values are at the end
    % of the low fidelity values.
    if dim == 1
        [Xcc, ic] = setdiff(Xc,Xe);
    else
        [Xcc, ic] = setdiff(Xc,Xe,'rows');
    end
    Ycc = Yc(ic);
    % reassign Xc and Yc so that the high fidelity values are at the bottom
    Xc = [Xcc;Xe];
    Yc = [Ycc;Yce];
    %Estimate the Expensive log parameters
    
    params = simulannealbnd(@likelihoodd,xc0,lbd,ubd);
    Te = params(1:dim);
    pd = 1;
    rho = params(dim+1);
    %Build the Model
    [model,flag] = buildmodel(Yc, Xc, Tc, pc, Ye, Xe, rho, Te, pd);
    i = i + 1;
end

%-----------------------------------------------------------------------

    function[NegLnLike] = likelihoodd(x)
        % Calculates the negative of the concentrated ln-likelihood.
        %
        % Inputs:
        %       x - vector of log(parameters) to be maximised
        %       Xe - The input data for the high fidelity function
        %       Ye - The Output for the high fidelity function
        %       Yc - The value of the low fidelity function value associated with
        %       the data from Xe
        %
        % Outputs:
        %       NegLnLike - concentrated ln - likelihood * -1 for minimizing
        
        
        r = x(dim+1);
        [ne,k] = size(Xe);
        theta = 10.^x(1:k);
        p = 2;
        one = ones(ne,1);
        
        
        % Pre-allocate memory
        PsidXe = zeros(ne,ne);
        % Building the upper half of the correlation matrix.
        for i2 = 1:ne
            for j = i2+1:ne
                PsidXe(i2,j) = exp(-sum(theta.*abs(Xe(i2,:)-Xe(j,:)).^p));
            end
        end
        
        % Adding upper and lower halves and diagonal of ones plus
        % a small number to reduce ill conditioning
        PsidXe = PsidXe + PsidXe' + eye(ne) + eye(ne).*eps;
        
        % Cholesky factorization
        [U,f] = chol(PsidXe);
        
        % Penalty if ill-conditioned
        if f>0
            NegLnLike = 1e4;
        else
            
            % Sum lns of diagonal to find ln(det(psi))
            LnDetPsidXe = 2*sum(log(abs(diag(U))));
            
            % Difference vector
            d = Ye - r.*Yce;
            
            % Back-substitutio of Cholesky instead of inverse
            mud = (one'*(U\(U'\d)))/(one'*(U\(U'\one)));
            SigmaSqrd = (d - one*mud)'*(U\(U'\(d-one*mud)))/ne;
            NegLnLike = -1*(-(ne/2)*log(SigmaSqrd) - 0.5*LnDetPsidXe);
        end
    end

%--------------------------------------------------------------

    function[NegLnLike] = likelihoodc(x)
        % Calculates the negative of the concentrated ln-likelihood.
        %
        % Inputs:
        %       x - vector of log(parameters) to measure ln-likelihood
        %       X - the Ipredictor data
        %       Y - The response data
        %
        % Outputs:
        %       NegLnLike - concentrated ln - likelihood * -1 for minimizing
        
        
        theta = 10.^x(1:dim);
        p = 2;
        one = ones(n2,1);
        
        
        % Pre-allocate memory
        Psi = zeros(n2,n2);
        % Building the upper half of the correlation matrix.
        for i3 = 1:n2
            for j = i3+1:n2
                Psi(i3,j) = exp(-sum(theta.*abs(Xc(i3,:)-Xc(j,:)).^p));
            end
        end
        
        % Adding upper and lower halves and diagonal of ones plus
        % a small number to reduce ill conditioning
        Psi = Psi + Psi' + eye(n2) + eye(n2).*eps;
        
        % Cholesky factorization
        [U,f] = chol(Psi);
        
        % Penalty if ill-conditioned
        if f>0
            NegLnLike = 1e4;
        else
            
            % Sum lns of diagonal to find ln(det(psi))
            LnDetPsi = 2*sum(log(abs(diag(U))));
            
            % Back-substitutio of Cholesky instead of inverse
            mu = (one'*(U\(U'\Yc)))/(one'*(U\(U'\one)));
            SigmaSqr = (Yc - one*mu)'*(U\(U'\(Yc-one*mu)))/n2;
            NegLnLike = -1*(-(n2/2)*log(SigmaSqr) - 0.5*LnDetPsi);
        end
    end
end
