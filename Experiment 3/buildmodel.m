function [model,f] = buildmodel(Yc, Xc, Tc, pc, Ye, Xe, r, Te, pd)
% Function for building the Covariance matrix of the Cokriging model
%
% Inputs:
%       Yc - The sample from the cheap function
%       Xc - The values associated with the cheap evaluation
%       Tc - The estimate of the scalar parameters associated with the cheap function
%       pc - The estimate of the exponential parameters associated with the
%       cheap function
%       Ye - The sample from the expensive function
%       Xe - The values associated with the expensive evaluation
%       r - The estiamte associated with the regresssian parameter
%       Te - The estimate associated with the scalar parameters for
%       the expensive function
%       pd - The estimate associated with the exponential parameter for the
%       expensive function
%
%
% Outputs:
%       An object called model containing all the data for prediction

% Start building the model by defining what is already defined.
T = 10.^Tc;
p = 2.^pc;
model.Y = Yc;
model.X = Xc;
model.T = Tc;
model.p = pc;

% Defining constants
n = size(Xc,1);
onec = ones(n,1);

% Pre-allocate memory
Psi = zeros(n,n);
% Building the upper half of the correlation matrix.
for i = 1:n
    for j = i+1:n
        Psi(i,j) = exp(-sum(T.*abs(Xc(i,:)-Xc(j,:)).^p));
    end
end

% Adding upper and lower halves and diagonal of ones plus
% a small number to reduce ill conditioning
Psi = Psi + Psi' + eye(n) + eye(n).*eps;

% Cholesky factorization
[Uc,f] = chol(Psi);
if f == 0
    % Back-substitutio of Cholesky instead of inverse
    muc = (onec'*(Uc\(Uc'\Yc)))/(onec'*(Uc\(Uc'\onec)));
    SigmaSqrc = (Yc - onec*muc)'*(Uc\(Uc'\(Yc-onec*muc)))/n;
    
    % the data required for evaluation.
    model.Uc = Uc;
    model.Sc = SigmaSqrc;
    model.mu = muc;
    model.type = "Kriging";
    if nargin > 4
        %define varaibles for the expensive function for the model
        model.Xe = Xe;
        model.Yc  = Yc;
        model.Ye = Ye;
        model.Te = Te;
        model.pd = pd;
        model.r = r;
        
        rho = r;
        T = 10.^Te;
        p = 2.^pd;
        
        ne = size(Xe,1);
        nc = n - ne;
        oned = ones(ne,1);
        
        % Pre-allocate memory
        PsidXe = zeros(ne,ne);
        % Building the upper half of the correlation matrix.
        for i = 1:ne
            for j = i+1:ne
                PsidXe(i,j) = exp(-sum(T.*abs(Xe(i,:)-Xe(j,:)).^p));
            end
        end
        
        % Adding upper and lower halves and diagonal of ones plus
        % a small number to reduce ill conditioning
        PsidXe = PsidXe + PsidXe' + eye(ne) + eye(ne).*eps;
        
        % Cholesky factorization
        [Ud,f] = chol(PsidXe);
        if f == 0
            % Difference vector
            d = Ye - rho.*Yc(nc+1:n);
            
            % Back-substitutio of Cholesky instead of inverse
            mud = (oned'*(Ud\(Ud'\d)))/(oned'*(Ud\(Ud'\oned)));
            SigmaSqrd = (d - oned*mud)'*(Ud\(Ud'\(d-oned*mud)))/ne;
            
            % Evaluate Covaraince Matrix for Cokriging
            C = [SigmaSqrc*Psi(1:nc,1:nc), rho * SigmaSqrc*Psi(1:nc,nc+1:n);
                rho * SigmaSqrc*Psi(nc+1:n,1:nc), rho^2*SigmaSqrc*Psi(nc+1:n,nc+1:n) + SigmaSqrd*PsidXe];
            % Cholesky factorization
            [U,f] = chol(C);
            if f== 0
                [U1,f] = chol(rho^2*SigmaSqrc*Psi(nc+1:n,nc+1:n) + SigmaSqrd*PsidXe);
            end
            if f == 0
                
                %Calculate mu
                mu = (oned'*(U1\(U1'\Ye)))/(oned'*(U1\(U1'\oned)));
                model.Y = [Yc(1:nc);Ye];
                model.mu = mu;
                model.Sd = SigmaSqrd;
                model.U = U;
                model.type = "Cokriging";
            end
        end
    end
end

