function y = rastrigin(x,phi,error, d)
% The Rastrigin benchmark optimisation problem with error functions from A
% Generic Test Suite for Evolutionary Multi-fidelity Optimization. This
% versi
%
% Input:
%   x - the point you wish to calculate the rastrigin value for.
%   phi - the cost of the evaluaiation, a number between 0 and 10000 with
%   10000 being the objective function. In our tests we used 1000 as the
%   low fidelity function
%   error - either 1,2 or 3 this values determines the error function given
%   a resolution error. It uses the resolution error function 1,2 and 4.
%   d - The number of dimension of the array either 1 or 2
%
% Output
%   y - the function output given the error function and the true rastrigin
%   function value
%
if nargin < 4, d = 2; end
% In case no error function is used the default is the error function 1
if nargin < 3, error = 1; end


% True function output 
fe = sum(x.^2 + 1 - cos(4*pi.*x),d);

% Calculate the value of the error function based of the error input
% vaaraible
if error == 2
    theta = exp(-0.00025*phi)-exp(-2.5);
    error = sum(theta*cos(4*pi*theta.*x + 0.5*pi*theta + pi),d);
elseif error == 3
    theta = 1 - 0.0001*phi;
    error = sum(theta.*(1-abs(x)).*cos(4*pi*theta.*x + 0.5*pi*theta + pi),d);
else
    theta = 1 - 0.0001*phi;
    error = sum(theta*cos(4*pi*theta.*x + 0.5*pi*theta + pi),d);
end

% Calculate the output.
y = fe + error;
end