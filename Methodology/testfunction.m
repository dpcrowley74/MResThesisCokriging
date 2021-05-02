function y = testfunction(x, fid)
% The test function from Advances in Surrogate-Based Optimisation
%
%   Input:
%   x - the input variable and an element of [0,1]
%   fid - The level of fidelity, 1 indicates low fidelity.
%
%   Output:
%   y - the function output either high or lowe fidleity

% If only one input assume high fidelity
if nargin < 2, fid = 0; end

% High Fidelity function ouput
y = (6.*x-2).^2.*sin(12.*x-4);

if fid == 1
    % Low fidelity function output
    y = y./2 + 5.*(2.*x - 1)-5;
end
end