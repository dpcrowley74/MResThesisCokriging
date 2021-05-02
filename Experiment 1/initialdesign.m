function X = initialdesign(n,dimension, a, b)
% Generates the initial sample for the expensive and cheap functions from
% a latin hypercube whose lower bound is defined by a and upper bound
% defined by b.
%
% Inputs:
%       n - The size of the sample
%       dimension - the dimension of the design
%       a - the lower bound of the latin hypercube
%       b - the upper bound of the latin hypercube
%
% Outputs:
%       X - the design sample

if nargin < 3
    a = 1; 
    b = 4; 
end
% if no upper or  lower bounds are entered

X = (b-a)*lhsdesign(n,dimension) + a;