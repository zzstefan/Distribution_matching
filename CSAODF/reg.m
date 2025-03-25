% The thresholding regularization for CSA-ODF computation, with the
% parameter delta. 
% 
% x = reg(x, delta)
%
% See also:  reconCSAODF, reconCSAODF3Q, CSAODF_CLI, EXAMPLE, EXAMPLE_CLI.

% Code by Iman Aganj.

function x = reg(x, delta)

if delta==0
    x = (x>=0 & x<1).*x + (x>=1);
elseif delta > 0
    x = (x<0).*(.5*delta) + (x>=0 & x<delta).*(.5*delta + (x.^2)*(.5/delta)) + (x>=delta & x<1-delta).*x + (x>=1-delta & x<1).*(1-.5*delta-((1-x).^2)*(.5/delta)) + (x>=1).*(1-.5*delta);
else
    error('''delta'' cannot be negative!')
end
