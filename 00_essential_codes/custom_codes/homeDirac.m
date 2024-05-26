function [y]= homeDirac(n,x)
if nargin==1
    y = dirac(x);
    idx = y == Inf; % find Inf
    y(idx) = 1; 
    idx = y == -Inf; % find Inf
    y(idx) = -1; 
else
    y = dirac(n,x);
    idx = y == Inf; % find Inf
    y(idx) = 1; 
    idx = y == -Inf; % find Inf
    y(idx) = -1; 
end
return