function [F] = Filter1D(dg_globals,Nc,s)

% function [F] = Filter1D(N,Nc,s)
% Purpose : Initialize 1D filter matrix of size N.
%           Order of exponential filter is (even) s with cutoff at Nc;
N=dg_globals.N;
r = JacobiGL(0,0,N);
V  = Vandermonde1D(N, r); invV = inv(V);


% Globals1D;
filterdiag = ones(N+1,1);
alpha = -log(eps);

% Initialize filter function
for i=Nc:N
    filterdiag(i+1) = exp(-alpha*((i-Nc)/(N-Nc))^s);
end;

F = V*diag(filterdiag)*invV;
return;