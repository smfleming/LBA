function f = LBA_wpdf(z, t, A, v, sv)
% Get PDF of activation values z of ith accumulator at time t in LBA model
% w/o boundaries
% f = LBA_wpdf(z, t, A, v, sv)
%
% See LBA_mle for definition of parameters
%
% SF 2012

g = normcdf((z-A)./t, v, sv);
h = normcdf(z./t, v, sv);

f = (h-g)./A;