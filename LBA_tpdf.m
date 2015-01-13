function f = LBA_tpdf(t, A, b, v, sv)
% Get PDF of first passage time of ith accumulator in LBA model
% F = LBA_tpdf(t, A, b, v, sv)
%
%
% SF 2012

g = (b-A-t.*v)./(t.*sv);
h = (b-t.*v)./(t.*sv);

f = (-v.*normcdf(g) + sv.*normpdf(g) + v.*normcdf(h) - sv.*normpdf(h))./A;