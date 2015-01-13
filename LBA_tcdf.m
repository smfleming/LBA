function F = LBA_tcdf(t, A, b, v, sv)
% Get CDF of first passage time of ith accumulator in LBA model
% F = LBA_tcdf(t, A, b, v, sv)
%
%
% SF 2012

g = (b-A-t.*v)./(t.*sv); % chizumax
h = (b-t.*v)./(t.*sv);  % chizu
i = b-t.*v; % chiminuszu
j = i-A; % xx

tmp1 = t.*sv.*(normpdf(g)-normpdf(h));
tmp2 = j.*normcdf(g) - i.*normcdf(h);

F = 1 + (tmp1 + tmp2)./A;