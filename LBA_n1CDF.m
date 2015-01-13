function cdf = LBA_n1CDF(A, b, v, sv)
% Generates choice probability for responses on node #1 by numerical
% integration of LBA_n1PDF
%
% pdf = LBA_n1CDF(A, b, v, sv)
%
% SF 2012

cdf = quad(@(t)LBA_n1PDF(t,A,A+b,v,sv),1,100000);