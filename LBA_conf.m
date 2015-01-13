function c = LBA_conf(t, A, b, u_v, sv)
% Get maximum likelihood estimate for confidence for response i at time t
% given set of LBA parameters
%
% function c = LBA_conf(t, A, b, v, sv)
%
% Model calculates balance of evidence between chosen and unchosen response
% as proxy for confidence (Vickers, 1979). Relies on estimating
% distribution of activation on unchosen integrator after
% chosen response is made
%
% Inputs:
%
% See LBA_mle for A, b, sv
%
% u_v = unchosen drift rate
%
% Maximum likelihood estimate is taken as reported confidence, but you
% could also hack this to work with the full distribution
%
% If t is a 2-element vector, returns confidence numerically integrated over
% this time window
%
% See LBA_mle for definition of parameters
%
% SF 2012

z = linspace(1,b,200);  % activation range over which to estimate pdf
dt = 1;     % millisecond bins over which to calculate integral

if length(t) == 1
    
    w = LBA_wpdf(z,t,A,u_v,sv);
    confpdf = w;
elseif length(t) == 2
    base = t(1):dt:t(2);
    
    for i = 1:length(base)
        w = LBA_wpdf(z,base(i),A,u_v,sv);
        confpdf(:,i) = w;
    end
    confpdf = mean(confpdf,2); % marginalise over time
else
    fprintf('\n\n\nBad input! See help LBA_conf.\n\n');
    return
end

try
    c = z(confpdf == max(confpdf));
catch
    c = median(z(confpdf == max(confpdf)));
end
c = b./(b + c);