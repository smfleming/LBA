function [conf SDconf] = LBA_bootConf(t, A, b, u_v, sv, Nboot)
% conf = LBA_bootConf(t, A, b, u_v, sv, Nboot)
%
% Bootstrap predictions for the unchosen integrator(s) from several draws from
% posterior over z (random draws of k and d) at time t
%
% Nboot is number of samples 
% u_v is a vector of unchosen drift rates
%
% sampled confidence is derived from the balance of evidence between the
% chosen and next-best unchosen integrator, normalised between 0.5 and 1
% conf is mean of samples, SDconf is SD
%
% SF 2012

ndrift = length(u_v);
for i = 1:Nboot
    sampOK = false;
    while ~sampOK
        k = rand(1,ndrift).*A;
        d = normrnd(u_v, sv);
        z = k + t.*d;
        samp(i) = b./(b + max(z));
        if max(z) < b
            sampOK = true;
        end
    end
end

conf = mean(samp);
SDconf = std(samp);