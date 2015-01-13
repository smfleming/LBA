function data = LBA_clean(data, cutoffs)
% Clean up RTs
% Excludes outliers based on cutoffs
% data = LBA_clean(data, [cutoffs])
%
% Default cutoffs are 200 and 2000 ms
% To change defaults specify cutoffs as optional argument (cutoffs = [LB UB])
%
%
% SF 2012

data = structfun(@(x) reshape(x,length(data.rt),1), data, 'UniformOutput', false);

if nargin == 1
    cutoffs = [200 2000];
end

exc = data.rt < cutoffs(1) | data.rt > cutoffs(2);

data = structfun(@(x) x(~exc,:), data, 'UniformOutput', false);