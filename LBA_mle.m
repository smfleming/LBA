function [params LL] = LBA_mle(data, model, pArray)
% Continuous maximum likelihood estimation of LBA model
% [params fVal] = LBA_mle(data, model, pArray)
%
% Provides maximum likelihood fits to vector of choice and RT data. User
% can set various parameterisations of model via "model" structure
%
% Usage:
%
% Inputs:
%
% data --> structure containing following fields, all vectors of length
% 1 x trials
% (RT data and missed trials should be cleaned before passing to LBA_mle, see LBA_clean)
% data.rt - RT in milliseconds
% data.cond - condition vector (e.g. 1=easy, 2=hard)
% if no conditions, specify all as "1"
% data.stim - stimulus code (e.g. 1=left, 2=right)
% data.response - response code (e.g. 1=left, 2=right)
% data.coherence - optional, if present code will estimate link
% parameter from stimulus coherence into drift rate (cf. Palmer et al.,
% 2005, J Vis)
% All fields in data should be columnar vectors/matrices
%
% model --> structure containing information on which parameters to
% share between conditions. Must contain fields for v, A, b, sv, t0. Each
% field must be a scalar equal to either 1 or Ncond.
% E.g. to share bounds between 3 conditions, but to keep drift rates
% constant, set:
% model.v = 1; model.A = 1; model.b = 3; model.sv = 1; model.t0 = 1;
%
% pArray --> vector of starting points for parameters to be estimated
% pArray = [v A b-A t0 sv]
% length of pArray corresponds to setup in model
%
% Outputs:
%
% params --> vector of fitted parameters in same order as pArray (note
% output of b-A, need to subtract condition-specific A to get b)
%
% LL --> log-likelihood of data given model
%
% SF 2012 sf102@nyu.edu

options = optimset('Display','iter','MaxFunEvals',100000);
LB = ones(1,length(pArray)).*1e-5;
UB = ones(1,length(pArray)).*Inf;
[params fVal] = fmincon(@fitfunc,pArray,[],[],[],[],LB,UB,[],options);

LL = -fVal;

    function negLL = fitfunc(pArray)
        
        Ncond = max(data.cond);
        Nresp = max(data.response);
        ntrials = length(data.response);
        
        % ensure data is in right format
        data = structfun(@(x) reshape(x,ntrials,1), data, 'UniformOutput', false);
        
        [v A b sv t0] = LBA_parse(model, pArray, Ncond);
        A = real(log(A));
        b = real(log(b));
        t0 = real(log(t0));
        
        %% Get likelihoods
        if isfield(data, 'coherence')
            vi = data.coherence;
            for t = 1:length(vi)  % which drift rate to place in 1st slot (based on subject's response) - faster way to do this??
                currFirst = vi(t,1);
                vi(t,1) = vi(t,data.response(t));
                vi(t,data.response(t)) = currFirst;
            end
            vi = repmat(v(data.cond),1,Nresp).*vi; % linear mapping into drift rate
            rtfit = data.rt - exp(t0(data.cond));
            % trial likelihoods
            p = LBA_n1PDF(rtfit, exp(A(data.cond)), exp(b(data.cond)) + exp(A(data.cond)), vi, sv(data.cond));
            
        elseif isfield(data, 'cond')
            
            % Get log-liks for these parameters
            cor = data.response == data.stim;
            if Ncond == 1
                vi = [repmat(v,length(data.cond),1) repmat(1-v,length(data.cond),1)];
            else
                vi = [v(data.cond) 1-v(data.cond)];
            end
            vi(logical(~cor),:) = vi(logical(~cor),[2 1]);    % flip if incorrect
            rtfit = data.rt - exp(t0(data.cond));
            
            % trial likelihoods
            p = LBA_n1PDF(rtfit, exp(A(data.cond)), exp(b(data.cond)) + exp(A(data.cond)), vi, sv(data.cond));
            
        else
            fprintf('\n\n\nBad input! See help LBA_mle.\n\n');
            return;
        end
        
        p(p<=1e-5) = 1e-5;  % avoid underflow
        negLL = -sum(log(p));
    end
end