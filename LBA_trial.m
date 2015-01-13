function [choice RT conf] = LBA_trial(A, b, v, t0, sv, N)
% Run a single trial of the LBA model (Brown & Heathcote, 2008, Cog
% Psychol)
%
% Usage: [choice RT conf] = LBA_trial(A, b, v, t0, sv, N)
%
% Inputs:
%
% A = range of uniform distribution U[0,A] from which starting point k is
% drawn
% b = bound
% v = vector of drift rates
% sv = standard deviation of drift rate
% t0 = non-decision time
% N = number of response options
%
% Outputs:
%
% choice = scalar from 1:N indicating response chosen by model
% RT = reaction time in ms
% confidence = confidence computed using balance of evidence rule (Vickers,
% 1979)
%
% SF 2012

trialOK = false;

while ~trialOK
    for i = 1:N
        
        % Get starting point
        k(i) = rand.*A;
        
        % Get drift rate
        d(i) = normrnd(v(i), sv);
        
        % Get time to threshold
        t(i) = (b-k(i))./d(i);
        
        % Add on non-decision time
        allRT(i) = t0 + t(i);
    end
    
    % Get choice and confidence
    [RT choice] = min(allRT);
    
    % Confidence is equal to threshold minus value of next best accumulator at decision
    % time
    j=1;
    if N == 1
        conf = NaN;
    else
        for i = 1:N
            if i ~= choice
                z(j) = t(choice).*d(i) + k(i);
                j=j+1;
            end
        end
        [nb i] = max(z);
        conf = b-nb;
    end
    
    % Check we have not sampled negative drift(s)
    if RT > 0 & nb > 0
        trialOK = true;
    end
end