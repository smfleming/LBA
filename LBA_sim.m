% Script to generate heatmap as in Nat Neuro paper but here using LBA model
%
% SF 2012

clear all

k = 1;  % parameter that maps stimulus values into drift rates
lum = rand(2,1000).*2; % two stimulus values

% LBA params
t0 = 425;
A = 380;
b = 600;
sv = 0.4;
N = 2;

for t = 1:length(lum)
    
    v = k.*lum(:,t)';
    [choice(t) RT(t) conf(t)] = LBA_trial(A, b, v, sv, t0, N);
    
end

% make heatmap
bins = 6;
dv = diff(lum);

surfDV = quantile(abs(dv),linspace(0,1,bins));
surfRT = quantile(RT,linspace(0,1,bins));
[X1, X2] = meshgrid(surfRT(1:end-1), surfDV);
for i = 1:length(surfDV)-1
    for j = 1:length(surfRT)-1
        surfConf(i,j) = mean(conf((RT>surfRT(j) & RT<=surfRT(j+1)) & (abs(dv)>surfDV(i) & abs(dv)<=surfDV(i+1))));
    end
end
imagesc(surfRT,surfDV,surfConf);
axis square