% lba_wrapper
%
% script that nests lba model fitting and graphical output for individual subject data
%
% SF 2012

clear all
close all

%% simulate data
ntrials = 500;
% lba params
t0 = 425;
a = 380;
b = 600;
sv = 0.4;
v = [0.3 0.7];    % "2" is the correct response (higher drift rate)
n = 2;
for t = 1:ntrials
    [data.response(t) data.rt(t) modconf(t)] = LBA_trial(a, b, v, t0, sv, n);
end
data.stim = ones(1,ntrials)*2;    % correct response was "2"
data.cond = ones(1,ntrials);    % single difficulty condition

%% clean up data
data = LBA_clean(data);

%% get fitted parameters
% #### need to change to have model/parray as cell arrays for comparing
% models
model.v = 1;
model.A = 1;
model.b = 1;
model.sv = 1;
model.t0 = 1;
pArray = [0.8 300 150 0.4 200];

[params ll] = LBA_mle(data, model, pArray);

Ncond = max(data.cond);
cor = data.response == data.stim;
[v a b sv t0] = LBA_parse(model, params, Ncond);

nbin = 10;
minrt = 100;
maxrt = max(data.rt);
bins = linspace(minrt, maxrt, nbin);
bindiff = bins(3)-bins(2);

figure;
set(gcf,'position',[150 500 450.*Ncond 250]);
for j = 1:Ncond
    subplot(1,Ncond,j);
    
    % corrects
    histdata = hist(data.rt(cor & data.cond == j), bins);
    vi = [v(j) 1-v(j)];
    histpred = LBA_n1PDF(bins-t0(j),a(j),b(j)+a(j),vi,sv(j));
    % histpred is likelihood of a single trial at this rt
    % need to scale to area of bar (bindiff*trials)
    histpred = histpred.*bindiff.*sum(data.cond == j);
    h = bar(bins, histdata);
    set(h, 'edgecolor','b','linewidth',2,'facecolor','w');
    hold on
    plot(bins, histpred, 'b','linewidth',2);
    
    % errors
    histdata = hist(data.rt(~cor & data.cond == j), bins);
    vi = [1-v(j) v(j)];
    histpred = LBA_n1PDF(bins-t0(j),a(j),b(j)+a(j),vi,sv(j));
    % histpred is likelihood of a single trial at this rt
    % need to scale to area of bar (bindiff*trials)
    histpred = histpred.*bindiff.*sum(data.cond == j);
    h = bar(bins+(bindiff/2), histdata);
    set(h, 'edgecolor','r','linewidth',2,'facecolor','w');
    hold on
    plot(bins+(bindiff/2), histpred, 'r','linewidth',2);
    xlabel('rt');
    ylabel('freq');
    
end
