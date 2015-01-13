function out = LBA_plot(data, params, model)
% LBA_plot(data, params)
%
% Produces plots of LBA model predictions and data to assess fit
%
% Can handle both fits by cond or coherence (see LBA_mle for details)
%
% Currently in sandbox form, needs tidying up
% SF 2012

Ncond = max(data.cond);
[v A b sv t0] = LBA_parse(model, params, Ncond);
rtQ = [0 0.1 0.3 0.5 0.7 0.9 1.0];
out = [];
Nboot = 100;
% Check what data we have
if isfield(data, 'coherence')
    %% 2-choice coherence design
    if max(data.response) == 2  % 2AFC, can plot psychometric function
        
        dv = abs(diff(data.coherence'))';
        
        % Simulate trials to get mean RT and confidence for each dv
        % combination (get accuracy analytically)
        for t = 1:length(data.coherence)
            if data.correct(t)
                vi = [max(data.coherence(t,:)) min(data.coherence(t,:))];
            else
                vi = [min(data.coherence(t,:)) max(data.coherence(t,:))];
            end
            for i = 1:Nboot
                [choice RT(t,i) conf(t,i)] = LBA_trial(A, A+b, v.*vi, t0, sv, 2);
            end
        end
        pred.conf = mean(conf');
        pred.rt = mean(RT');
        
        % Group coherence levels together for each response, get mean and
        % predicted RT
        bins = 5;
        cohQ = quantile(data.coherence(:),linspace(0,1,bins));
        for i = 1:length(cohQ)-1
            for j = 1:length(cohQ)-1
                
                % Observed
                meanRT(i,j) = mean(data.rt(data.correct == 1 & data.coherence(:,1) > cohQ(i) &  data.coherence(:,1) <= cohQ(i+1)...
                    & data.coherence(:,2) > cohQ(j) &  data.coherence(:,2) <= cohQ(j+1)));
                meanAcc(i,j) = mean(data.correct(data.coherence(:,1) > cohQ(i) &  data.coherence(:,1) <= cohQ(i+1)...
                    & data.coherence(:,2) > cohQ(j) &  data.coherence(:,2) <= cohQ(j+1)));
                
                % Predicted
                predRT(i,j) = mean(pred.rt(data.correct == 1 & data.coherence(:,1) > cohQ(i) &  data.coherence(:,1) <= cohQ(i+1)...
                    & data.coherence(:,2) > cohQ(j) &  data.coherence(:,2) <= cohQ(j+1)));
                currStim = [mean([cohQ(i) cohQ(i+1)]) mean([cohQ(j) cohQ(j+1)])];
                vi = [max(currStim) min(currStim)]; % Put highest drift rate in first place (for RT on correct trials)
                predAcc(i,j) = LBA_n1CDF(A,A+b,vi.*v,sv);
                
            end
        end
        
        figure;
        k=1;
        for i = 1:length(cohQ)-1
            subplot(2,bins-1,k);
            bar(meanRT(i,:),'LineWidth',2,'facecolor','w')
            hold on
            plot(predRT(i,:),'o ','MarkerSize',7);
            k=k+1;
        end
        for i = 1:length(cohQ)-1
            subplot(2,bins-1,k);
            bar(meanAcc(i,:),'LineWidth',2,'facecolor','w')
            hold on
            plot(predAcc(i,:),'o ','MarkerSize',7);
            k=k+1;
        end
        
        % RT-DV-conf heatmaps
        baseStim = 0.1;
        bins = 5;
        surfDV = quantile(dv,linspace(0,1,bins));
        surfRT = quantile(data.rt,linspace(0,1,bins));
        u_v = max(surfDV)-surfDV + baseStim;  % Derive unchosen integrator from DV values
        for i = 1:length(surfDV)-1
            for j = 1:length(surfRT)-1
                surfConf(i,j) = mean(data.conf((data.rt>surfRT(j) & data.rt<=surfRT(j+1)) & (dv>surfDV(i) & dv<=surfDV(i+1))));
                modelConf(i,j) = LBA_conf([surfRT(j) surfRT(j+1)]-t0, A, A+b, u_v(i).*v, sv);
            end
        end
        figure;
        subplot(1,2,1);
        imagesc(surfRT,surfDV,surfConf);
        axis square
        colormap gray
        title('Data');
        subplot(1,2,2);
        imagesc(surfRT,surfDV,modelConf);
        axis square
        colormap gray
        title('Model')
        
        % Package output
        out.meanAcc = meanAcc;
        out.predAcc = predAcc;
        out.meanRT = meanRT;
        out.predRT = predRT;
        out.meanConf = surfConf;
        out.predConf = modelConf;
        
    end
else
    %% 2-choice by condition design
    cor = data.response == data.stim;
    nbin = 8;
    bins = linspace(min(data.rt), max(data.rt), nbin);
    bindiff = bins(3)-bins(2);
    
    % RT distributions (analytic)
    figure;
    set(gcf,'Position',[150 500 450.*Ncond 250]);
    for j = 1:Ncond
        subplot(1,Ncond,j);
        
        % corrects
        histdata = hist(data.rt(cor & data.cond == j), bins);
        vi = [v(j) 1-v(j)]; % this assumes two choice, and unchosen response as 1-v
        histpred = LBA_n1PDF(bins-t0(j),A(j),b(j)+A(j),vi,sv(j));
        % histpred is likelihood of a single trial at this RT
        % need to scale to area of bar (bindiff*trials)
        histpred = histpred.*bindiff.*sum(data.cond == j);
        h = bar(bins, histdata);
        set(h, 'edgecolor','b','LineWidth',2,'facecolor','w');
        hold on
        plot(bins, histpred, 'b','LineWidth',2);
        
        % errors
        histdata = hist(data.rt(~cor & data.cond == j), bins);
        vi = [1-v(j) v(j)];
        histpred = LBA_n1PDF(bins-t0(j),A(j),b(j)+A(j),vi,sv(j));
        % histpred is likelihood of a single trial at this RT
        % need to scale to area of bar (bindiff*trials)
        histpred = histpred.*bindiff.*sum(data.cond == j);
        h = bar(bins+(bindiff/2), histdata);
        set(h, 'edgecolor','r','LineWidth',2,'facecolor','w');
        hold on
        plot(bins+(bindiff/2), histpred, 'r','LineWidth',2);
        xlabel('RT');
        ylabel('Freq');
        
    end
    
    
    % RT - conf relationship ### currently can't deal with more than one
    % condition
    if isfield(data, 'conf')
        figure;
        subplot(1,Ncond,j);
        
        % corrects
        vi = [v(j) 1-v(j)];
        for i = 1:length(bins)-1
            meanConfCor(i) = mean(data.conf(cor & data.cond == j & data.rt > bins(i) & data.rt <= bins(i+1)));
            [predConfCor(i) SD] = LBA_bootConf(bins(i)-t0(j), A(j), A(j)+b(j), vi(2), sv(j), 1000);
            x(i) = mean([bins(i+1) bins(i)]);
        end
        
        % Errors
        vi = [1-v(j) v(j)];
        for i = 1:length(bins)-1
            meanConfErr(i) = mean(data.conf(~cor & data.cond == j & data.rt > bins(i) & data.rt <= bins(i+1)));
            [predConfErr(i) SD] = LBA_bootConf(bins(i)-t0(j), A(j), A(j)+b(j), vi(2), sv(j), 1000);
        end
        
        zConf = zscore([predConfCor predConfErr]);
        
        plot(x, meanConfCor, 'bo ','MarkerSize',7);
        hold on
        plot(x, zConf(1:length(predConfCor)), 'b','LineWidth',2);
        plot(x, meanConfErr, 'ro ','MarkerSize',7);
        plot(x, zConf(length(predConfCor)+1:end), 'r','LineWidth',2);
        
        % Package output
        out.meanConf = [meanConfCor; meanConfErr];
        out.predConf = [zConf(1:length(predConfCor)); zConf(length(predConfCor)+1:end)];
    end
    
    
end