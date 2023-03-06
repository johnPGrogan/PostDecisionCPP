% CPPPeakLatency
% post-CPP peak surface plot
% e.g. Twomey et al., 2016 Fig 2E
% low-pass filter data (6Hz, zero phase-shift fourth-order Butterworth IIR filter
% then find peak latency within window (200-1600ms)
% smooth over bins of 100 trials with Gaussian-weighted moving average
% then plot heatmap of time vs amplitude, sorted by peak latency, with latency superimposed

clc; clear; close all;
if ~exist('filter_tf','file'); eeglab nogui;  end

%% load

load('./Saves/CPPAnalysis_CSD_.mat','cpp','behData','eeg','fileInfo','labels');

[nPP, nT, nTr] = size(cpp);

colours.certainty = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];
colours.cond = [0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
colours.CoM = [0 0 1; 1 0 0; .2 .2 .2]; %[0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2]; % CoM
colours.conf3 = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];%colours.certainty;
colours.stimDur = tail(crameri('lajolla',5),4);
colours.acc = [1 0 1; 0 .8 0; .2 .2 .2];
colours.confInR1 =  [flipud(crameri('roma',6)); .2 .2 .2];

names.confInR1 = 'confidence-in-initial-choice';
names.certainty = 'confidence-in-final-choice';
names.conf3 = 'confidence-in-initial-choice: binned';
names.CoM = 'change-of-mind';
names.acc = 'initial accuracy';
%% filter

loPass = 6; % Hz

isAllNaN = repmat(all(isnan(cpp),2),1,size(cpp,2),1);
cpp(isAllNaN) = 0;   

% filter
cppFilt = eegfilt(reshape(permute(cpp,[1,3,2]),fileInfo.nPP*fileInfo.maxTr,size(cpp,2)),eeg.fs,0, loPass)'; % low-pass
cppFilt = permute(reshape(cppFilt,size(cpp,2),fileInfo.nPP,fileInfo.maxTr),[2,1,3]);

% set zeros back to NAN
cppFilt(isAllNaN) = NaN;
cpp(isAllNaN) = NaN;

%% find latency of max ampl in window

peakWindow = [0 1000]; % ms

peakWindowInds = isBetween(eeg.respTimes, peakWindow);

[peakAmpls, peakInds] = max(cppFilt(:,peakWindowInds,:),[],2);
peakAmpls = sq(peakAmpls);


peakLat = sq(eeg.respTimes(peakInds + find(peakWindowInds,1,'first')-1)); % ms
peakLat(isnan(peakAmpls)) = NaN;
peakInds(isnan(peakAmpls)) = NaN;

% if peakLat is at upper limit, not really a peak is it?
% so just exclude it?

%% heatmap sorted
% x=time, y=single trials (sorted by peak latency)

% sort by peak lat, across pps 
cppSorted = reshape(permute(cppFilt, [1 3 2]), nPP*nTr, nT); %[pp*tr, T]
[~, latOrder] = sort(peakLat(:));

cppSorted = cppSorted(latOrder,:);
peakLatSorted = peakInds(latOrder);


% remove NaN
toRemove = all(isnan(cppSorted),2);% | peakLatSorted > 512 | peakLatSorted<=1; %exclue those at edges

cppSorted(toRemove,:) = [];
peakLatSorted(toRemove) = [];


% smooth over trials?
n=100;
cppSortedSmooth = movmean(cppSorted, n, 1);
peakLatSortedSmooth = movmean(peakLatSorted, n, 1);


plotWin = [-200 1000];%peakWindow;
xInds = isBetween(eeg.respTimes,plotWin);
x = eeg.respTimes(xInds);
[~, xTickInds] = min(abs(x - [-200:200:max(x)]'),[],2);

figure();
colormap('jet');%(crameri('bam'));
imagesc(cppSortedSmooth(:,xInds),[-150 150]);
set(gca,'YDir','normal');
colorbar;
xlabel('time from init RT (ms)')
xticks(xTickInds); xticklabels(round(x(xTickInds),1,'significant'));
xline(find(x>0,1,'first'),':');
ylabel('trials');
yticks(5000:5000:length(peakLatSorted));

hold on;
plot(peakLatSortedSmooth + round((peakWindow(1)-plotWin(1))*.512), 1:numel(peakLatSorted), '-k');

title('post-choice CPP peak latency');

%% split by cond
figure();

cond = col(behData.cond);
cond(toRemove)=[];

plotWin = [-200 1000];%peakWindow;
xInds = isBetween(eeg.respTimes,plotWin);
x = eeg.respTimes(xInds);
[~, xTickInds] = min(abs(x - [-200:200:max(x)]'),[],2);


for iC = 2
    subplot(2,2,1);
    cppSortedCondSmooth = movmean(cppSorted(cond==iC,:), n, 1);
    peakLatSortedCondSmooth = movmean(peakLatSorted(cond==iC,:), n, 1);


    colormap('jet');%(crameri('bam'));
    imagesc(cppSortedCondSmooth(:,xInds),[-150 150]);
    set(gca,'YDir','normal');
    colorbar;
    xlabel('time from init RT (ms)')
    xticks(xTickInds); xticklabels(round(x(xTickInds),1,'significant'));
    xline(find(x>0,1,'first'),':');
    ylabel('trials');
    yticks(1000:1000:length(peakLatSortedCondSmooth));
    
    hold on;
    plot(peakLatSortedCondSmooth + round((peakWindow(1)-plotWin(1))*.512), 1:numel(peakLatSortedCondSmooth), '-k');
    
    title(labels.cond(iC));
end


%% plot peak latency by cond + confInR1

factors = {'confInR1','conf3','certainty'};
nF = length(factors);
% figure(); t=tiledlayout(nF,1); ax=[];
for iF = 1
    factor = factors{iF};
%     ax(iF) = nexttile(t);
    subplot(2,2,2);
        
    set(gca,'ColorOrder',colours.(factor),'nextplot','replacechildren');
    peakLatByCondFac = reshape(groupMeans(peakLat,2,behData.cond + behData.(factor)*10),30,2,[]); %[pp cond fac]
    h = errorBarPlot(peakLatByCondFac,'type','bar');
    xticks(1:2); xticklabels(labels.cond);
    ylabel('post-choice peak CPP latency (ms)');
    ylim(minMax([h.YData]) .* [.9 1.1])
    legend(labels.(factor),'Location','Best');
    title(names.(factor));
end

%% plot conf-RT split by confInR1 + cond

factor = 'confInR1';

behDataByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),behData,'UniformOutput',0); %[pp cond tr]
confRTByCondFac = groupMeans(behDataByCond.confRT,3,behDataByCond.(factor)); %[pp cond conf]

subplot(2,2,3);
% set(gca,'ColorOrder',colours.cond,'nextplot','replacechildren');
% h = errorBarPlot(permute(confRTByCondFac,[1 3 2]),'type','bar','plotargs',{'LineWidth',2});
% xticks(1:6); xticklabels(labels.(factor));
set(gca,'ColorOrder',colours.(factor),'nextplot','replacechildren');
h = errorBarPlot(confRTByCondFac,'type','bar');
xticks(1:2); xticklabels(labels.cond);
ylabel('final-response RT (ms)');
legend(h, labels.(factor), 'Location','Best');

%% plot beta lat here 

% in BetaLatStuff, using endLat (time of n.s. difference from -1, only for
% CoM trials - on group average waveforms)


%% do stats?

if ~exist('regTab','var'); load('./Saves/CPPAnalysis_CSD_.mat','regTab'); end
regTab.peakLat = nanzscore(col(peakLat));

% factors = {'confInR1','certainty'};
% nF = length(factors);
for iF = 1:nF

    f{iF,1} = fitglme(regTab,sprintf('peakLat ~ 1 + cond*%s+ (1 | pp)',factors{iF}));
    
    % split by cond
    for iC = 1:2
        f1{iF,iC} = fitglme(regTab((regTab.cond>=0)==(iC-1),:),sprintf('peakLat ~ 1 + %s+ (1 | pp)',factors{iF}));
    end
end

stats = StatsTableFromLMECells(f, factors);
stats2 = StatsTableFromLMECells(f1(:,1), factors,'Estimate',[],labels.cond(1)); % effect of fac
stats2 = [stats2; StatsTableFromLMECells(f1(:,2), factors,'Estimate',[],labels.cond(2))];

bic = cellfun(@(x) x.ModelCriterion.BIC, f);
bic2 = cellfun(@(x) x.ModelCriterion.BIC, f1);

%% include as covar on ampl~conf

% peakLat correl with amplPost i.e. late peak = high ampl, early peak = neg ampl
fitglme(regTab(regTab.cond>=0,:),'amplPost ~ 1 + peakLat*confInR1 + (1 | pp)');
fitglme(regTab(regTab.cond>=0,:),'amplPost ~ 1 + peakLat*certainty + (1 | pp)'); % interacts with certainty

% what is the interaction due to
u=unique(regTab.certainty(~isnan(regTab.certainty)));
for i = 1:3 % each certainty
    f3{1,i} = fitglme(regTab(regTab.cond>=0 & regTab.certainty==u(i),:),'amplPost ~ 1 + peakLat + (1 | pp)');
end
cellfun(@(x) x.Coefficients.Estimate(2), f3);
% steeper peak effect as certainty increases
% i.e. quicker decay with confidence?

%% compare continued RT to extinguished
% we said that 'certain' Cont was almost as fast as Ext, but was that
% overall or just certain Ext?

% i.e. compare cont vs ext per confinr1, which are n.s.?
u = unique(regTab.confInR1(~isnan(regTab.confInR1)));
for i = 1:6
    f4{i} = fitglme(regTab(regTab.confInR1==u(i),:), 'confRTLog ~ 1 + cond + (1 | pp)');
end
cellfun(@(x) x.Coefficients.pValue(2), f4); % all are sig diff
cellfun(@(x) x.Coefficients.Estimate(2), f4); % 

% or can compare each cont confinr1 to all-ext?
for i = 1:6
    % all Ext, plus continued just for this confInR1
    r = [regTab(regTab.cond<=0,:); regTab(regTab.cond>=0 & regTab.confInR1==u(i),:)];
    f5{i} = fitglme(r, 'confRTLog ~ 1 + cond + (1 | pp)');
end
cellfun(@(x) x.Coefficients.pValue(2), f5); % all are sig diff
cellfun(@(x) x.Coefficients.Estimate(2), f5); % 


%% we would predict that earlier peak latency is correlated with faster RT
factor = 'certainty';

minLat = 100; % exclude below these
peakLatCutOff = peakLat;
peakLatCutOff(peakLat<minLat) = NaN;
confRT = behData.confRT;
confRT(peakLat<minLat) = NaN;

figure();

subplot(2,2,1); % overall
conditionalPlot(peakLatCutOff', confRT');
xlabel('CPP peak latency after init RT (ms)');
ylabel('final RT (ms)');

% split by cond
subplot(2,2,2); % overall
x = groupMeans(peakLatCutOff,2,behData.cond,'dim');
y = groupMeans(confRT,2,behData.cond,'dim');
[~,~,~,~,h] = conditionalPlot(permute(x,[3,1,2]), permute(y,[3,1,2]));
xlabel('CPP peak latency after init RT (ms)');
ylabel('final RT (ms)');
legend(h(2:2:end), labels.cond,'Location','Best');

% split by confInR1

facByCond = groupMeans(behData.(factor),2,behData.cond,'dim');
x1 = groupMeans(x,3,facByCond,'dim');
y1 = groupMeans(y,3,facByCond,'dim');
for iC = 1:2
    subplot(2,2,2+iC); 
    set(gca,'ColorOrder',colours.(factor),'nextplot','replacechildren');
    conditionalPlot(permute(x1(:,iC,:,:),[4,1,3,2]), permute(y1(:,iC,:,:),[4,1,3,2]));
    h = findobj(gca,'Type','Patch');
    xlabel('CPP peak latency after init RT (ms)');
    ylabel('final RT (ms)');
    if iC==1; legend(flip(h), labels.(factor),'Location','Best'); end
    title(labels.cond(iC));
end

%% and greater decay (steeper neg-slope or lower final amplitudes)

% take ampl in final window?
finWin = [900 1000];
finalAmpl = sq(nanmean(cppFilt(:,isBetween(eeg.respTimes, finWin),:),2));
finalAmpl(peakLat<minLat) = NaN;

% or use slope?


figure();

subplot(2,2,1); % overall
conditionalPlot(peakLatCutOff', finalAmpl');
xlabel('CPP peak latency after init RT (ms)');
ylabel('final post-choice CPP amplitude(ms)');

% split by cond
subplot(2,2,2); % overall
x = groupMeans(peakLatCutOff,2,behData.cond,'dim');
y = groupMeans(finalAmpl,2,behData.cond,'dim');
[~,~,~,~,h] = conditionalPlot(permute(x,[3,1,2]), permute(y,[3,1,2]));
xlabel('CPP peak latency after init RT (ms)');
ylabel('final post-choice CPP amplitude(ms)');
legend(h(2:2:end), labels.cond,'Location','Best');

% split by confInR1
facByCond = groupMeans(behData.(factor),2,behData.cond,'dim');
x1 = groupMeans(x,3,facByCond,'dim');
y1 = groupMeans(y,3,facByCond,'dim');
for iC = 1:2
    subplot(2,2,2+iC); 
    set(gca,'ColorOrder',colours.(factor),'nextplot','replacechildren');
    conditionalPlot(permute(x1(:,iC,:,:),[4,1,3,2]), permute(y1(:,iC,:,:),[4,1,3,2]));
    h = findobj(gca,'Type','Patch');
    xlabel('CPP peak latency after init RT (ms)');
    ylabel('final post-choice CPP amplitude(ms)');
    if iC==1; legend(flip(h), labels.(factor),'Location','Best'); end
    title(labels.cond(iC));
end