function PlotNewFig6()
% make the post-choice CPP peak latency figures
% along with confidence-RT figure
% and motor lateralisation delay figure
% Requires:
% ContinuedAccumulation/Analysis/Saves/CPPAnalysis_CSD_.mat

%% defaults
gf = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',14);
if ~exist('filter_tf','file'); eeglab nogui;  end


%% set options

opts.useCSD = 1;
opts.excludeBadPps = 1; % remove pps with <640 good trials?
opts.excludeTooFew = 1; % remove pps with <20 per conf3
opts.excludeByRT = 1; % remove trials outside [100 1500] ms
opts.doFilt = 1; % whether to plot with loPass filtered data
opts.excludeCoMFromCert = 0; % remove CoM trials from behData.certainty

opts.outFolder = './Saves';

%% load

opts.saveOpts = {'Volt','CSD'; '', 'ExclCoMFromCert'};
opts.saveName = sprintf('CPPAnalysis_%s_%s.mat', opts.saveOpts{1,opts.useCSD+1}, opts.saveOpts{2, opts.excludeCoMFromCert+1});

optsNames = fieldnames(opts);
data = load(fullfile(opts.outFolder, opts.saveName), optsNames{:}, ...
    'cpp','behData','eeg','fileInfo','labels','factors');

%% check things match

dataNames = fieldnames(data);
data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));

if ~equals(opts, data1)
    warning('loaded data and options do not match');
    keyboard;
end

%%


data.respWin = [-100 1000];
data.showErrBars = 1; % error bars on traces?
data.faceAlpha = .1;
data.yl = [-30 35];
data.rInds = isBetween(data.eeg.respTimes, data.respWin);

colours.certainty = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];
% colours.certainty = [0.3731    0.8864    0.7382; 0.6115    0.5872    0.1107; 0.5505    0.0068    0.4520; .2 .2 .2];
colours.cond = [0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
colours.CoM = [0 0 1; 1 0 0; .2 .2 .2]; %[0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2]; % CoM
colours.conf3 = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];%colours.certainty;
% colours.conf3 = [0.0937    0.3636    0.6709; 0.6796    0.8444    0.6806; 0.5829    0.2780    0.0740; .2 .2 .2]; % av within confinr1 pairs
colours.confInR1 =  [flipud(crameri('roma',6)); .2 .2 .2];
colours.acc=colours.cond;
colours.cmap = crameri('vik');
data.colours = colours;    

widths.certainty = [2 2 2 2];
styles.certainty = {'-','-','-','--'};
widths.conf3 = [2 2 2 2];
styles.conf3 = {'-','-','-','--'};
widths.CoM = [2 2 2];
styles.CoM = {'-','-','--'};
widths.cond = [2 2 2];
styles.cond = {'-','-','--'};
widths.acc = widths.CoM;
styles.acc=  styles.CoM;
widths.confInR1 = [2 2 2 2 2 2];
styles.confInR1 = {'-','-','-','-','-','-', '--'};
data.lines.widths = widths;
data.lines.styles = styles;
data.alphas = table2struct(array2table([.2 .2 .2 .1 .2],'VariableNames',data.factors));

data.names.confInR1 = 'confidence-in-initial-choice';
data.names.certainty = 'final-certainty';
data.names.conf3 = 'confidence-in-initial-choice: binned';
data.names.CoM = 'change-of-mind';
data.names.acc = 'initial accuracy';



%% plot peak latency sorted + bar graph

subplotInds = [2 1 1; 2 2 3; 2 2 4]; 

plotExt = 1; % don't show extinguished in bar plots

figure();
PlotCPPPeakLat(data, subplotInds(1:2,:), plotExt);


%% plot conf-RT bar graph

PlotConfRT(data, subplotInds(3,:), plotExt);


%% plot beta latencies
% not used now
% PlotBetaLatency(data, subplotInds(4,:))


%%  
set(0,'DefaultAxesFontSize',gf); % reset

end

function PlotCPPPeakLat(data, subplotInds, plotExt)

struct2workspace(data);


[nPP, nT, nTr] = size(cpp);

%% filter

loPass = 6; % Hz

isAllNaN = repmat(all(isnan(cpp),2),1,size(cpp,2),1);
cpp(isAllNaN) = 0;   

% filter
cppFilt = eegfilt(reshape(permute(cpp,[1,3,2]),fileInfo.nPP*fileInfo.maxTr,size(cpp,2)),eeg.fs,0, loPass)'; % low-pass
cppFilt = permute(reshape(cppFilt,size(cpp,2),fileInfo.nPP,fileInfo.maxTr),[2,1,3]);

% set zeros back to NAN
cppFilt(isAllNaN) = NaN;
% cpp(isAllNaN) = NaN;

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



% smooth over trials
n=100;

%% split by cond

cond = col(behData.cond);
cond(toRemove)=[];

plotWin = [-200 1000];%peakWindow;
xInds = isBetween(eeg.respTimes,plotWin); % x-ticks to show
x = eeg.respTimes(xInds);
[~, xTickInds] = min(abs(x - [-200:200:max(x)]'),[],2);


iC = 2; % continued only


cppSortedCondSmooth = movmean(cppSorted(cond==iC,:), n, 1);
peakLatSortedCondSmooth = movmean(peakLatSorted(cond==iC,:), n, 1);


subplot(subplotInds(1,1),subplotInds(1,2),subplotInds(1,3));

colormap('jet');%(crameri('bam'));
imagesc(cppSortedCondSmooth(:,xInds),[-150 150]);
set(gca,'YDir','normal');
c = colorbar;
% c.Label.String = '\muV/m^2';
set(get(c,'Title'),'String', '\muV/m^2');

xlabel('time from init RT (ms)')
xticks(xTickInds); xticklabels(round(x(xTickInds),1,'significant'));
xline(find(x>0,1,'first'),':');
ylabel('trials');
yticks(1000:1000:length(peakLatSortedCondSmooth));

hold on;
plot(peakLatSortedCondSmooth + round((peakWindow(1)-plotWin(1))*.512), 1:numel(peakLatSortedCondSmooth), '-k');

title({'post-choice CPP Peak Latency';'Continued Evidence'});

%% plot peak latency by cond + confInR1


factor = 'confInR1';

subplot(subplotInds(2,1),subplotInds(2,2),subplotInds(2,3));

set(gca,'ColorOrder',colours.(factor),'nextplot','replacechildren');
peakLatByCondFac = reshape(groupMeans(peakLat,2,behData.cond + behData.(factor)*10),30,2,[]); %[pp cond fac]

if plotExt
    h = errorBarPlot(peakLatByCondFac,'type','bar');
    xticks(1:2); xticklabels(labels.cond);
    ylabel('post-choice CPP peak latency (ms)');
    ylim(minMax([h.YData]) .* [.9 1.1])
    legend(labels.(factor),'Location','Best');
    title('Peak Latency');

else % just continued
    h = errorBarPlot(sq(peakLatByCondFac(:,2,:)),'type','bar');
    h.FaceColor = 'flat';
    h.CData = colours.(factor)(1:end-1,:);
%     xticks(mean(xticks)); xticklabels(names.(factor));
    xticklabels(labels.(factor)); xtickangle(45);
    ylabel('post-choice CPP peak latency (ms)');
    ylim(minMax([h.YData]) .* [.9 1.1])
    
    c = cell(length(labels.(factor)),1);
    for i = 1:length(labels.(factor))
        c{i} = {'FaceColor', colours.(factor)(i,:)};
    end
    emptyLegend(length(labels.(factor)), c, {}, labels.(factor), {'Location','Best'}, @bar)
%     legend(labels.(factor),'Location','Best');
    title('Peak Latency');

end
     
end


function PlotConfRT(data, subplotInds, plotExt)

struct2workspace(data);

%% plot peak latency by cond + confInR1

factor = 'confInR1';

behDataByCond.confRT = groupMeans(behData.confRT, 2, behData.cond, 'dim');
behDataByCond.(factor) = groupMeans(behData.(factor),2,behData.cond,'dim');

confRTByCondFac = groupMeans(behDataByCond.confRT,3,behDataByCond.(factor)); %[pp cond conf]

subplot(subplotInds(1,1),subplotInds(1,2),subplotInds(1,3));

% % line plot
% set(gca,'ColorOrder',colours.cond,'nextplot','replacechildren');
% h = errorBarPlot(permute(confRTByCondFac,[1 3 2]),'type','bar','plotargs',{'LineWidth',2});
% xticks(1:6); xticklabels(labels.(factor));

set(gca,'ColorOrder',colours.(factor),'nextplot','replacechildren');

if plotExt
    h = errorBarPlot(confRTByCondFac,'type','bar');
    xticks(1:2); xticklabels(labels.cond);
    ylabel('final-response RT (ms)');
    % legend(h, labels.(factor), 'Location','Best');

else % just continued
    h = errorBarPlot(sq(confRTByCondFac(:,2,:)),'type','bar');
%     xticks(mean(xticks)); xticklabels(names.(factor));
    xticklabels(labels.(factor)); xtickangle(45);
    ylabel('final-response RT (ms)');
    % legend(h, labels.(factor), 'Location','Best');
    h.FaceColor = 'flat';
    h.CData = colours.(factor)(1:end-1,:);

end
title('final-response RT');

end

function PlotBetaLatency(data, subplotInds)
% extract just some
colours = data.colours;
behData = data.behData;
labels = data.labels;

% load others
load('./Saves/MuBetaAnalysis_CSD__.mat','betas','f','factors');

%%


% invert CoM trials so that it reflects final response
betas1 = reshape(betas,[],3);
c = col(repmat(permute(behData.CoM==1,[1,3,2]),1,length(f.respWindows),1)); % in same shape
betas1(c,[1 2]) = betas1(c,[2 1]); % swap hemis
betas1 = reshape(betas1, size(betas));
betas1(:,:,:,3) = diff(betas1(:,:,:,1:2),[],4); % contra-ipsi
betas = betas1;
clear betas1;


rInds = isBetween(f.respWindows,[-500 1000]);
nT = sum(rInds);
   

iF=4; % confinr1
facByCond = groupMeans(behData.(factors{iF}), 2, behData.cond,'dim');
facByCond = repmat(permute(facByCond,[1,4,2,5,3]),1,nT,1,3,1);
   
resplocked = groupMeans(betas(:,rInds,:,:),3,repmat(permute(behData.cond,[1,3,2]),1,nT,1,3),'dim'); %[pp t cond 3 tr]
resplocked1 = groupMeans(resplocked,5,facByCond,'dim'); %[pp t cond hemi fac tr]

% use av traces per pp for now
resplockedMean = nanmean(resplocked1,6); 

t = f.respWindows(rInds);
t2 = t(t>=0);



% find at mean condition level, across pps, using t-tests
endLat = NaN(3,2);
v = endLat;
for iC = 1:2
    for iL = 1:size(endLat,1)

       % so take mean in last 100ms
       % do t-tests to see earliest consec point of difference
       y = sq(resplockedMean(:,t>=0,iC,3,iL)); %[t tr]
%        [~,p] = ttest(y, -1, 'tail','both'); % [1 t]

       v(iL,iC) = -0.9881; % mean from figure
       % or vs mean final value?
%        v(iL,iC) = nanmean(nanmean(y(:,t(t>=0)>=900),2),1); % mean final value for this level
%        v(iL,iC) = nanmean(nanmean(nanmean(resplockedMean(:,isBetween(t, [-150 -50]),iC,iH,4:6),5),2),1); % mean initResp level
       [~,p] = ttest(y, v(iL,iC), 'tail','both'); % [1 t] 

       pReg = findregions(p>.05); % find is sig
       d = find(diff(pReg,[],2)>=1,1,'first'); % length
       tInd = pReg(d,1); % start of that period (10 samples x 20ms = 200ms)       
%            tInd = find(p>.05,1,'first');
       if tInd
           endLat(iL,iC) = t2(tInd); % find last difference
       else
           if all(p<=.05)
               endLat(iL,iC) = NaN; %t2(end);
           elseif all(p>.05)
               endLat(iL,iC) = NaN; %t2(1);
           end
       end
   end
end

% seems too short
% % or find at cond-level within pps 
% endLatCondConf = sq(apply(2,@(X) max([NaN find(X<=-1,1,'first')]), resplockedMean(:,t>0,:,3,:))); %[pp cond factor]
% endLat = t2(round(sq(nanmean(endLatCondConf(:,:,1:3),1))'));


subplot(subplotInds(1,1),subplotInds(1,2),subplotInds(1,3));

set(gca,'ColorOrder',colours.confInR1,'nextplot','replacechildren');
h = errorBarPlot(sq(resplockedMean(:,:,2,3,1:3)),'area',1,'xaxisvalues',t);
xline(0);yline(0);
xlabel('time from init RT');
ylabel('beta lateralisation');

% if overlap, separate by 5ms for plotting
if length(unique(endLat(1:3,2))) < 3
    endLat(1:3,2) = endLat(1:3,2) + [0 -5 5]';
end
h1 = xlines(endLat(1:3,2),'--','LineWidth',2);
for i = 1:3
    h1(i).Color = colours.confInR1(i,:);
end

if length(unique(v(~isnan(v))))>1 % if diff levels used
    h2 = ylines(v(1:3,2),':','LineWidth',1.5);
    for i = 1:3
        h2(i).Color = colours.confInR1(i,:);
    end
else
    yline(unique(v(~isnan(v))),':k');
end

legend(([h{:,1}]),labels.confInR1,'Location','Best');
title('\beta Lateralisation latency');

end

