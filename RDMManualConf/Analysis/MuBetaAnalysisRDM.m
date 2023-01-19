% MuBetaAnalysis
% run GetMuBEta first

if ~exist('topoplot','file'); eeglab nogui; end
%%
clc; clear; close all

% set these
useCSD = 1;
excludeBadPps = 1;
excludeTooFew = 1;
excludeByRT = 1;

amplWindows = [-250 -100; ]; % resplocked. respcuelocked was -200:0 & -300:0
slopeWindows = [-600 -300; ];
cmap = 'jet';%crameri('vik');

% load
outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));
load(fullfile(outFolder, 'FlagArtefacts.mat'),'isFlagged','ppsToExclude');
load(fullfile(outFolder, 'BehDataLoad.mat'),'behData','labels');

if useCSD
    f = load('./Saves/DoFFTCSD.mat');
    load(fullfile(outFolder, 'GetMuBetaCSD.mat'));
%     fileInfo.fftFolder = 'D:\TCD\Projects\RDMManualConf\Data\FFT_CSD';
    load(fullfile(fileInfo.outFolder, 'PreRespMeanBetaCSD.mat'),'betaTopo');
    mapLims = [-1.3 1.3];

else
    f = load('./Saves/DoFFT.mat');
    load(fullfile(outFolder, 'GetMuBeta.mat'));
%     fileInfo.fftFolder = 'D:\TCD\Projects\RDMManualConf\Data\FFT';
    load(fullfile(fileInfo.outFolder, 'PreRespMeanBeta.mat'),'betaTopo');
    mapLims = [-.15 .15];

end

labels.hemispheres = {'ipsi','contra','contra-ipsi'};

%% exlude bad?
toRemove = sq(any(isnan(betas(:,:,:,3)),2)); % flagged trials or cpp outliers

if excludeByRT
    rtLims = [100 1500];
    toRemove(~isBetween(behData.RT, rtLims)) = 1;
end

if excludeBadPps
    toRemove(ppsToExclude,:) = 1;
end

if excludeTooFew
    % only remove those trials, not the entire person
    behData.certainty(toRemove) = NaN;
    toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.stimDur);

    toRemove(toRemove2) = 1;
end

% actually exclude now
betas(repmat(permute(toRemove,[1,3,2]),1,size(betas,2),1,3)) = NaN;
betaStims(repmat(permute(toRemove,[1,3,2]),1,size(betaStims,2),1,3)) = NaN;
betaRespCues(repmat(permute(toRemove,[1,3,2]),1,size(betaRespCues,2),1,3)) = NaN;
betaTopo(:, all(toRemove,2)) = NaN;


%% set


factors = {'acc','certainty','confInCorr','confResp'};
nF = length(factors);

cols = {[0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2];
        crameri('roma',6);
        crameri('roma',6)};


%% calculate amplitudes and slopes 

betaVars.betaAmpl = sq(nanmean(betas(:,isBetween(f.respWindows,amplWindows(1,:)),:,:),2));

betaVars.betaSlope = FitCPPSlope(betas, slopeWindows(1,:), f.respWindows); % get slopes

betaVarNames = fieldnames(betaVars);
nBetaVars = length(betaVarNames);
%% does single-trial lme have better sensitivity?

behVarNames = {'pp','stimDur','corrLR','respLR','acc','RT','confResp','confInCorr','certainty'};
nBehVars = length(behVarNames);
behNames = fieldnames(behData); % names of all behData

behTabNames = behVarNames;
behTab = struct2table(structfun(@col, rmfield(behData, behNames(~ismember(behNames, behTabNames))), 'UniformOutput',0));

% vars
matVars = reshape(struct2array(structfun(@(x) col(x)', betaVars, 'UniformOutput',0)),numel(behData.RT),[]);
matVars(:,[1 2]) = log(matVars(:,[1 2])); % ampl contra+ipsi are log-normal distrib
% matVars = nanzscore(matVars);

regTab = horzcat(behTab, array2table(matVars));
betaNames = strcat(col(repmat({'ampl','slope'},3,1))', col(repmat({'Ipsi','Contra','Lat'}',1,2))');
regTab.Properties.VariableNames(end-5:end) = betaNames;
regNames = regTab.Properties.VariableNames;

% re-order names 
[~,b] = ismember([behTabNames, betaNames],regNames);
regTab = regTab(:, b);

regNames = regTab.Properties.VariableNames; % udpate

regTab.accLogistic = regTab.acc;

isLogistic = ismember(regNames, {'accLogistic'}); % which are [0 1] for logistic regressions
regTab(:,~isLogistic) = varfun(@nanzscore, regTab(:,~isLogistic));

glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression, if behVars are DV


%% paper stats

fitglmeCell = @(x) fitglme(regTab, x);
% for each acc sep, i=0 | 1
fitglmeCellSep = @(x, i) fitglme(regTab(regTab.acc>0 == i,:), x); 

paperFormulae = {'amplLat ~ 1 + stimDur*acc + (1 | pp)';
                 'amplLat ~ 1 + stimDur*certainty + (1 | pp)';
                 'amplLat ~ 1 + stimDur*certainty*acc + (1 | pp)';
                 };
paperFormulaeSep = {'amplLat ~ 1 + stimDur*certainty + (1 | pp)'}; % in each cond

% run
paperFits = cellfun(fitglmeCell, paperFormulae, 'UniformOutput',0);
paperSepFits = cellfun(fitglmeCellSep, repmat(paperFormulaeSep,2,1), {0;1}, 'UniformOutput',0);


%% save now

saveOpts = {'Volt','CSD'};
saveName = sprintf('MuBetaAnalysis_%s.mat', saveOpts{1,useCSD+1});

save(fullfile(outFolder, saveName));

% use ../../ContinuedAccumulation/Analysis/PlotPaperFigures.m to plot


%% plot per pp

figure();
[r,c] = GetSubPlotShape(fileInfo.nPP);
for iPP = 1:fileInfo.nPP
    
    subplot(r,c,iPP);
    errorBarPlot(permute(betas(iPP,:,:,3),[3,2,1]), 'area',1, 'xaxisvalues',f.respWindows);
    
    xline(0,':k');
    yline(0, ':k');
    
    title(sprintf('%s: %s', fileInfo.ppID{iPP}, betaChans{iPP}));
    
end

%% plot all included pps together

figure();
h = errorBarPlot(permute(betas(:,:,:,3),[3,2,1]),'area',1,'xaxisvalues',f.respWindows);

xline(0,':k');
yline(0, ':k');
ylabel('\muV'); xlabel('time from resp 1 (ms)');
box off;

title('grand mean');



%% average over pps

figure();

subplot(2,2,1);
h = errorBarPlot(sq(nanmean(betaStims(:,:,:,:),3)),'area',1,'xaxisvalues',f.stimWindows);
xline(0,':k');
yline(0, ':k');
xlines(amplWindows, '-k');
xlines(slopeWindows, '-r');
ylabel('\muV'); xlabel('time from stim(ms)');
box off;

title('grand mean');
legend([h{:,2}], labels.hemispheres, 'Location','Best');


subplot(2,2,2);
h = errorBarPlot(sq(nanmean(betaRespCues(:,:,:,:),3)),'area',1,'xaxisvalues',f.respCueWindows);
xline(0,':k');
yline(0, ':k');
xlines(amplWindows, '-k');
xlines(slopeWindows, '-r');
ylabel('\muV'); xlabel('time from resp cue (ms)');
box off;
title('respcue');
legend([h{:,2}], labels.hemispheres, 'Location','Best');

subplot(2,2,3);
h = errorBarPlot(sq(nanmean(betas(:,:,:,:),3)),'area',1,'xaxisvalues',f.respWindows);
xline(0,':k');
yline(0, ':k');
xlines(amplWindows, '-k');
xlines(slopeWindows, '-r');
ylabel('\muV'); xlabel('time from resp 1 (ms)');
box off;
title('resp locked');
legend([h{:,2}], labels.hemispheres, 'Location','Best');

%% Fig 3.6 

figure();


subplot(2,2,1); % topoplot of activity at RT

% try removing outliers etc
topoplot(nanmean(betaTopo,2),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'mapLimits',mapLims,...
    'emarker',{'.','k',10,1});
title('at RT')

subplot(2,2,2); % stimlocked trace by cond
win = [-200 1200];
stimlocked = groupMeans(betaStims,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(betaStims,2),1,3),'dim'); %[pp t cond 3 tr]
stimlocked = permute(stimlocked,[1,2,3,5,4]); %[pp t dur tr 3]
f.stimTimes = f.stimWindows(isBetween(f.stimWindows, win));
h = plot(f.stimTimes, reshape(nanmean(nanmean(stimlocked(:,isBetween(f.stimWindows,win),:,:,1:2),4),1),length(f.stimTimes),[]),'LineWidth',2);
[h(1:3).LineStyle] = deal('--');
[h(1:3:end).Color] = deal(h(1).Color);
[h(2:3:end).Color] = deal(h(2).Color);
[h(3:3:end).Color] = deal(h(3).Color);
xline(0);
hold on;
h1 = emptyLegend(5, {{'--k'}, {'-k'}, {'-','Color', h(1).Color}, {'-', 'Color', h(2).Color}, {'-', 'Color', h(3).Color}}, {'LineWidth',2}, [{'ipsi','contra'} labels.stimDur], {'Location','Best'});
title('evidence onset');
xlim(win);

subplot(2,2,3);
data = groupMeans(betaRespCues,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(betaRespCues,2),1,3),'dim'); %[pp t cond 3 tr]
data = permute(data,[1,2,3,5,4]);

win = [-500 500];
x = isBetween(f.respCueWindows,win);
h = plot(f.respCueWindows(x), reshape(nanmean(nanmean(data(:,x,:,:,1:2),4),1),sum(x),[]),'LineWidth',2);
[h(1:3).LineStyle] = deal('--');
[h(1:3:end).Color] = deal(h(1).Color);
[h(2:3:end).Color] = deal(h(2).Color);
[h(3:3:end).Color] = deal(h(3).Color);
xline(0);
hold on;
h1 = emptyLegend(5, {{'--k'}, {'-k'}, {'-','Color', h(1).Color}, {'-', 'Color', h(2).Color}, {'-', 'Color', h(3).Color}}, {'LineWidth',2}, [{'ipsi','contra'} labels.stimDur], {'Location','Best'});
title('response cue');


subplot(2,2,4);
data = groupMeans(betas,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(betas,2),1,3),'dim'); %[pp t cond 3 tr]
data = permute(data,[1,2,3,5,4]);

win = [-500 1000];
x = isBetween(f.respWindows,win);
h = plot(f.respWindows(x), reshape(nanmean(nanmean(data(:,x,:,:,1:2),4),1),sum(x),[]),'LineWidth',2);
[h(1:3).LineStyle] = deal('--');
[h(1:3:end).Color] = deal(h(1).Color);
[h(2:3:end).Color] = deal(h(2).Color);
[h(3:3:end).Color] = deal(h(3).Color);
xline(0);
hold on;
h1 = emptyLegend(5, {{'--k'}, {'-k'}, {'-','Color', h(1).Color}, {'-', 'Color', h(2).Color}, {'-', 'Color', h(3).Color}}, {'LineWidth',2}, [{'ipsi','contra'} labels.stimDur], {'Location','Best'});
title('response');


%% plot lateralisation as above

figure();


subplot(2,2,1); % topoplot of activity at RT

% try removing outliers etc
topoplot(nanmean(betaTopo,2),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'mapLimits',mapLims,...
    'emarker',{'.','k',10,1});
title('at RT')

subplot(2,2,2); % stimlocked trace by cond
win = [-200 1200];
stimlocked = groupMeans(betaStims,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(betaStims,2),1,3),'dim'); %[pp t cond 3 tr]
stimlocked = permute(stimlocked,[1,2,3,5,4]); %[pp t dur tr 3]
f.stimTimes = f.stimWindows(isBetween(f.stimWindows, win));
h = plot(f.stimTimes, reshape(nanmean(nanmean(stimlocked(:,isBetween(f.stimWindows,win),:,:,3),4),1),length(f.stimTimes),[]),'LineWidth',2);
xline(0);
hold on;
legend(labels.stimDur,'Location','Best');
title('evidence onset');
xlim(win);
ylabel('beta lateralisation');

subplot(2,2,3);
data = groupMeans(betaRespCues,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(betaRespCues,2),1,3),'dim'); %[pp t cond 3 tr]
data = permute(data,[1,2,3,5,4]);
win = [-500 500];
x = isBetween(f.respCueWindows,win);
h = plot(f.respCueWindows(x), reshape(nanmean(nanmean(data(:,x,:,:,3),4),1),sum(x),[]),'LineWidth',2);
xline(0);
hold on;
title('response cue');
ylabel('beta lateralisation');

subplot(2,2,4);
data = groupMeans(betas,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(betas,2),1,3),'dim'); %[pp t cond 3 tr]
data = permute(data,[1,2,3,5,4]);
win = [-500 1000];
x = isBetween(f.respWindows,win);
h = plot(f.respWindows(x), reshape(nanmean(nanmean(data(:,x,:,:,3),4),1),sum(x),[]),'LineWidth',2);
xline(0);
hold on;
title('response');
ylabel('beta lateralisation');
%% make fig3.7


factors = {'acc','certainty'};
nF = length(factors);

respWin = [-500 1000];
x = isBetween(f.respWindows,respWin);
ylims = [.78 1.1; .78 1.1; -.15 .1];
figure();
cols2 = {'b','r','y'};
lineStyles = {'--','-'};
lineWidths = {[2 .5]; linspace(.5,3,6); linspace(.5,3,6); [.5 1.2 2];[.5 1.2 2];};

% split resplocked by cond2
resplocked = groupMeans(betas,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(betas,2),1,3),'dim'); %[pp t cond 3 tr]
nT = size(resplocked,2);
for i = 1:nF
     % split by this factor
    facByDur = groupMeans(behData.(factors{i}), 2, behData.stimDur,'dim');
    facByDur = repmat(permute(facByDur,[1,4,2,5,3]),1,nT,1,3,1);

    resplocked1 = groupMeans(resplocked,5,facByDur); %[pp t stimDur 3 cond]
    
    for j = 1:3
        subplot(nF,4,i*4 - (4-j));

        % plot
        h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,j,1:2,:),1),sum(x),[]), 'Color', cols2{1});
        [h(1:2:end).LineStyle] = deal(lineStyles{1});
        for k = 1:(length(h)/2); [h([k*2-1 k*2]).LineWidth] = deal(lineWidths{i}(k));end

        xline(0);%yline(0);
    %     xlines(amplWindows, '-k');
    %     xlines(slopeWindows, '-r');
        if j==1; legend(h(2:2:end), labels.(factors{i}),'Location','Best'); end
        if i==1; title(labels.stimDur{j});end
        if i==nF; xlabel('time from resp (ms'); end
        xlim(respWin); %ylim(ylims(2,:));
        ylabel('Beta (uV)');
    end
    
    
    
    
    subplot(nF,4,i*4)
    % contra-ipsi
    h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,:,3,:),1),sum(x),[]));
    for j = 1:3; [h(j:3:end).Color] = deal(cols2{j}); end
    for j = 1:(length(h)/3); [h([j*2-1 j*2]).LineWidth] = deal(lineWidths{i}(j));end
    
    xline(0);yline(0);
%     xlines(amplWindows, '-k');
%     xlines(slopeWindows, '-r');
    if i==1
        title('Contra-ipsi');
        hold on;
        emptyLegend(3,{{'-b'},{'-r'},{'-','Color',h(3).Color}},{'LineWidth',2},labels.stimDur,{'Location','Best'});
    end
    if i==nF; xlabel('time from resp (ms'); end
    xlim(respWin); %ylim(ylims(3,:));
    
    ylabel('Beta (uV)');


%     makeSubplotScalesEqual(nF,4,(i*4-3):(i*4-1));
end

% makeSubplotScalesEqual(nF,4,4:4:(4*nF));

%% do similar figure for conf6

lineWidths = {[.5 2]};
c = crameri('roma',6);
figure();


% split by this factor
facByDur = groupMeans(behData.confInCorr, 2, behData.stimDur,'dim');
facByDur = repmat(permute(facByDur,[1,4,2,5,3]),1,nT,1,3,1);

resplocked = groupMeans(betas,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(betas,2),1,3),'dim'); %[pp t dur 3 tr]
resplocked1 = groupMeans(resplocked,5,facByDur); %[pp t dur 3 cond]

% plot
for iD = 1:3
    subplot(2,2,iD);
    set(gca,'ColorOrder', c, 'nextplot','replacechildren');
    h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,iD,3,:),1),sum(x),[]),'LineWidth',1.5);%, 'Color', cols{1});
%     [h(1:2:end).LineStyle] = deal(lineStyles{1});
%     for j = 1:(length(h)/2); [h([j*2-1 j*2]).Color] = deal(c(j,:)); end

    xline(0);%yline(0);
    if iD==1; legend(h, labels.confInCorr,'Location','Best'); end
    title(labels.stimDur{iD});
    xlabel('time from first resp (ms');
    xlim(respWin); %ylim(ylims(2,:));
    ylabel('Beta lateralisation (uV)');
end


% subplot(2,2,4);
% % contra-ipsi
% h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,:,3,:),1),sum(x),[]),'LineWidth',2);
% for j = 1:3; [h(j:3:end).LineStyle] = deal(lineStyles{j}); end
% for j = 1:(length(h)/3); [h([j*2-1 j*2]).Color] = deal(c(j,:)); end
% 
% 
% xline(0);yline(0);
% title('Contra-ipsi');
% hold on;
% emptyLegend(2, { {lineStyles{1}}, {lineStyles{2}}},{'Color','k','LineWidth',2}, labels.stimDur, {'Location','Best'})
% xlabel('time from resp (ms');
% xlim(respWin); %ylim(ylims(3,:));
% 
% ylabel('Beta (uV)');
% 
% 
makeSubplotScalesEqual(2,2,1:3);

%% do ipsi/contra separately

figure();
for j = 1:2
    for i = 1:3
        subplot(2,3,(j-1)*3+i);
        set(gca,'ColorOrder',c,'nextplot','replacechildren');
%         h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,i,j,:),1),sum(x),[]),'LineWidth',2);
        h = errorBarPlot(sq(resplocked1(:,x,i,j,:)),'area',1,'xaxisvalues',f.respWindows(x));
        
        xline(0);%yline(1);
        title([labels.hemispheres{j} ' ' labels.stimDur{i}]);
        
        xlabel('time from first resp (ms');
        xlim(respWin); %ylim(ylims(3,:));
        
        ylabel('Beta (uV)');
%         if i==2 && j==2; legend(h, labels.confInR1,'Location','Best'); end
        if i==2 && j==2; legend([h{:,1}], labels.confInCorr,'Location','Best'); end
    
    end
end
makeSubplotScalesEqual(2,3);



%% contra-ipsi by cond/acc/cert


% split by this factor
factor = 'certainty'; % adjust here
facByDur = groupMeans(behData.(factor), 2, behData.stimDur,'dim');
facByDur = repmat(permute(facByDur,[1,4,2,5,3]),1,nT,1,3,1);


resplocked1 = groupMeans(resplocked,5,facByDur,'dim'); %[pp t dur 3 fac tr]

accByDur = groupMeans(behData.acc,2,behData.stimDur,'dim');
accByDurFac = groupMeans(accByDur,3,sq(facByDur(:,1,:,1,:)),'dim');

resplocked1 = groupMeans(resplocked1, 6, repmat(permute(accByDurFac,[1,5,2,6,3,4]),1,size(resplocked1,2),1,3,1,1)); %[pp time dur 3 conf acc]

figure();
for j = 1:2
    for i = 1:3
        subplot(2,3,(j-1)*3+i);
        set(gca,'ColorOrder',[0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410],'nextplot','replacechildren');
        
        % contra-ipsi
        h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,i,3,:,j),1),sum(x),[]),'LineWidth',2);
%         h = errorBarPlot(sq(resplocked1(:,x,i,3,:,j)),'area',1,'xaxisvalues',f.respWindows(x));
                
        xline(0);yline(0);
        title([labels.stimDur{i} ' ' labels.acc{j}]);
        
        xlabel('time from first resp (ms');
        xlim(respWin); %ylim(ylims(3,:));
        
        ylabel('Beta lateralisation');
        if i==2 && j==2; legend(h, labels.(factor),'Location','Best'); end
%         if i==2 && j==2; legend([h{:,1}], labels.(factor),'Location','Best'); end
    
    end
end
makeSubplotScalesEqual(2,3);


%% collapse across durations

figure();

for iF = 1:nF-1
    subplot(2,2,iF);
    set(gca,'ColorOrder', cols{iF}, 'nextplot', 'replacechildren');

    resplocked1 = groupMeans(betas(:,:,:,3),3,repmat(permute(behData.(factors{iF}),[1,3,2]),1,size(betas,2))); % [pp time fac]
    h = errorBarPlot(resplocked1, 'area',1, 'xaxisvalues',f.respWindows);

    % remove error bars
%     for i = 1:size(h,1); h{i,2}.Visible = 'off';  end

    xline(0);
    yline(0);

    xlabel('time to resp (ms)');
    ylabel('beta lateralisation');

    legend([h{:,1}], labels.(factors{iF}), 'Location','Best');


end

%% now do contra-ipsi sep


figure();
for iF = 1:nF-1
    
    resplocked1 = groupMeans(betas(:,:,:,1:2),3,repmat(permute(behData.(factors{iF}),[1,3,2]),1,size(betas,2),1,2)); % [pp time fac hemi]
    for iH = 1:2
        subplot(nF-1,2,(iF-1)*2+iH);
        set(gca,'ColorOrder', cols{iF}, 'nextplot', 'replacechildren');
        h = errorBarPlot(resplocked1(:,:,:,iH), 'area',1, 'xaxisvalues',f.respWindows);

        % remove error bars
        for i = 1:size(h,1); h{i,2}.Visible = 'off';  end
    
        xline(0);
    
        xlabel('time to resp (ms)');
        ylabel('beta lateralisation');
    
        if iH==1; legend([h{:,1}], labels.(factors{iF}), 'Location','Best'); end

        if iF==1; title(labels.hemispheres{iH}); end        

    end
    makeSubplotScalesEqual(nF-1,2,iF*2-1:iF*2);
end



%% figs - means

betaVarsByDur = structfun(@(x) groupMeans(x, 2, repmat(behData.stimDur,1,1,3),'dim'),betaVars,'UniformOutput',0);

figure();    
for iF = 1:nF
    % split by dur %[pp dur 3 tr]
    factor = repmat(permute(groupMeans(behData.(factors{iF}),2,behData.stimDur, 'dim'),[1,2,4,3]),1,1,3);

    for iV = 1:nBetaVars

        % split by factor
        data = groupMeans(betaVarsByDur.(betaVarNames{iV}), 4, factor); %[pp dur 3 factor]
        % plot means
        
        subplot(nBetaVars, nF, (iV-1)*nF+iF);
        h = errorBarPlot(sq(data(:,:,3,:)), 'type','bar');
        for i = 1:length(h)
            h(i).FaceColor = cols{iF}(i,:);
        end

        set(gca,'XTick', 1:3, 'XTickLabel', labels.stimDur);
%         ylim(ylims(iV,:));
        if iF==1; ylabel(betaVarNames{iV});end
        if iV==1; title(factors{iF}); else; legend(h, labels.(factors{iF}),'Location','Best');end
    end
    
end

for iV = 1:nBetaVars
    makeSubplotScalesEqual(2,nF,(iV*nF-(nF-1)):iV*nF);
end

%% figs - ipsi & contra

nF=2;
figure();
if useCSD
    ylims = [4 6; -.004 .0]; 
else
    ylims = [.7 1; .7 1.1; -.0004 .00002; -.0002 .0001]; 
end
for iF = 1:nF
    % split by dur %[pp dur 3 tr]
    factor = repmat(permute(groupMeans(behData.(factors{iF}),2,behData.stimDur, 'dim'),[1,2,4,3]),1,1,3);

    for iV = 1:nBetaVars

        % split by factor
        data = groupMeans(betaVarsByDur.(betaVarNames{iV}), 4, factor); %[pp cond 3 factor]
        % plot means
        
        subplot(nBetaVars, nF, (iV-1)*nF+iF);
        h = errorBarPlot(reshape(permute(data(:,:,1:2,:),[1,2,4,3]),fileInfo.nPP,3,[]), 'type','bar'); % [cond [ipsi/contra factor]
        n = length(h)/2;
        [h(1:n).LineStyle] = deal('--'); % ipsi
        [h(1:n).FaceAlpha] = deal(0.5);
        for i = 1:length(h)
            h(i).FaceColor = cols{iF}(mod(i-1,n)+1,:);
        end

        set(gca,'XTick', 1:3, 'XTickLabel', labels.stimDur);
        ylim(ylims(iV,:));
        if iF==1; ylabel(betaVarNames{iV});end
        if iV==1; title(factors{iF}); legend(h((n+1):end), labels.(factors{iF}),'Location','Best');end
        if iV==2; legend(h([1 n+1]), {'ipsi','contra'},'Location','Best'); end
    end
    
end

for iV = 1:nBetaVars
    makeSubplotScalesEqual(2,nF,(iV*nF-(nF-1)):iV*nF);
end

%% av over stimdur

nF=2;
figure();
if useCSD
    ylims = [4 6; -.004 .0]; 
else
    ylims = [.7 1; .7 1.1; -.0004 .00002; -.0002 .0001]; 
end
for iF = 1:nF
    % split by dur %[pp dur 3 tr]
%     factor = repmat(permute(groupMeans(behData.(factors{iF}),2,behData.stimDur, 'dim'),[1,2,4,3]),1,1,3);

    for iV = 1:nBetaVars

        % split by factor
        data = groupMeans(betaVars.(betaVarNames{iV})(:,:,1:2), 2, repmat(behData.(factors{iF}),1,1,2)); %[pp cond 3 factor]
        % plot means
        
        subplot(nBetaVars, nF, (iV-1)*nF+iF);
        set(gca,'ColorOrder', cols{iF},'nextplot','replacechildren');
        h = errorBarPlot(permute(data,[1,3,2]), 'type','bar'); % [cond [ipsi/contra factor]
        set(gca,'XTick', 1:2, 'XTickLabel', labels.hemispheres);
        ylim(ylims(iV,:));
        if iF==1; ylabel(betaVarNames{iV});end
        if iV==1; title(factors{iF}); legend(h, labels.(factors{iF}),'Location','Best');end
    end
    
end

for iV = 1:nBetaVars
    makeSubplotScalesEqual(2,nF,(iV*nF-(nF-1)):iV*nF);
end

%% acc * cert


betaVarsByAcc = structfun(@(x) groupMeans(x, 2, repmat(behData.acc,1,1,3),'dim'),betaVars,'UniformOutput',0); %[pp acc 3 tr
certByAcc = repmat(permute(groupMeans(behData.certainty,2,behData.acc, 'dim'),[1,2,4,3]),1,1,3);

figure();
for iV = 1:nBetaVars

    % split by factor
    data = groupMeans(betaVarsByAcc.(betaVarNames{iV}), 4, certByAcc); %[pp acc 3 factor]
    % plot means
    
    subplot(1, nBetaVars, iV);
    h = errorBarPlot(sq(data(:,:,3,:)), 'type','bar');
    for i = 1:length(h)
        h(i).FaceColor = cols{iF}(i,:);
    end

    set(gca,'XTick', 1:3, 'XTickLabel', labels.acc);
%         ylim(ylims(iV,:));
    if iV==1; ylabel(betaVarNames{iV});end
    if iV==1; title('certainty'); else; legend(h, labels.certainty,'Location','Best');end
end
    




%% do stats

[anovas, anovasLat, anovasContra, anovasIpsi] = deal(struct());
for iF = 1:nF
    % split by cond %[pp cond 3 tr]
    factor = repmat(permute(groupMeans(behData.(factors{iF}),2,behData.stimDur, 'dim'),[1,2,4,3]),1,1,3);

    for iV = 1:nBetaVars

        % split by factor
        data = groupMeans(betaVarsByDur.(betaVarNames{iV}), 4, factor); %[pp dur 3 factor]
        anovas.(factors{iF}).(betaVarNames{iV}) = rmanova(data(:,:,1:2,:),{'pp','stimDur','hemi',factors{iF}},'categorical',2:4,'DummyVarCoding','effects');

        anovasLat.(factors{iF}).(betaVarNames{iV}) = rmanova(sq(data(:,:,3,:)),{'pp','stimDur',factors{iF}},'categorical',2:3,'DummyVarCoding','effects');
        
        anovasIpsi.(factors{iF}).(betaVarNames{iV}) = rmanova(sq(data(:,:,1,:)),{'pp','stimDur',factors{iF}},'categorical',2:3,'DummyVarCoding','effects');
        anovasContra.(factors{iF}).(betaVarNames{iV}) = rmanova(sq(data(:,:,2,:)),{'pp','stimDur',factors{iF}},'categorical',2:3,'DummyVarCoding','effects');
    end
    
    anovaTabs.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovas.(factors{iF}),'UniformOutput',0),'rowNames',anovas.(factors{iF}).betaAmpl.Term(2:end));
    anovaTabsLat.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovasLat.(factors{iF}),'UniformOutput',0),'rowNames',anovasLat.(factors{iF}).betaAmpl.Term(2:end));
    anovaTabsIpsi.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovasIpsi.(factors{iF}),'UniformOutput',0),'rowNames',anovasIpsi.(factors{iF}).betaAmpl.Term(2:end));
    anovaTabsContra.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovasContra.(factors{iF}),'UniformOutput',0),'rowNames',anovasContra.(factors{iF}).betaAmpl.Term(2:end));
end


