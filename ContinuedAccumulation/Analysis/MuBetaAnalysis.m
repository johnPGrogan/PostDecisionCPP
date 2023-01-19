% MuBetaAnalysis
% run GetMuBEta first

if ~exist('topoplot','file'); eeglab nogui; end
%%
clc; clear; close all

% set these
useCSD = 1;
excludeBadPps = 1; % <640 trials
excludeTooFew = 1; % remove if <20 trials per conf3
excludeByRT = 1; % remove outside [100 1500]
excludeCoMFromCert = 0; % remove CoM trials from behData.certainty
doBaseline = 0; % baseline to lsat window before evonset, log

amplWindows = [-300 0; 700 1000]; % ms
slopeWindows = [-340 -150; 640 850];
cmap = 'jet';%crameri('vik');

% load
outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));
load(fullfile(outFolder, 'FlagArtefacts3.mat'),'isFlagged','ppsToExclude');
load(fullfile(outFolder, 'BehDataLoad.mat'),'behData','labels','behTab');

if useCSD
    f = load('./Saves/DoFFTCSD.mat');
    load(fullfile(outFolder, 'GetMuBetaCSD.mat'));
    fileInfo.fftFolder = './Saves/';%'D:\TCD\Projects\ContinuedAccumulation\Data\FFT_CSD';
    load(fullfile(fileInfo.fftFolder, 'PreRespMeanBetaCSD.mat'),'betaTopo');
    mapLims = [-1.3 1.3];

else
    f = load('./Saves/DoFFT.mat');
    load(fullfile(outFolder, 'GetMuBeta.mat'));
    fileInfo.fftFolder = './Saves';%'D:\TCD\Projects\ContinuedAccumulation\Data\FFT';
    load(fullfile(fileInfo.fftFolder, 'PreRespMeanBeta.mat'),'betaTopo');
    mapLims = [-.15 .15];

end

labels.hemispheres = {'ipsi','contra','contra-ipsi'};


%% exlude bad?

toRemove = sq(any(isnan(betas(:,:,:,3)),2)); % flagged trials or cpp outliers

if excludeBadPps
    toRemove(ppsToExclude,:) = 1;
end


if excludeByRT
    rtLims = [100 1540];
    toRemove(~isBetween(behData.RT, rtLims)) = 1;
end

if excludeCoMFromCert
    %%%% set any certainty with CoM to NaN - will be ignored everywhere
    behData.certainty(behData.CoM==1) = NaN;
end

if excludeTooFew

    % only remove those trials, not the entire person
    behData.certainty(toRemove) = NaN;
    toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.cond);

    toRemove(toRemove2) = 1;

%     % remove those people entirely, not just for that cond
%     toRemove(any(CountUnique(groupMeans(behData.certainty,2,behData.cond,'dim'),3)<10,[2 3]),:) = 1;


end


% actually exclude now
betas(repmat(permute(toRemove,[1,3,2]),1,size(betas,2),1,3)) = NaN;
betaStims(repmat(permute(toRemove,[1,3,2]),1,size(betaStims,2),1,3)) = NaN;
betaTopo(:, all(toRemove,2)) = NaN;

%% baseline?

if doBaseline
    baselineInds = find(f.stimWindows <= -250,1,'last');%isBetween(windows.stim, [-150 0]);
    baseline = nanmean(betaStims(:, baselineInds,:,:),2);

    betaStims = log( betaStims ./ baseline );
    betas = log( betas ./ baseline ); 

    % recalc lat ind
    betas(:,:,:,3) = diff(betas(:,:,:,1:2),[],4);
    betaStims(:,:,:,3) = diff(betaStims(:,:,:,1:2),[],4);

end


%%

factors = {'acc','certainty','CoM','confInR1','conf3'};
nF = length(factors);
cols = {[0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2];
        crameri('roma',6)};

%% calculate amplitudes and slopes 

betaVars.betaAmpl = sq(nanmean(betas(:,isBetween(f.respWindows,amplWindows(1,:)),:,:),2));
betaVars.betaAmplPost = sq(nanmean(betas(:,isBetween(f.respWindows,amplWindows(2,:)),:,:),2));

betaVars.betaSlope = FitCPPSlope(betas, slopeWindows(1,:), f.respWindows); % get slopes
betaVars.betaSlopePost = FitCPPSlope(betas, slopeWindows(2,:), f.respWindows); % get slopes

betaVarNames = fieldnames(betaVars);
nBetaVars = length(betaVarNames);

% invert lateralisation when they change their minds - recodes to be for
% final button pressed, not initial button

% x = cat(3, false(30,1280,2), behData.CoM==1); % only flip latInd, not ipsi or contra
% betaVars.betaAmplPost(x) = -betaVars.betaAmplPost(x);
% betaVars.betaSlopePost(x) = -betaVars.betaSlopePost(x);

% need to swap ipsi + contra rather than just inverting them all (won't
% affect lateralisation score)

for i = 1:nBetaVars
    y = betaVars.(betaVarNames{i}); % copy
    y = reshape(y,[],3); % reshape so can index across pp/tr in each hemi
    y(col(behData.CoM==1),[1 2]) = y(col(behData.CoM==1),[2 1]); % swap hemis for CoM
    
    betaVars1.(betaVarNames{i}) = reshape(y, size(betaVars.(betaVarNames{i}))); % back into shape

    % re-calc contra-ipsi
    betaVars1.(betaVarNames{i})(:,:,3) = diff(betaVars1.(betaVarNames{i})(:,:,1:2),[],3);
end

% should I also flip betas?
betas1 = betas;
betas1 = reshape(betas,[],3);
c = col(repmat(permute(behData.CoM==1,[1,3,2]),1,length(f.respWindows),1)); % in same shape
betas1(c,[1 2]) = betas1(c,[2 1]); % swap hemis
betas1 = reshape(betas1, size(betas));
betas1(:,:,:,3) = -diff(betas1(:,:,:,1:2),[],4); % contra-ipsi






%% does single-trial lme have better sensitivity?

if excludeCoMFromCert
    behTab.certainty(behTab.CoM>0) = NaN;
end

% remove trials
behTab(col(toRemove),:) = array2table(NaN(sum(toRemove,'all'),width(behTab)));

% vars
matVars = reshape(struct2array(structfun(@(x) col(x)', betaVars1, 'UniformOutput',0)),numel(behData.RT),[]);
if ~doBaseline
    matVars(:,[1 2 4 5]) = log(matVars(:,[1 2 4 5])); % ampl contra+ipsi are log-normal distrib
end
% matVars = nanzscore(matVars);

regTab = horzcat(behTab, array2table(matVars));
betaNames = strcat(col(repmat({'ampl','amplPost','slope','slopePost'},3,1))', col(repmat({'Ipsi','Contra','Lat'}',1,4))');
regTab.Properties.VariableNames(end-11:end) = betaNames;
regNames = regTab.Properties.VariableNames;

isLogistic = cellRegexpi(regNames, 'Logistic')>0;
regTab(:,~isLogistic) = varfun(@nanzscore, regTab(:,~isLogistic));

glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression, if behVars are DV
nF=3 - excludeCoMFromCert;

%% stats in the paper


fitglmeCell = @(x) fitglme(regTab, x);
% for eachcond sep, i=0 | 1
fitglmeCellSep = @(x, i) fitglme(regTab(regTab.cond>0 == i,:), x); 

if ~excludeCoMFromCert
    paperFormulae = {'%s ~ 1 + cond*acc + (1 | pp)';
                     '%s ~ 1 + cond*certainty + (1 | pp)';
                     '%s ~ 1 + cond*CoM + (1 | pp)';
                     '%s ~ 1 + cond*acc*certainty + (1 | pp)';
                     '%s ~ 1 + cond*confInR1 + (1 | pp)';};
    paperFormulaeSep = {'%s ~ 1 + confInR1 + (1 | pp)'}; % in each cond

    for i = 1:3
        dvName = betaNames{i+3}; % amplPostIpsi,contra,lat
        paperFormulae1 = cellfun(@(x) sprintf(x, dvName), paperFormulae, 'UniformOutput',0);
        paperFormulaeSep1 = cellfun(@(x) sprintf(x, dvName), paperFormulaeSep, 'UniformOutput',0);
    
        paperFits(:,i) = cellfun(fitglmeCell, paperFormulae1, 'UniformOutput',0);

        paperSepFits(:,i) = cellfun(fitglmeCellSep, repmat(paperFormulaeSep1,2,1), {0;1}, 'UniformOutput',0);
    end

else % do first-order (for supplementary)

    paperFormulae = {'%s ~ 1 + cond*acc + (1 | pp)';
                     '%s ~ 1 + cond*certainty + (1 | pp)';
                     '%s ~ 1 + cond*acc*certainty + (1 | pp)';};
    
    paperFormulaeSep = {'%s ~ 1 + acc*certainty + (1 | pp)'}; % in each cond
    for i = 1:3
        dvName = betaNames{i+3}; % amplPostIpsi,contra,lat
        paperFormulae1 = cellfun(@(x) sprintf(x, dvName), paperFormulae, 'UniformOutput',0);
        paperFormulaeSep1 = cellfun(@(x) sprintf(x, dvName), paperFormulaeSep, 'UniformOutput',0);
  
        paperPreFits(:,i) = cellfun(fitglmeCell, paperFormulae1, 'UniformOutput',0);

        paperPreSepFits(:,i) = cellfun(fitglmeCellSep, repmat(paperFormulaeSep1,2,1), {0;1}, 'UniformOutput',0);
    end

end

%% save now

saveOpts = {'Volt','CSD'; '', 'ExclCoMFromCert'; '', 'BS'};
saveName = sprintf('MuBetaAnalysis_%s_%s_%s.mat', saveOpts{1,useCSD+1}, saveOpts{2, excludeCoMFromCert+1}, saveOpts{3, doBaseline+1});

save(fullfile(outFolder, saveName));

%%%% don't run the below figures, just call PlotPaperFigures.m

%% plot per pp

% figure();
% [r,c] = GetSubPlotShape(fileInfo.nPP);
% for iPP = 1:fileInfo.nPP
%     
%     subplot(r,c,iPP);
%     errorBarPlot(permute(betas(iPP,:,:,3),[3,2,1]), 'area',1, 'xaxisvalues',f.respWindows);
%     
%     xline(0,':k');
%     yline(0, ':k');
%     
%     title(sprintf('%s: %s', fileInfo.ppID{iPP}, betaChans{iPP}));
%     
% end
% 
% %% plot all included pps together
% 
% figure();
% h = errorBarPlot(permute(betas(:,:,:,3),[3,2,1]),'area',1,'xaxisvalues',f.respWindows);
% 
% xline(0,':k');
% yline(0, ':k');
% ylabel('\muV'); xlabel('time from resp 1 (ms)');
% box off;
% 
% title('grand mean');
% 
% 

%% average over pps

figure();

h = errorBarPlot(sq(nanmean(betas(:,:,:,:),3)),'area',1,'xaxisvalues',f.respWindows);

xline(0,':k');
yline(0, ':k');
xlines(amplWindows, '-k');
xlines(slopeWindows, '-r');
ylabel('\muV'); xlabel('time from resp 1 (ms)');
box off;

title('grand mean');
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
stimlocked = groupMeans(betaStims,3,repmat(permute(behData.cond2,[1,3,2]),1,size(betaStims,2),1,3),'dim'); %[pp t cond 3 tr]
stimlocked = permute(stimlocked,[1,2,3,5,4]); %[pp t cond tr 3]
f.stimTimes = f.stimWindows(isBetween(f.stimWindows, win));

win = [-200 1200]; % for plotting
h = plot(f.stimTimes, reshape(nanmean(nanmean(stimlocked(:,isBetween(f.stimWindows,win),:,:,1:2),4),1),length(f.stimTimes),[]),'LineWidth',2);
[h(1:2).LineStyle] = deal('--');
[h(1:2:end).Color] = deal(h(1).Color);
[h(2:2:end).Color] = deal(h(2).Color);
xline(0);
hold on;
h1 = emptyLegend(4, {{'--k'}, {'-k'}, {'-','Color', h(1).Color}, {'-', 'Color', h(2).Color}}, {'LineWidth',2}, {'ipsi','contra','continue','interrupt'}, {'Location','Best'});
title('evidence onset');
xlim(win);

subplot(2,2,3);
data = groupMeans(betas,3,repmat(permute(behData.cond2,[1,3,2]),1,size(betas,2),1,3),'dim'); %[pp t cond 3 tr]
data = permute(data,[1,2,3,5,4]);

win = [-500 1000];
x = isBetween(f.respWindows,win);
h = plot(f.respWindows(x), reshape(nanmean(nanmean(data(:,x,:,:,1:2),4),1),sum(x),[]),'LineWidth',2);
[h(1:2).LineStyle] = deal('--');
[h(1:2:end).Color] = deal(h(1).Color);
[h(2:2:end).Color] = deal(h(2).Color);
xline(0);
hold on;
h1 = emptyLegend(4, {{'--k'}, {'-k'}, {'-','Color', h(1).Color}, {'-', 'Color', h(2).Color}}, {'LineWidth',2}, {'ipsi','contra','continue','interrupt'}, {'Location','Best'});
title('first-order response');

% inset - contra minus ipsi for resplocked, by cond2
subplot(2,2,4);
win = [-500 1000];
x = isBetween(f.respWindows,win);
h = plot(f.respWindows(x), reshape(nanmean(nanmean(data(:,x,:,:,3),4),1),sum(x),[]),'LineWidth',2);
xline(0);
legend(h, labels.cond2,'Location','Best')
title('first-order response: contra-ipsi');

%% make fig3.7

% split [stimlocked cpp] by cond2, then [conf; acc; CoM]
% plot [continue interrupted [contra-ipsi both]



respWin = [-500 1000];
x = isBetween(f.respWindows,respWin);
ylims = [.78 1.1; .78 1.1; -.15 .1];
figure();
colours = {'b','r'};
lineStyles = {'--','-'};
lineWidths = {[2 .5];[.5 1.2 2];  [2 .5]; [.5 1.2 2]};

% split resplocked by cond2
resplocked = groupMeans(betas,3,repmat(permute(behData.cond2,[1,3,2]),1,size(betas,2),1,3),'dim'); %[pp t cond 3 tr]
nT = size(resplocked,2);
for i = 1:nF-1
    subplot(nF-1,3,i*3-2);
    
    % split by this factor
    facByCond = groupMeans(behData.(factors{i}), 2, behData.cond2,'dim');
    facByCond = repmat(permute(facByCond,[1,4,2,5,3]),1,nT,1,3,1);
    
    resplocked1 = groupMeans(resplocked,5,facByCond); %[pp t cond2 3 cond]
    
    % plot
    h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,1,1:2,:),1),sum(x),[]), 'Color', colours{1});
    [h(1:2:end).LineStyle] = deal(lineStyles{1});
    for j = 1:(length(h)/2); [h([j*2-1 j*2]).LineWidth] = deal(lineWidths{i}(j));end
    
    xline(0);%yline(0);
%     xlines(amplWindows, '-k');
%     xlines(slopeWindows, '-r');
    legend(h(2:2:end), labels.(factors{i}),'Location','Best');
    if i==1; title('Continue');end
    if i==nF; xlabel('time from first resp (ms'); end
    xlim(respWin); %ylim(ylims(2,:));
    ylabel('Beta (uV)');
    
    
    
    subplot(nF-1,3,i*3-1);
    % interrupted
    
    h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,2,1:2,:),1),sum(x),[]), 'Color', colours{2});
    [h(1:2:end).LineStyle] = deal(lineStyles{1});
    for j = 1:(length(h)/2); [h([j*2-1 j*2]).LineWidth] = deal(lineWidths{i}(j));end
    
    xline(0);%yline(0);
%     xlines(amplWindows, '-k');
%     xlines(slopeWindows, '-r');
    if i==1
        title('Interrupt');
        hold on;
        emptyLegend(2,{{'--k'},{'-k'}},{'LineWidth',2},{'ipsi','contra'},{'Location','Best'});
    end
    if i==nF; xlabel('time from first resp (ms'); end
    xlim(respWin); %ylim(ylims(2,:));
    ylabel('Beta (uV)');
    
    
    
    subplot(nF-1,3,i*3)
    % contra-ipsi
    h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,:,3,:),1),sum(x),[]));
    for j = 1:2; [h(j:2:end).Color] = deal(colours{j}); end
    for j = 1:(length(h)/2); [h([j*2-1 j*2]).LineWidth] = deal(lineWidths{i}(j));end
    
    xline(0);yline(0);
%     xlines(amplWindows, '-k');
%     xlines(slopeWindows, '-r');
    if i==1
        title('Contra-ipsi');
        hold on;
        emptyLegend(2,{{'-b'},{'-r'}},{'LineWidth',2},labels.cond2,{'Location','Best'});
    end
    if i==nF; xlabel('time from first resp (ms'); end
    xlim(respWin); %ylim(ylims(3,:));
    
    ylabel('Beta (uV)');


    makeSubplotScalesEqual(nF-1,3,(i*3-2):(i*3-1));
end

makeSubplotScalesEqual(nF-1,3,3:3:(3*(nF-1)));

%% do similar figure for conf6

lineWidths = {[.5 2]};
c = crameri('roma',6);
figure();

subplot(1,3,1);

% split by this factor
facByCond = groupMeans(behData.confInR1, 2, behData.cond,'dim');
facByCond = repmat(permute(facByCond,[1,4,2,5,3]),1,nT,1,3,1);

resplocked = groupMeans(betas,3,repmat(permute(behData.cond,[1,3,2]),1,size(betas,2),1,3),'dim'); %[pp t cond 3 tr]
resplocked1 = groupMeans(resplocked,5,facByCond); %[pp t cond 3 fac]

% plot
h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,1,1:2,:),1),sum(x),[]),'LineWidth',1.5);%, 'Color', colours{1});
[h(1:2:end).LineStyle] = deal(lineStyles{1});
for j = 1:(length(h)/2); [h([j*2-1 j*2]).Color] = deal(c(j,:)); end

xline(0);%yline(0);
legend(h(2:2:end), labels.confInR1,'Location','Best');
title('Continue');
xlabel('time from first resp (ms');
xlim(respWin); %ylim(ylims(2,:));
ylabel('Beta (uV)');



subplot(1,3,2);
% interrupted

h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,2,1:2,:),1),sum(x),[]),'LineWidth',1.5);%, 'Color', colours{1});
[h(1:2:end).LineStyle] = deal(lineStyles{1});
for j = 1:(length(h)/2); [h([j*2-1 j*2]).Color] = deal(c(j,:)); end

xline(0);%yline(0);
title('Interrupt');
hold on;
emptyLegend(2,{{'--'},{'-'}},{'LineWidth',2,'Color','k'},{'ipsi','contra'},{'Location','Best'});
xlabel('time from first resp (ms');
xlim(respWin); %ylim(ylims(2,:));
ylabel('Beta (uV)');



subplot(1,3,3);
% contra-ipsi
h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,:,3,:),1),sum(x),[]),'LineWidth',2);
for j = 1:2; [h(j:2:end).LineStyle] = deal(lineStyles{j}); end
for j = 1:(length(h)/2); [h([j*2-1 j*2]).Color] = deal(c(j,:)); end


xline(0);yline(0);
title('Contra-ipsi');
hold on;
emptyLegend(2, { {lineStyles{1}}, {lineStyles{2}}},{'Color','k','LineWidth',2}, labels.cond, {'Location','Best'})
xlabel('time from first resp (ms');
xlim(respWin); %ylim(ylims(3,:));

ylabel('Beta (uV)');


makeSubplotScalesEqual(1,3,1:2);

%% cont/int 
figure();
for i = 1:2
    subplot(1,2,i);
    set(gca,'ColorOrder',c,'nextplot','replacechildren');
    % contra-ipsi
    h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,i,3,:),1),sum(x),[]),'LineWidth',2);
%     h = errorBarPlot(sq(resplocked1(:,x,i,3,:)),'area',1,'xaxisvalues',f.respWindows(x));


%     for j = 1:(length(h)); [h(j).Color] = deal(c(j,:)); end
    
    
    xline(0);yline(0);
    title(['Contra-ipsi: ' labels.cond{i}]);
    
    xlabel('time from first resp (ms');
    xlim(respWin); %ylim(ylims(3,:));
    
    ylabel('Beta (uV)');
    if i==1
        legend(h, labels.confInR1,'Location','Best'); 
%         legend([h{:,1}], labels.confInR1,'Location','Best'); 
    end

end
makeSubplotScalesEqual(1,2,1:2);

%% do ipsi/contra separately

figure();
for j = 1:2
    for i = 1:2
        subplot(2,2,(i-1)*2+j);
%         subplot(1,2,j);
        set(gca,'ColorOrder',c,'nextplot','replacechildren');
%         h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,i,j,:),1),sum(x),[]),'LineWidth',2);
        h = errorBarPlot(sq(resplocked1(:,x,i,j,:)),'area',1,'xaxisvalues',f.respWindows(x),'plotargs',{'LineWidth',2});
        for ii = 1:size(h,1)
            h{ii,2}.Visible = 'off';
        end
        xline(0); yline(0);
        title([labels.hemispheres{j} ' ' labels.cond{i}]);
        
        xlabel('time from first resp (ms');
        xlim(respWin); %ylim(ylims(3,:));
        
        ylabel('Beta (uV)');
%         if i==2 && j==2; legend(h, labels.confInR1,'Location','Best'); end
        if i==2 && j==2; legend([h{:,1}], labels.confInR1,'Location','Best'); end
    
    end
end
makeSubplotScalesEqual(2,2);

%% contra/ipsi on same fig, for each cond. average over everything else

% resplocked is [pp t cond 3 tr]
% just use that 
lineStyles = {'--','-'};

figure();
% for j = 1:2 % hemi
    for i = 1:2 % cond
%         subplot(2,2,(j-1)*2+i);
        subplot(1,2,i);
%         set(gca,'ColorOrder',c,'nextplot','replacechildren');
%         h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,i,j,:),1),sum(x),[]),'LineWidth',2);
        h = errorBarPlot(sq(nanmean(resplocked(:,x,i,1:2,:),5)),'area',1,'xaxisvalues',f.respWindows(x));
        h{1,1}.LineStyle = lineStyles{1};
        
        xline(0);
        title(labels.cond{i});
        
        xlabel('time from first resp (ms');
        xlim(respWin); %ylim(ylims(3,:));
        
        ylabel('Beta (uV)');
%         if i==2 && j==2; legend(h, labels.confInR1,'Location','Best'); end
        if i==1; legend([h{:,1}], labels.hemispheres,'Location','Best'); end
    
    end
% end
makeSubplotScalesEqual(1,2);

%% then split by cert (but not CoM)

% resplocked is [pp t cond 3 tr]
% split by cert

% split by this factor
facByCond = groupMeans(behData.certainty, 2, behData.cond,'dim');
facByCond = repmat(permute(facByCond,[1,4,2,5,3]),1,nT,1,3,1);

resplockedByCert = groupMeans(resplocked,5,facByCond); %[pp t cond 3 cert]
lineStyles = {'--','-'};
c = [0    0.4470    0.8; 0.8500    0.3250    0.0980; .0    0.7510         0;];
figure();
for i = 1:2 % cond
%         subplot(2,2,(j-1)*2+i);
    subplot(1,2,i);
    for j = 1:2 % hemi
        set(gca,'ColorOrder',c,'nextplot','replacechildren');
%         h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,i,j,:),1),sum(x),[]),'LineWidth',2);
        hold on;
        h = errorBarPlot(sq(resplockedByCert(:,x,i,j,:)),'area',1,'xaxisvalues',f.respWindows(x),'plotargs',{'LineStyle',lineStyles{j},'LineWidth',2});
        for ii = 1:size(h,1)
            h{ii,2}.Visible = 'off';
        end
    end
    xline(0); yline(0);
    title(labels.cond{i});
    
    xlabel('time from first resp (ms');
    xlim(respWin); %ylim(ylims(3,:));
    
    ylabel('Beta (uV)');
%         if i==2 && j==2; legend(h, labels.confInR1,'Location','Best'); end
    if i==1
        legend([h{:,1}], labels.certainty,'Location','Best');
    else 
        emptyLegend(2, { {'--k'}, {'-k'}}, {'LineWidth', 2}, labels.hemispheres, {'Location','Best'}); 
    end
    
end
makeSubplotScalesEqual(1,2);

%% contra-ipsi by cond/acc/conf


% split by this factor
factor = 'confInR1'; % adjust here
facByCond = groupMeans(behData.(factor), 2, behData.cond,'dim');
facByCond = repmat(permute(facByCond,[1,4,2,5,3]),1,nT,1,3,1);


resplocked1 = groupMeans(resplocked,5,facByCond,'dim'); %[pp t cond 3 fac tr]

accByCond = groupMeans(behData.acc,2,behData.cond,'dim');
accByCondFac = groupMeans(accByCond,3,sq(facByCond(:,1,:,1,:)),'dim');

resplocked1 = groupMeans(resplocked1, 6, repmat(permute(accByCondFac,[1,5,2,6,3,4]),1,size(resplocked1,2),1,3,1,1)); %[pp time cond 3 conf acc]

figure();
for j = 1:2
    for i = 1:2
        subplot(2,2,(j-1)*2+i);
        set(gca,'ColorOrder',crameri('roma',6),'nextplot','replacechildren');
        
        % contra-ipsi
        h = plot(f.respWindows(x), reshape(nanmean(resplocked1(:,x,i,3,:,j),1),sum(x),[]),'LineWidth',2);
%         h = errorBarPlot(sq(resplocked1(:,x,i,3,:,j)),'area',1,'xaxisvalues',f.respWindows(x));
                
        xline(0);yline(0);
        title([labels.cond{i} ' ' labels.acc{j}]);
        
        xlabel('time from first resp (ms');
        xlim(respWin); %ylim(ylims(3,:));
        
        ylabel('Beta lateralisation');
        if i==2 && j==2; legend(h, labels.(factor),'Location','Best'); end
%         if i==2 && j==2; legend([h{:,1}], labels.(factor),'Location','Best'); end
    
    end
end
makeSubplotScalesEqual(2,2);


%% split by conf*acc*com*cert?

iF = find(strcmp(factors,'confInR1')); % confInR1, will split into cert*CoM later

% resplocked1 is pp time cond hemi conf6 acc

% need to change order of conf6
resplocked2 = permute(resplocked1,[1,2,3,6,5,4]); % pp time cond acc conf6 hemi
resplocked2 = resplocked2(:,x,:,:,[4 5 6 3 2 1],:);
resplocked2 = reshape(resplocked2,30,[],2,2,3,2,3); %[pp time cond acc cert com hemi


iF = 2; % reset to cert
figure();
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        set(gca,'ColorOrder',cols{iF}(1:end-1,:),'nextplot','replacechildren');
        hold on;
        for k = 1:2
%             h = errorBarPlot(sq(resplocked2(:,:,i,j,:,k,3)), 'area',1,'xaxisvalues',f.respWindows(x),'plotargs',{'LineStyle',lineStyles{k}});
            h = plot(f.respWindows(x), sq(nanmean(resplocked2(:,:,i,j,:,k,3),1)), 'LineWidth',2,'LineStyle',lineStyles{k});
        end
        
        title([labels.cond{i}, ', ' labels.acc{j}]);
        xline(0);yline(0);
        ylabel(factors{iF}); 

%         % perm testing - get corr-Err, Com-none, interaction
%         clear diffWaves;
%         diffWaves(:,:,1) = sq(diff(nanmean(cppFiltCondAccCertCoM(:,:,j,:,1:2),5),[],4)); % av over fac, effect of acc
%         diffWaves(:,:,2) = sq(diff(nanmean(cppFiltCondAccCertCoM(:,:,j,:,1:2),4),[],5)) .* diffFlips(iF); % av over acc, effect of fac
%         diffWaves(:,:,3) = sq(diff(diff(cppFiltCondAccCertCoM(:,:,j,:,1:2),[],5),[],4)) .* diffFlips(iF); % diff over fac, then acc
%         p = [];
%         for k = 1:3
%             [~,p(k,:)] = permutationOLS(diffWaves(:,:   ,k),[],[],[],clusterArgs{:});
%             pbar(p(k,:), 'xVals', eeg.respTimes(tInds), 'yVal',  min(ylim) + k/5, 'plotargs',{'LineWidth',5,'Color',circshift([1 0 0],k-1)});
%         end
        
    end
end

clear legStyles;
for i = 1:length(h)
    legStyles{1,i} = {'-','Color',cols{iF}(i,:)};
end
legStyles = [legStyles {{lineStyles{1},'Color', [.5 .5 .5]}, {lineStyles{2},'Color', [.5 .5 .5]}, {'-k'}} ];
legNames = [labels.(factors{iF}), labels.CoM ];
emptyLegend(length(legNames), legStyles, {'LineWidth',2}, legNames, {'Location','Best'});   

makeSubplotScalesEqual(2,2);

% will draw the means within window below




%% figs - means

% Certainty: no-CoM trials
% for iV = 1:nVars
%     varsNoCoM.(varNames{iV}) = vars.(varNames{iV});
%     varsNoCoM.(varNames{iV})(repmat(behData.CoM==1,1,1,3)) = NaN;
% end
% 
% varsNoCoMByCond = structfun(@(x) groupMeans(x, 2, repmat(behData.cond,1,1,3),'dim'), varsNoCoM, 'UniformOutput',0);

nBetaVars = 2; % no slopes

betaVarsByCond = structfun(@(x) groupMeans(x, 2, repmat(behData.cond,1,1,3),'dim'),betaVars,'UniformOutput',0);

% factors{5} = 'confInR1';
nF = 3;%length(factors);
figure();
ylims = [-.1 0; -.15 .1; 0.8 1; 0.9 1.1];
for iF = 1:nF
    % split by cond %[pp cond 3 tr]
    factor = repmat(permute(groupMeans(behData.(factors{iF}),2,behData.cond, 'dim'),[1,2,4,3]),1,1,3);

    for iV = 1:nBetaVars%1:nBetaVars

        % split by factor
        data = groupMeans(betaVarsByCond.(betaVarNames{iV}), 4, factor); %[pp cond 3 factor]
        % plot means
        
        subplot(nBetaVars, nF, (iV-1)*nF+iF);
        h = errorBarPlot(sq(data(:,:,3,:)), 'type','bar');
        for i = 1:length(h)
            h(i).FaceColor = cols{iF}(i,:);
        end

        set(gca,'XTick', 1:2, 'XTickLabel', labels.cond);
%         ylim(ylims(iV,:));
        if iF==1; ylabel(betaVarNames{iV});end
        if iV==1; title(factors{iF}); legend(h, labels.(factors{iF}),'Location','Best');end
    end
    
end

for iV = 1:nBetaVars
    makeSubplotScalesEqual(nBetaVars,nF,(iV*nF-(nF-1)):iV*nF);
end

%% figs - ipsi & contra

% Certainty: no-CoM trials
% for iV = 1:nVars
%     varsNoCoM.(varNames{iV}) = vars.(varNames{iV});
%     varsNoCoM.(varNames{iV})(repmat(behData.CoM==1,1,1,3)) = NaN;
% end
% 
% varsNoCoMBycond = structfun(@(x) groupMeans(x, 2, repmat(behData.cond,1,1,3),'dim'), varsNoCoM, 'UniformOutput',0);

figure();

if useCSD
    ylims = [4.5 7.5; 4.5 7.5; -.004 .0002; -.003 .001]; 
    if doBaseline; ylims = [-.3 0; -.3 0; -.003 0; -.003 0]; end
else
    ylims = [.7 1; .7 1.1; -.0004 .00002; -.0002 .0001]; 
end
for iF = 1:nF
    % split by cond %[pp cond 3 tr]
    factor = repmat(permute(groupMeans(behData.(factors{iF}),2,behData.cond, 'dim'),[1,2,4,3]),1,1,3);

    for iV = 1:nBetaVars

        % split by factor
        data = groupMeans(betaVarsByCond.(betaVarNames{iV}), 4, factor); %[pp cond 3 factor]
        % plot means
        
        subplot(nBetaVars, nF, (iV-1)*nF+iF);
        h = errorBarPlot(reshape(permute(data(:,:,1:2,:),[1,2,4,3]),fileInfo.nPP,2,[]), 'type','bar'); % [cond [ipsi/contra factor]
        n = length(h)/2;
        [h(1:n).LineStyle] = deal('--'); % ipsi
        [h(1:n).FaceAlpha] = deal(0.5);
        for i = 1:length(h)
            h(i).FaceColor = cols{iF}(mod(i-1,n)+1,:);
        end

        set(gca,'XTick', 1:2, 'XTickLabel', labels.cond);
        ylim(ylims(iV,:));
        if iF==1; ylabel(betaVarNames{iV});end
        if iV==1; title(factors{iF}); legend(h((n+1):end), labels.(factors{iF}),'Location','Best');end
        if iV==2; legend(h([1 n+1]), {'ipsi','contra'},'Location','Best'); end
    end
    
end

for iV = 1:nBetaVars
    makeSubplotScalesEqual(nBetaVars,nF,(iV*nF-(nF-1)):iV*nF);
end

%% acc*fac


behDataByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),behData, 'UniformOutput',0);
betaVarsByCondAcc = structfun(@(x) groupMeans(x, 4, repmat(permute(behDataByCond.acc,[1,2,4,3]),1,1,3),'dim'),betaVarsByCond,'UniformOutput',0);
behDataByCondAcc = structfun(@(x) groupMeans(x, 3, behDataByCond.acc,'dim'),behDataByCond, 'UniformOutput',0);
%[pp cond 3 acc tr]


for iV = 1:4
    figure();
    for iF = 1:2 % cert + CoM
        for iC = 1:2 % cond
            
            subplot(2,2,(iF-1)*2+iC);
            set(gca,'ColorOrder',cols{iF+1},'nextplot','replacechildren');
            data = groupMeans(betaVarsByCondAcc.(betaVarNames{iV}), 5, repmat(permute(behDataByCondAcc.(factors{iF+1}),[1,2,5,3,4]),1,1,3)); %[pp cond 3 acc fac]
    
            h = errorBarPlot(sq(data(:,iC,3,:,:)),'type','bar');
            xticklabels(labels.acc);
            if iC==1; ylabel(factors{iF+1}); end
            if iF==1; title(labels.cond{iC}); end
        end
        legend(labels.(factors{iF+1}));
    end
    SuperTitle(betaVarNames{iV});
end

%% cond*acc*cert*com

% only post-dec
iF=2;

for iV = 2:2:4

    data = groupMeans(betaVarsByCondAcc.(betaVarNames{iV}), 5, repmat(permute(behDataByCondAcc.CoM,[1,2,5,3,4]),1,1,3),'dim'); %[pp cond 3 acc CoM tr]

    figure();
    factor = groupMeans(behDataByCondAcc.(factors{iF}), 4, behDataByCondAcc.CoM, 'dim'); %[pp cond acc Com tr]

    data2 = groupMeans(data, 6, repmat(permute(factor,[1,2,6,3,4,5]),1,1,3)); %[pp cond 3 acc com cert]

    for iC = 1:2 % cond
        for iA = 1:2 % acc
        
            subplot(2,2,(iA-1)*2+iC);

            set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
    
            h = errorBarPlot(sq(data2(:,iC,3,iA,:,:)),'type','bar');
            xticklabels(labels.CoM);
            if iC==1; ylabel(factors{iF}); end
            title([labels.cond{iC} ', ' labels.acc{iA}]);
        end
    end
    legend(labels.(factors{iF}));
    SuperTitle(betaVarNames{iV});

end

%% do stats

[anovas, anovasLat, anovasContra, anovasIpsi] = deal(struct());
for iF = 1:nF
    % split by cond %[pp cond 3 tr]
    factor = repmat(permute(groupMeans(behData.(factors{iF}),2,behData.cond, 'dim'),[1,2,4,3]),1,1,3);

    for iV = 1:nBetaVars

        % split by factor
        data = groupMeans(betaVarsByCond.(betaVarNames{iV}), 4, factor); %[pp cond 3 factor]
        anovas.(factors{iF}).(betaVarNames{iV}) = rmanova(data(:,:,1:2,:),{'pp','cond','hemi',factors{iF}},'categorical',2:4,'DummyVarCoding','effects');

        anovasLat.(factors{iF}).(betaVarNames{iV}) = rmanova(sq(data(:,:,3,:)),{'pp','cond',factors{iF}},'categorical',2:3,'DummyVarCoding','effects');
        
        anovasIpsi.(factors{iF}).(betaVarNames{iV}) = rmanova(sq(data(:,:,1,:)),{'pp','cond',factors{iF}},'categorical',2:3,'DummyVarCoding','effects');
        anovasContra.(factors{iF}).(betaVarNames{iV}) = rmanova(sq(data(:,:,2,:)),{'pp','cond',factors{iF}},'categorical',2:3,'DummyVarCoding','effects');
    end
    
    anovaTabs.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovas.(factors{iF}),'UniformOutput',0),'rowNames',anovas.(factors{iF}).betaAmpl.Term(2:end));
    anovaTabsLat.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovasLat.(factors{iF}),'UniformOutput',0),'rowNames',anovasLat.(factors{iF}).betaAmpl.Term(2:end));
    anovaTabsIpsi.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovasIpsi.(factors{iF}),'UniformOutput',0),'rowNames',anovasIpsi.(factors{iF}).betaAmpl.Term(2:end));
    anovaTabsContra.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovasContra.(factors{iF}),'UniformOutput',0),'rowNames',anovasContra.(factors{iF}).betaAmpl.Term(2:end));
end

% Certainty:
% ampl: hemi*certainty - not found
% amplPost: hemi*certainty - yes
% amplPost: cond*certainty - yes

% slopePost: hemi*cert*cond - no (but cert is sig

% acc:
% amplPost & slopePost: cond*hemi*acc - no (amplPost has hemi*acc tho)
% lat amplPost & slopePost: cond*acc

% CoM:
% ampl & slope: hemi*CoM - amplPost has hemi*CoM, slopes has CoM only
% Lat: amplPost & slopePost: ???




