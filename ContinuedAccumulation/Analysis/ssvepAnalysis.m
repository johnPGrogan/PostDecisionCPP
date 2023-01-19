%% ssvepAnalysis

%% set folders etc

clc; clear all; close all;

useCPPSSVEP = 0; % load up SSVEP in CPP chans. otherwise use Oz

useCSD = 1; % CSD or voltage
excludeBadPps = 1; % remove pps with <640 good trials?
excludeTooFew = 1; % remove pps with <20 per conf3
excludeByRT = 1; % remove trials outside [100 1500] ms
excludeCoMFromCert = 0; % remove CoM trials from behData.certainty

outFolder = './Saves';
load('./Saves/ExtractEpochs.mat')
load('./Saves/FlagArtefacts3.mat','isFlagged','ppsToExclude');
load('./Saves/BehDataLoad.mat','behData','labels','behTab'); % for l/r

if useCPPSSVEP
    loadName1 = 'GetCPPSSVEP';
else
    loadName1 = 'GetSSVEP';
end
if useCSD
    load(sprintf('./Saves/%sCSD.mat', loadName1)); % get FFT info
else
    load(sprintf('./Saves/%s.mat', loadName1)); % get FFT info
end


lockNames{3} = 'preResp'; % for regressions
%% remove trials/pps

toRemove = sq(any(isnan(ssvep.resp),2)); % flagged trials or cpp outliers

if excludeBadPps
    toRemove(ppsToExclude,:) = 1;     
end

if excludeTooFew
    behData.certainty(toRemove) = NaN;
    toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.cond);

    toRemove(toRemove2) = 1;
end

if excludeByRT
    rtLims = [100 1540];
    toRemove(~isBetween(behData.RT, rtLims)) = 1;
end

%%%% 
if excludeCoMFromCert
    %%%% set any certainty with CoM to NaN - will be ignored everywhere
    behData.certainty(behData.CoM==1) = NaN;
end

%% actually exclude
for i = 1:2
    ssvep.(lockNames{i})(repmat(permute(toRemove,[1,3,2]),1,nTimes.(lockNames{i}),1)) = NaN;
end


%% topoplot stuff
% % uncomment to run
% 
% load('./Saves/DoSSVEPTopo15.mat');
% % load('./Saves/DoFFT.mat','stimWindows','respWindows'); % get FFT info - only for 30Hz
% 
% % remove flagged
% ssvepTopo(repmat(permute(toRemove,[1,3,4,2]),1,eeg.nChans,length(stimWindows),1)) = NaN;
% ssvepTopoResp(repmat(permute(toRemove,[1,3,4,2]),1,eeg.nChans,length(respWindows),1)) = NaN;
% 
% % un-baseline?
% if doBaseline==1
%     ssvepTopo = exp(ssvepTopo) .* baseline;
%     ssvepTopoResp = exp(ssvepTopoResp) .* baseline;
% end   
% 
% ssvepTopoRespCont = groupMeans(ssvepTopoResp,4, repmat(permute(behData.cond,[1,3,4,2]),1,eeg.nChans, length(respWindows),1));
% 
% % % animate topoplots
% % AnimateTopoplots(sq(nanmean(ssvepTopo,[4 1])), stimWindows, 1:length(stimWindows), eeg.chanlocs, [], [0 1], './Figs/SSVEP15TopoStim.gif');
% % 
% % % resp-locked
% % AnimateTopoplots(sq(nanmean(ssvepTopoResp,[4 1])), respWindows, 1:length(respWindows), eeg.chanlocs, [], [0 1], './Figs/SSVEP15TopoResp.gif');
% % 
% % % continued only
% % AnimateTopoplots(sq(nanmean(ssvepTopoRespCont,1)), respWindows, 1:length(respWindows), eeg.chanlocs, labels.cond, [0 1], './Figs/SSVEP15TopoRespCont.gif');
% 
% % draw a few times only
% t = [ -510 10 510];
% figure();
% for i = 1:length(t)
%     for j = 1:2
%         subplot(2,length(t),i +(length(t) * (j-1)))
%         topoplot(nanmean(ssvepTopoRespCont(:,:,respWindows==t(i),j),1),eeg.chanlocs, 'maplimits', [0 1]);
%         if j==1; title(t(i));end
%     end
% end
% c = colorbar('Location','NorthOutside');
% c.Position = [0.2893 0.4730 0.5036 0.0556];

%% set post-dec interrupted ev resp to NaN? (as there is no stimulus so is just noise)

ssvep.stim(repmat(permute(behData.cond==1,[1,3,2]),1,size(ssvep.stim,2))) = NaN;
ssvep.resp(repmat(permute(behData.cond==1,[1,3,2]),1,size(ssvep.resp,2))) = NaN;

%% take mean within windows?

wins.stim = [150 350];
wins.resp = [-400 -200; 700 1000];

ssvepMeans.stim = sq(nanmean(ssvep.stim(:,isBetween(windows.stim, wins.stim),:),2));
ssvepMeans.preResp = sq(nanmean(ssvep.resp(:,isBetween(windows.resp, wins.resp(1,:)),:),2));
ssvepMeans.resp = sq(nanmean(ssvep.resp(:,isBetween(windows.resp, wins.resp(2,:)),:),2));

%% set

factors = {'acc','certainty','CoM','confInR1','conf3'};
nF = length(factors);
cols = {[0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2];
        [crameri('roma',6); .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2];};

%% regressions

if excludeCoMFromCert
    behTab.certainty(behTab.CoM>0) = NaN;
end

% remove trials
behTab(col(toRemove),:) = array2table(NaN(sum(toRemove,'all'),width(behTab)));

% vars
ssvepTab = struct2table(structfun(@(x) col(x), ssvepMeans, 'UniformOutput',0));

regTab = horzcat(behTab, ssvepTab);
regNames = regTab.Properties.VariableNames;

isLogistic = cellRegexpi(regNames, 'Logistic')>0;
regTab(:,~isLogistic) = varfun(@nanzscore, regTab(:,~isLogistic));
glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression, if behVars are DV

%% paper stats

paperFormulae = {'preResp ~ 1 + certainty + (1 | pp)';
                 'resp ~ 1 + certainty + (1 | pp)';
                 'resp ~ 1 + confInR1 + (1 | pp)';
                 'resp ~ 1 + acc*confInR1 + (1 | pp)';};

paperFormulaeSep = {'resp ~ 1 + confInR1 + (1 | pp)'}; % to run in err/corr sep

fitglmeCell = @(x) fitglme(regTab, x);
% for eachcond sep, i=0 | 1
fitglmeCellSep = @(x, i) fitglme(regTab(regTab.acc>0 == i,:), x); 

paperFits = cellfun(fitglmeCell, paperFormulae,'UniformOutput',0);
paperSepFits = cellfun(fitglmeCellSep, repmat(paperFormulaeSep,2,1), {0;1}, 'UniformOutput',0);

%% save

saveOpts = {'Volt','CSD'; '', 'CPPChans'};
saveName = sprintf('SSVEPAnalysis_%s_%s.mat', saveOpts{1,useCSD+1}, saveOpts{2, useCPPSSVEP+1});

save(fullfile(outFolder, saveName));


%%%% don't run the below figures, just call PlotPaperFigures.m

%% plot grand average

figure();
for i = 1:2
    subplot(1,2,i)
    errorBarPlot(nanmean(ssvep.(lockNames{i}),3),'area',1,'xaxisvalues',windows.(lockNames{i}));
    xline(0);
    xlabel(sprintf('time from %s (ms)', lockNames{i}));
    ylabel('SSVEP amplitude');
    title(lockNames{i});

end

%% split by cond


figure();
for i = 1:2
    subplot(1,2,i)
    ssvepByCond.(lockNames{i}) = groupMeans(ssvep.(lockNames{i}), 3, repmat(permute(behData.cond,[1,3,2]),1,nTimes.(lockNames{i}),1), 'dim'); %[pp t cond tr]
    
    h = errorBarPlot(nanmean(ssvepByCond.(lockNames{i}),4),'area',1,'xaxisvalues',windows.(lockNames{i}));
    xline(0);
    xlabel(sprintf('time from %s (ms)', lockNames{i}));
    ylabel('SSVEP amplitude');
    legend([h{:,1}], labels.cond,'Location','Best');
    title(lockNames{i});
end
%% plot each pp - av over all trials

figure();
for i = 1:2
    subplot(1,2,i)
    errorBarPlot(permute(ssvep.(lockNames{i}),[3 2 1]),'area',1,'xaxisvalues',windows.(lockNames{i}));
    xline(0);
    xlabel(sprintf('time from %s (ms)', lockNames{i}));
    ylabel('SSVEP amplitude');
    title(lockNames{i})
end


%% subplot each pp by cond
i = 2; % just resp
figure();
for iPP = 1:fileInfo.nPP
    subplot(5,6,iPP);
    h = errorBarPlot(permute(ssvepByCond.(lockNames{i})(iPP,:,:,:),[4,2,3,1]),'area',1,'xaxisvalues',windows.(lockNames{i}));
    xline(0);
    xlabel(sprintf('time from %s (ms)', lockNames{i}));
    ylabel('SSVEP amplitude');
    title(iPP);
    if iPP==3; legend([h{:,1}], labels.cond,'Location','Best'); end
end

%% show some individual trials

iPP = 30;
figure(); 
for i = 1:2
    subplot(1,2,i);
    hold on;
    cols2 = get(gca,'ColorOrder');
    clear h;
    for j = 1:2
        h(:,j) = plot(windows.(lockNames{i}),sq(ssvepByCond.(lockNames{i})(iPP,:,j,1:50)),'Color',cols2(j,:));
    end
    xlabel(sprintf('time from %s (ms)', lockNames{i}));
    ylabel('SSVEP amplitude');
    xline(0);
    legend(h(1,:), labels.cond,'Location','Best');

end


%% RT median split

% 2 or 3 speed bins? uncomment one below to choose
binNames = {'fast','slow'};
% binNames = {'fast','medium','slow'};
p = linspace(0,100, length(binNames)+1);
p(1) = [];
prc = prctile(behData.RT, p,2);
binInds = sum(behData.RT > permute(prc,[1,3,2]),3);


figure();
for i = 1:2

    subplot(1,2,i);
    ssvepByRT.(lockNames{i}) = groupMeans(ssvep.(lockNames{i}), 3, repmat(permute(binInds,[1,3,2]),1,nTimes.(lockNames{i}),1),'dim'); %[npp t bins tr]

    h = errorBarPlot(nanmean(ssvepByRT.(lockNames{i}),4),'area',1,'xaxisvalues',windows.(lockNames{i}));
    xline(0,':k');
    yline(0, ':k');
    ylabel('\muV'); 
    xlabel(sprintf('time from %s (ms)', lockNames{i}));
    box off;
    title('evidence locked');
end


%% check a few other conds


figure();
for iF = 1:nF
    for i = 1:2
        subplot(nF,2,(iF-1)*2+i)
        set(gca,'ColorOrder', cols{iF}, 'nextplot','replacechildren');
        ssvepByFac = groupMeans(ssvep.(lockNames{i}), 3, repmat(permute(behData.(factors{iF}),[1,3,2]),1,nTimes.(lockNames{i}),1)); %[pp t cond tr]
        
        h = errorBarPlot(ssvepByFac,'area',1,'xaxisvalues',windows.(lockNames{i}));
        xline(0);
        xlabel(sprintf('time from %s (ms)', lockNames{i}));
        if i==1; ylabel('SSVEP amplitude'); end
        if i==1; legend([h{:,1}], labels.(factors{iF}),'Location','Best'); end
        title(factors{iF});
    end
end

%% split by acc/cond

accByCond = groupMeans(behData.acc,2,behData.cond,'dim');%[pp cond tr]

ssvepByCondAcc.resp = groupMeans(ssvepByCond.resp,4,repmat(permute(accByCond,[1,4,2,3]),1,nTimes.resp,1,1),'dim'); %[pp time cond acc tr]

figure();
h = errorBarPlot(reshape(nanmean(ssvepByCondAcc.resp,5),30,[],4),'area',1,'xaxisvalues',windows.resp);
xline(0);
xlabel(sprintf('time from %s (ms)', lockNames{2}));
ylabel('SSVEP amplitude'); 
legend([h{:,1}], col(strcat(repmat(labels.cond',1,2), repmat(labels.acc,2,1))),'Location','Best'); 


%% split that by CoM or conf6

% split factor by accByCond, split ssvep by that

iF = 4; 
factor = factors{iF};

facByCond = groupMeans(behData.(factor),2,behData.cond,'dim');
facByCondAcc = groupMeans(facByCond,3,accByCond,'dim');

ssvepByCondAccFac.resp = groupMeans(ssvepByCondAcc.resp,5,repmat(permute(facByCondAcc,[1,5,2,3,4]),1,nTimes.resp,1,1,1),'dim');
%[pp time cond acc fac tr]

figure();
for iC = 2 % cond
    for iA = 1:2
        subplot(1,2,iA);
        set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
        h = errorBarPlot(reshape(nanmean(ssvepByCondAccFac.resp(:,:,iC,iA,:,:),6),30,nTimes.resp,[]),'area',1,'xaxisvalues',windows.resp);
        xline(0); yline(0);
        xlim([0 1000]);
        xlabel(sprintf('time from %s (ms)', lockNames{2}));
        ylabel('SSVEP amplitude'); 
        legend([h{:,1}], labels.(factor),'Location','Best'); 
        
        title([labels.cond{iC}, ', ' labels.acc{iA}]);
    end
end

makeSubplotScalesEqual(1,2);

%% cond plot mean ssvep vs RTs
ssvepMeansByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), ssvepMeans,'UniformOutput',0);
behDataByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), behData,'UniformOutput',0);

% remove all interrupted?
for i = 1:2 %ignore preResp for now
    ssvepMeansByAcc.(lockNames{i})(behDataByAcc.cond==1) = NaN;
end

condNames = {'stim','resp','resp','resp';
    'RT','CoM','confRT','confInR1'};

% figure();
% 
% for i = 1:4
%     subplot(2,2,i)
%     [~,~,~,~,h] = conditionalPlot(permute(ssvepMeansByAcc.(condNames{1,i}),[3,1,2]), permute(behDataByAcc.(condNames{2,i}),[3,1,2]));
%     xlabel(sprintf('mean ssvep in %s', condNames{1,i}));
%     ylabel(condNames{2,i}); 
%     legend(h(2:2:end), labels.acc,'Location','Best'); 
% end

%% bar plot ssvep means by acc + factors - set resp-int to NaN

for i = 1:3 % time
figure();
[r,c] = GetSubPlotShape(nF);
for iF = 1:nF
    subplot(r,c,iF);
    set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
    respByAccFac = groupMeans(ssvepMeansByAcc.(lockNames{i}),3,behDataByAcc.(factors{iF})); % [pp acc fac]
    h = errorBarPlot(respByAccFac,'type','bar');
    xticks(1:2); xticklabels(labels.acc);
    ylabel(['SSVEP: ' lockNames{i}]);
    title(factors{iF});
end
SuperTitle(lockNames{i});
end

%% now by acc + fac - continued only
% same as above if excluding resp-int

% figure();
% [r,c] = GetSubPlotShape(nF);
% for iF = 1:nF
%     subplot(r,c,iF);
%     set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
%     respByAccFac = groupMeans(ssvepMeansByAcc.resp,3,behDataByAcc.(factors{iF})); % [pp cond fac]
%     h = errorBarPlot(respByAccFac,'type','bar');
%     xticks(1:2); xticklabels(labels.acc);
%     ylabel('SSVEP post-response');
%     title(factors{iF});
% end


%% acc*cert*CoM

ssvepMeansByAccCert = structfun(@(x) groupMeans(x,3,behDataByAcc.certainty,'dim'), ssvepMeansByAcc, 'UniformOutput',0);
behDataByAccCert = structfun(@(x) groupMeans(x, 3, behDataByAcc.certainty, 'dim'), behDataByAcc, 'UniformOutput',0);

ssvepMeansByAccCertCoM = structfun(@(x) groupMeans(x,4,behDataByAccCert.CoM,'dim'),ssvepMeansByAccCert, 'UniformOutput',0);
%[pp acc cert com tr]

iF=2;
figure();
for iV = 1:3
    for iA = 1:2
        subplot(3,2,(iV-1)*2+iA);
        set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');

        errorBarPlot(permute(nanmean(ssvepMeansByAccCertCoM.(lockNames{iV})(:,iA,:,:,:),5),[1,4,3,2,5]), 'type', 'bar');
        
        xticks(1:2); xticklabels(labels.CoM);
        xlabel('change of mind');
        if iA==1
            ylabel(['SSVEP :' lockNames{iV}]); 
            if iV==1; legend(labels.certainty); end
        end
        if iV==1; title(labels.acc(iA)); end
        
    end
end
        
    

    







