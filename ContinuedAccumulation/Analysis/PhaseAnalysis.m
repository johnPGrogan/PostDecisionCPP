% PhaseAnalysis
% Load GetPhaseTagging
% analyse like cpp

clc; clear; close all;

%% set stuff

useCSD = 1;
excludeBadPps = 1; % these are already excluded, as are isFlagged, during GetPhaseTagging.m
excludeTooFew = 1; % remove if <10 trials per cert*cond
excludeByRT = 1; % bad pps and flagged trials already gone
% bad pps are already excluded

outFolder = './Saves';
load('./Saves/ExtractEpochs.mat');

load('./Saves/BehDataLoad.mat','behData', 'labels','behTab');

load('./Saves/FlagArtefacts3.mat','isFlagged','ppsToExclude');

if useCSD
    loadName = fullfile(outFolder, 'GetPhaseTagging_CSD2.mat');
    topoName = fullfile(outFolder, 'GetPhaseTaggingTopo_CSD.mat');
else
    loadName = fullfile(outFolder, 'GetPhaseTagging.mat');
    topoName = fullfile(outFolder, 'GetPhaseTaggingTopo.mat');
end

load(loadName, 'phase','phaseTimes','nTimes');
phaseNames = {'resp'}; % can change this if discarding
nPhase = length(phaseNames); 

% set the post-response interrupted condition to NaN
% as there is no flicker, it can only be noise
phase.resp(repmat(permute(behData.cond==1,[1,3,2]),1,nTimes.resp)) = NaN;


%% other exclusions?
toRemove = sq(isnan(phase.resp(:,258,:))); % flagged trials or cpp outliers

if excludeBadPps
    toRemove(ppsToExclude,:) = 1;
end


if excludeByRT
    rtLims = [100 1540];
    toRemove(~isBetween(behData.RT, rtLims)) = 1;
end


if excludeTooFew
    % only remove those trials, not the entire person
    x = behData.certainty;
    x(behData.cond==1) = NaN; % no interrupteds
    x(toRemove)=NaN;
    toRemove2 = GetCellsWithLowN(10, x);

    toRemove(toRemove2) = 1;

%     % remove those people entirely, not just for that cond
%     toRemove(any(CountUnique(x,2)<10,2),:) = 1;


end


% exclude
for iPh = 1:nPhase
    phase.(phaseNames{iPh})(repmat(permute(toRemove,[1,3,2]),1,nTimes.(phaseNames{iPh})))  = NaN;
end

%% 


cmap = crameri('vik');

factors = {'acc','certainty','CoM','confInR1','conf3'};
nF = length(factors);
nF=3;
labels.diffNames = {'Error - Correct' ,'Low - High', 'CoM - NoCoM','Low - High','Low - High'};
diffFlips = [-1 -1 1 -1 -1]; % invert some to match order

cols = {[0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2];
        [flipud(crameri('roma',6)); .2 .2 .2]
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2]; };


% for taking means
windows.stim = [200 400];
windows.resp = [700 1000];
windows.preResp = [-300 -100];

% calc mean
for iPh = 1:nPhase
    meanPhase.(phaseNames{iPh}) = sq(nanmean(phase.(phaseNames{iPh})(:,isBetween(phaseTimes.(phaseNames{iPh}), windows.(phaseNames{iPh})),:),2));
end

% calc slopes
slopeWindow = [400 1000];
slopePhase.resp = FitCPPSlope(phase.resp, slopeWindow, phaseTimes.resp); % get slopes

%% set up regressions

% can have phase as a DV, or a covariate (if it improves fit)
behVarNames = {'acc','RT','confAcc','confRT','confInR1','certainty','CoM','confInCorr'};
nBehVars = length(behVarNames);

% remove trials
behTab(col(toRemove),:) = array2table(NaN(sum(toRemove,'all'),width(behTab)));

% put mean phases in
phaseTab = struct2table(structfun(@col, meanPhase, 'UniformOutput',0));

% combine
regTab = horzcat(behTab, phaseTab);
regTab.slopeResp = col(slopePhase.resp);
regNames = regTab.Properties.VariableNames; % get names

isLogistic = cellRegexpi(regNames, 'Logistic')>0;
regTab(:,~isLogistic) = varfun(@nanzscore, regTab(:,~isLogistic));
glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression, if behVars are DV

% remove nans?
toRemoveCol = true(height(regTab),1);
for iPh = 1:nPhase
    toRemoveCol = isnan(regTab.(phaseNames{iPh})) & toRemoveCol;
end
% regTab(toRemoveCol,:) = [];

%% stats in the paper

paperFormulae = {'resp ~ 1 + acc + (1 | pp)';
                 'resp ~ 1 + certainty + (1 | pp)';
                 'resp ~ 1 + CoM + (1 | pp)';
                 'resp ~ 1 + acc*CoM + (1 | pp)';
                 'resp ~ 1 + acc*certainty + (1 | pp)';
                 'resp ~ 1 + confInR1 + (1 | pp)';
                 'resp ~ 1 + confInR1*acc + (1 | pp)';};

paperFormulaeSep = {'resp ~ 1 + confInR1 + (1 | pp)';
                    'resp ~ 1 + certainty + (1 | pp)';
                    'resp ~ 1 + CoM + (1 | pp)';}; % to run in err/corr sep

fitglmeCell = @(x) fitglme(regTab, x);
% for eachcond sep, i=0 | 1
fitglmeCellSep = @(x, i) fitglme(regTab(regTab.acc>0 == i,:), x); 

paperFits = cellfun(fitglmeCell, paperFormulae,'UniformOutput',0);
paperSepFits = cellfun(fitglmeCellSep, repmat(paperFormulaeSep,1,2), repmat({0,1},3,1), 'UniformOutput',0);

%% save now

saveOpts = {'Volt','CSD'; '', 'ExclRT'};
saveName = sprintf('PhaseAnalysis_%s_%s.mat', saveOpts{1,useCSD+1}, saveOpts{2, excludeByRT+1});

save(fullfile(outFolder, saveName));

%%%%%%%%%%%%%%% don't run the figures below, just do PlotPaperFigures.m

%% split by things
% do all the splits here, for later use

% behData ones
behDataByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), behData, 'UniformOutput',0);
behDataByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),behData,'UniformOutput',0); %[pp cond tr]
behDataByCondAcc = structfun(@(x) groupMeans(x,3,behDataByCond.acc,'dim'), behDataByCond, 'UniformOutput',0);
behDataByAccCert = structfun(@(x) groupMeans(x, 3, behDataByAcc.certainty, 'dim'), behDataByAcc, 'UniformOutput',0);

% phase data ones
phaseByAcc = structfun(@(x) groupMeans(x, 3, repmat(permute(behData.acc,[1,3,2]),1,size(x,2)), 'dim'), phase, 'UniformOutput',0); %[pp t acc tr]
phaseByCond = structfun(@(x) groupMeans(x, 3, repmat(permute(behData.cond,[1,3,2]),1,size(x,2)), 'dim'), phase, 'UniformOutput',0); %[pp t cond tr]
phaseByCondAcc = structfun(@(x) groupMeans(x, 4, repmat(permute(behDataByCond.acc,[1,4,2,3]),1,size(x,2)), 'dim'), phaseByCond, 'UniformOutput',0); %[pp t cond tr]
phaseByAccCert = structfun(@(x) groupMeans(x,4, repmat(permute(behDataByAcc.certainty,[1,4,2,3]),1,size(x,2)),'dim'), phaseByAcc, 'UniformOutput',0); %[pp t acc cert tr]

% mean phase
meanPhaseByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), meanPhase, 'UniformOutput',0);
meanPhaseByCond = structfun(@(x) groupMeans(x,2,behData.cond,'dim'), meanPhase, 'UniformOutput',0);
meanPhaseByCondAcc = structfun(@(x) groupMeans(x,3,behDataByCond.acc,'dim'), meanPhaseByCond, 'UniformOutput',0);
meanPhaseByAccCert = structfun(@(x) groupMeans(x, 3, behDataByAcc.certainty, 'dim'), meanPhaseByAcc, 'UniformOutput',0);

slopePhaseByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), slopePhase, 'UniformOutput',0);

%% topo stuff - uncomment to run
% 
% load(topoName, 'projWavePreRespMeanWindows','projWaveRespMeanWindows','preRespWins','respWins')
% 
%  
% % remove flagged
% projWavePreRespMeanWindows(repmat(permute(toRemove,[1,3,4,2]),1,length(preRespWins),eeg.nChans,1)) = NaN;
% projWaveRespMeanWindows(repmat(permute(toRemove,[1,3,4,2]),1,length(respWins),eeg.nChans,1)) = NaN;
% 
% % set after int to NAN
% projWaveRespMeanWindows(repmat(permute(behData.cond==1,[1,3,4,2]),1,length(respWins),eeg.nChans,1)) = NaN;
% 
% 
% % % un-baseline?
% % if doBaseline==1
% %     projWavePreRespMeanWindows = exp(projWavePreRespMeanWindows) .* baseline;
% %     projWaveRespMeanWindows = exp(projWaveRespMeanWindows) .* baseline;
% % end   
% 
% 
% % animate topoplots
% AnimateTopoplots(sq(nanmean(projWavePreRespMeanWindows,[4 1]))', preRespWins, 1:length(preRespWins), eeg.chanlocs, [], [0 800], './Figs/PhaseTopoPreResp.gif');
% 
% % continued only
% AnimateTopoplots(sq(nanmean(projWaveRespMeanWindows,[4 1]))', respWins, 1:length(respWins), eeg.chanlocs, labels.cond(2), [0 800], './Figs/PhaseTopoRespCont.gif');
% 
% phaseTopoRespAcc = groupMeans(projWaveRespMeanWindows,4, repmat(permute(behData.acc,[1,3,4,2]),1,length(respWins),eeg.nChans, 1));
% % draw a few times only
% t = [100 400 700];
% figure();
% for i = 1:length(t)
%     for j = 1:2
%         subplot(2,length(t),i +(length(t) * (j-1)))
%         topoplot(nanmean(phaseTopoRespAcc(:,respWins==t(i),:,j),1),eeg.chanlocs, 'maplimits', [0 800]);
%         if j==1; title(t(i));end
%     end
% end
% c = colorbar('Location','NorthOutside');
% c.Position = [0.2893 0.4730 0.5036 0.0556];

%% grand mean

figure();
[r,c] = GetSubPlotShape(nPhase);

for iPh = 1:nPhase
    subplot(r,c,iPh);
    h = errorBarPlot(nanmean(phase.(phaseNames{iPh}),3), 'area',1, 'xaxisvalues', phaseTimes.(phaseNames{iPh}));
    xlabel(phaseNames{iPh});
    ylabel('phase tagging');
    xline(0);
    yline(0);

end


%% plot the traces for each condition


figure();
for iPh = 1:nPhase
    subplot(r,c,iPh);
    h = errorBarPlot(nanmean(phaseByCond.(phaseNames{iPh}),4), 'area',1, 'xaxisvalues', phaseTimes.(phaseNames{iPh}));
    xlabel(phaseNames{iPh});
    ylabel('phase tagging');
    xline(0);
    yline(0);
end
legend([h{:,1}], labels.cond, 'Location','Best');

%% now split by accuracy


figure();
for iPh = 1:nPhase
    subplot(r,c,iPh);
    h = errorBarPlot(reshape(nanmean(phaseByCondAcc.(phaseNames{iPh}),5),fileInfo.nPP,[],4), 'area',1, 'xaxisvalues', phaseTimes.(phaseNames{iPh}));
    xlabel(phaseNames{iPh});
    ylabel('phase tagging');
    xline(0);
    yline(0);

    if iPh==1; legend([h{:,1}], col(strcat(repmat(labels.cond',1,2), repmat(labels.acc,2,1))), 'Location','Best'); end
end

%% just fac

figure();

for iF = 1:nF

    for iPh = 1:nPhase
        subplot(nF,nPhase,(iF-1)*nPhase+iPh);
        set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
    
        data = groupMeans(phase.(phaseNames{iPh}),3, repmat(permute(behData.(factors{iF}),[1,3,2]),1,nTimes.(phaseNames{iPh}))); %[pp t f]

        h = errorBarPlot(data, 'area',1, 'xaxisvalues', phaseTimes.(phaseNames{iPh}));
        xlabel(phaseNames{iPh});
        ylabel('phase tagging');
        xline(0);
        yline(0);

        if iPh==1
            legend([h{:,1}], labels.(factors{iF}), 'Location','Best')
            ylabel('phase tagging');
        end
    end
end


%% cond and fac - not needed if excluding resp-interrupted
% 
% for iPh = 1:nPhase % epoch
%     figure();
%     for iF = 1:nF
%         data = groupMeans(phaseByCond.(phaseNames{iPh}),4, repmat(permute(behDataByCond.(factors{iF}),[1,4,2,3]),1,nTimes.(phaseNames{iPh})));
%         for iCond = 1:2
%             subplot(nF,2,(iF-1)*2+iCond);
%             set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
%         
%             h = errorBarPlot(sq(data(:,:,iCond,:)), 'area',1, 'xaxisvalues', phaseTimes.(phaseNames{iPh}));
%             xlabel(phaseNames{iPh});
%             xline(0);
%             yline(0);
%     
%             if iCond==1 && iPh~=3
%                 legend([h{:,1}], labels.(factors{iF}), 'Location','Best')
%                 ylabel('phase tagging');
%             end
%             if iF==1; title(labels.cond{iCond}); end
% 
%         end
%     end
%     SuperTitle(phaseNames{iPh});
% 
% end

%% stim+preResp - acc*fac 


for iPh = 1:nPhase
figure();

for iF = 1:nF-1 % ignore acc

    data = groupMeans(phaseByAcc.(phaseNames{iPh}), 4, repmat(permute(behDataByAcc.(factors{iF+1}),[1,4,2,3]),1,nTimes.(phaseNames{iPh}))); %[pp t ac fac]

    for iA = 1:2
        subplot(nF-1,2,(iF-1)*2+iA);
%         subplot(1,2,iA) % if using just one factor
        set(gca,'ColorOrder',cols{iF+1},'nextplot','replacechildren');
        h = errorBarPlot(sq(data(:,:,iA,:)), 'area',1, 'xaxisvalues', phaseTimes.(phaseNames{iPh}));
        xlabel(phaseNames{iPh});
        xline(0);
        yline(0);
    
        if iA==1
            legend([h{:,1}], labels.(factors{iF+1}), 'Location','Best')
            ylabel('phase tagging');
        end
        if iF==1; title(labels.acc{iA}); end
    end
end
% makeSubplotScalesEqual(2,2)
SuperTitle(phaseNames{iPh});
end

%% acc*fac for resp continued
% % (same as above if excluding resp-int)
% 
% figure();
% 
% iPh = 3;
% for iF = 1:nF-1 % ignore acc
% 
%     data = groupMeans(phaseByCondAcc.(phaseNames{iPh}), 5, repmat(permute(behDataByCondAcc.(factors{iF+1}),[1,5,2,3,4]),1,nTimes.(phaseNames{iPh}))); %[pp t c ac fac]
% 
%     for iA = 1:2
%         subplot(nF-1,2,(iF-1)*2+iA);
% %         subplot(1,2,iA) % if using just one factor
%         set(gca,'ColorOrder',cols{iF+1},'nextplot','replacechildren');
%         h = errorBarPlot(sq(data(:,:,2,iA,:)), 'area',1, 'xaxisvalues', phaseTimes.(phaseNames{iPh}));
%         xlabel(phaseNames{iPh});
%         xline(0);
%         yline(0);
%     
%         if iA==1
%             legend([h{:,1}], labels.(factors{iF+1}), 'Location','NorthWest')
%             ylabel('phase tagging');
%         end
%         if iF==1; title(labels.acc{iA}); end
%     end
% end
% SuperTitle(phaseNames{iPh});

%% cert for acc/com separately

iF = 2;

for iPh = 1:nPhase
figure();

data = groupMeans(phaseByAccCert.(phaseNames{iPh}), 5, repmat(permute(behDataByAccCert.CoM,[1,5,2,3,4]),1,nTimes.(phaseNames{iPh}))); %[pp t ac cert com]

for iCoM = 1:2
    for iA = 1:2
        subplot(nF-1,2,(iCoM-1)*2+iA);
%         subplot(1,2,iA) % if using just one factor
        set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
        h = errorBarPlot(sq(data(:,:,iA,:,iCoM)), 'area',1, 'xaxisvalues', phaseTimes.(phaseNames{iPh}));
        xlabel(phaseNames{iPh});
        xline(0);
        yline(0);
    
        if iA==1
            legend([h{:,1}], labels.(factors{iF}), 'Location','Best')
            ylabel('phase tagging');
        end
%         if iCoM==1; title(labels.acc{iA}); end
        title([labels.acc{iA}, ', ' labels.CoM{iCoM}]);
    end
end
makeSubplotScalesEqual(2,2)
SuperTitle(phaseNames{iPh});
end

%% plot means - acc*fac
figure();
for iPh=1:nPhase
%     figure(); % just resp
    [r,c] = GetSubPlotShape(nF);
    for iF = 1:nF
%         subplot(r,c,iF);
        subplot(nPhase,nF,(iPh-1)*nF+iF)
        set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
    
        data = groupMeans(meanPhaseByAcc.(phaseNames{iPh}), 3, behDataByAcc.(factors{iF})); %[pp a f]
        h = errorBarPlot(data,'type','bar');
        ylim(minMax([h.YData],'all') .* [.8 1.2]');
        xticks(1:2); xticklabels(labels.acc)
        xlim([.5 2.5]);
        title(factors{iF});
        ylabel([phaseNames{iPh} ' mean phase']);
        if iPh==1; legend(labels.(factors{iF}),'Location','Best'); end
    end
       
%     SuperTitle(phaseNames{iPh});
end

%% for resp, just do continued (same as above if excluding resp-int)
% iPh=3;
% figure(); % just resp
% iC = 2;
% for iF = 1:nF
%     
%     subplot(r,c,iF);
%     set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
% 
%     data = groupMeans(meanPhaseByCondAcc.(phaseNames{iPh}), 4, behDataByCondAcc.(factors{iF})); %[pp c a f]
%     h = errorBarPlot(permute(data(:,iC,:,:),[1,3,4,2]),'type','bar');
%     ylim(minMax([h.YData],'all') .* [.8 1.2]');
%     xticks(1:2); xticklabels(labels.acc)
%     xlim([.5 2.5]);
%     title(factors{iF});
%     ylabel('mean phase tagging');
%     legend(labels.(factors{iF}),'Location','Best');
%     
% end
%    
% %     makeSubplotScalesEqual(nF,2);
% SuperTitle(phaseNames{iPh});

%% conditional plots

% maybe first plot the means against RT?
% do some conditional plots

figure();
v = {'acc', 'RT'};
for iV = 1:length(v)
    for iPh = 1:nPhase
        subplot(length(v),nPhase,(iV-1)*nPhase+iPh);
        conditionalPlot(meanPhase.(phaseNames{iPh})', behData.(v{iV})');
        
        if iV==length(v); xlabel('mean phase'); end
        if iPh==1; ylabel(v{iV}); end
        if iV==1; title(phaseNames{iPh}); end

    end
    makeSubplotScalesEqual(length(v), nPhase, (1:nPhase) + (iV-1)*nPhase);
end

%% now split by acc

figure();
v = {'RT', 'confAcc', 'confRT', 'certainty', 'confInR1'};
for iV = 1:length(v)
    for iPh = 1:nPhase
        subplot(length(v),nPhase,(iV-1)*nPhase+iPh);
        conditionalPlot(permute(meanPhaseByAcc.(phaseNames{iPh}),[3,1,2]), permute(behDataByAcc.(v{iV}),[3,1,2]));
        
        if iV==length(v); xlabel('mean phase'); end
        if iPh==1; ylabel(v{iV}); end
        if iV==1; title(phaseNames{iPh}); end
        if iV==1 && iPh==1
            h = flip(findobj(gca,'Type','patch'));
            legend(h, labels.acc, 'Location','Best');
        end
        
    end
    makeSubplotScalesEqual(length(v), nPhase, (1:nPhase) + (iV-1)*nPhase);
end

%% for resp, just do cont (resp is same as above if interrupted excluded)
% % but will reduce power of stim+preResp by half
% 
% figure();
% v = {'RT', 'confAcc', 'confRT', 'certainty', 'confInR1'};
% [r,c] = GetSubPlotShape(length(v));
% for iV = 1:length(v)
%     y = permute(behDataByCondAcc.(v{iV})(:,2,:,:),[4,1,3,2]);
% 
%     for iPh = 3
%         subplot(r,c,iV);
%         x = permute(meanPhaseByCondAcc.(phaseNames{iPh})(:,2,:,:),[4,1,3,2]);
%         conditionalPlot(x, y);
%         
%         if iV==length(v); xlabel('mean phase'); end
%         if iPh; ylabel(v{iV}); end
%         if iV==1; title(phaseNames{iPh}); end
%         if iV==1
%             h = flip(findobj(gca,'Type','patch'));
%             legend(h, labels.acc, 'Location','Best');
%         end
%         
%         
%     end
% %     makeSubplotScalesEqual(length(v), nPhase, (1:nPhase) + (iV-1)*nPhase);
% end

%% acc*cert*CoM

figure();
for iPh = 1:nPhase
    data = groupMeans(meanPhaseByAccCert.(phaseNames{iPh}), 4, behDataByAccCert.CoM); % [pp acc cert com]
    for iA = 1:2
        subplot(nPhase,2,(iPh-1)*2+iA);
        set(gca,'ColorOrder',cols{2},'nextplot','replacechildren');

        errorBarPlot(permute(data(:,iA,:,:), [1,4,3,2]),'type','bar');

        xticks(1:2); xticklabels(labels.CoM);
%         xlabel('change of mind');
        if iA==1
            ylabel(['SSVEP :' phaseNames{iPh}]);
            if iPh==1; legend(labels.certainty); end
        end
        if iPh==1; title(labels.acc(iA)); end
        
    end
end
        
        





