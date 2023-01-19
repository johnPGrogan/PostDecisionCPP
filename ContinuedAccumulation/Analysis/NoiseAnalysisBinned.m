% NoiseAnalysisBinned
% correl/regress binned noise slopes/ampl against cpp slopes/ampl
% in different time windows

clc; clear; close all;

%%

outFolder = './Saves';% where to load/save
useCSD = 1;

recodeToChosen = 1; % nosie relative to chosen option

doNorm = 1; % normalise noise by deltaC?
excludeInterrupted = 0; % remove interrupted to give similar number of pre/post?

excludeBadPps = 1; % remove pps with <640 good trials?
excludeTooFew = 1; % remove pps with <20 per conf3
excludeByRT = 1; % remove trials outside [100 1500] ms
excludeCoMFromCert = 0; % remove CoM trials from behData.certainty
doFilt=1; % use filtered cpp

measToUse = 'sum'; % slope or sum
cppVarName = 'sum'; % which cpp measure to use?
if strcmp(measToUse,'slope')
    winNames = {'preSlope','postSlopeEarly';'preSlopeStrong','postSlopeEarlyStrong'};
elseif strcmp(measToUse,'sum')
    winNames = {'preSum','postSumEarly';'preSumStrong','postSumEarlyStrong'};
else
    error('measToUse must be sum or slope');
end


load('./Saves/NoiseData.mat','noiseData');
% load('./Saves/BehDataLoad.mat','behData','labels','fileInfo','behTab');
load('./Saves/FlagArtefacts3.mat', 'isFlagged','ppsToExclude');

% load CPP from CPPAnalysis - already has pps excluded etc
saveOpts = {'Volt','CSD'; '', 'ExclCoMFromCert'};
saveName = sprintf('CPPAnalysis_%s_%s.mat', saveOpts{1,useCSD+1}, saveOpts{2, excludeCoMFromCert+1});
disp(saveName);

opts = struct('useCSD',useCSD, 'excludeBadPps', excludeBadPps, 'excludeTooFew',excludeTooFew,...
    'excludeByRT',excludeByRT, 'excludeCoMFromCert', excludeCoMFromCert,'doFilt',doFilt);
optsNames = fieldnames(opts);

data = load(fullfile(outFolder, saveName), optsNames{:}, ...
    'behData', 'cppFilt', 'factors','cols','nF', 'labels',...
    'eeg','amplWindows','regTab','cppVars');

% check they match
dataNames = fieldnames(data);
data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));

if ~equals(opts, data1)
    warning('loaded data and options do not match');
    keyboard;
end
struct2workspace(data);

labels.noise = {'weak','strong'};
labels.wins = {'pre','post'};
cols = [1 0 1; 0.4940 0.1840 0.5560];

colours.certainty = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];
colours.cond = [0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
colours.CoM = [0 0 1; 1 0 0; .2 .2 .2]; %[0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2]; % CoM
colours.conf3 = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];%colours.certainty;
colours.confInR1 =  [flipud(crameri('roma',6)); .2 .2 .2];
colours.confInCorr = colours.confInR1;
colours.cmap = crameri('vik');


%% need to remove bad trials/pps
% this is already done in CPPAnalysis, so just copy those

toRemove = sq(any(isnan(cppFilt),2)); % 

% % ppsToExclude = [2 16]; % can manually remove these?
% % ppsToExclude(16)=1; % noisy low numbers of low conf3
% if excludeBadPps
%     toRemove(ppsToExclude,:) = 1;     
% end
% 
% if excludeByRT
%     rtLims = [100 1540];
%     toRemove(~isBetween(behData.RT, rtLims)) = 1;
% end
% 
% if excludeCoMFromCert
%     %%%% set any certainty with CoM to NaN - will be ignored everywhere
%     behData.certainty(behData.CoM==1) = NaN; % will have to update in behTab
%     toRemove(behData.CoM==1) = 1;
% end
% 
% if excludeTooFew
%     
%     % only remove those trials, not the entire person
%     behData.certainty(toRemove) = NaN;
%     toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.cond);
% 
%     toRemove(toRemove2) = 1;
% % 
% %     % remove those people entirely, not just for that cond
% %     toRemove(any(CountUnique(groupMeans(behData.certainty,2,behData.cond,'dim'),3)<10,[2 3]),:) = 1;
%   
% end

if excludeInterrupted
    toRemove(behData.cond==1) = 1;
end

% actually remove from data now
cppFilt(repmat(permute(toRemove,[1,3,2]),1,size(cppFilt,2),1)) = NaN;
% cppStimiFilt(repmat(permute(toRemove,[1,3,2]),1,size(cppStimFilt,2),1)) = NaN;
% cppMeanTopo(repmat(toRemove,1,1,128,size(cppMeanTopo,4))) = NaN;

% also need to exclude from behData, as will use that to exclude low cell counts
behNames = fieldnames(behData);
for i = 1:length(behNames)
    behData.(behNames{i})(toRemove) = NaN;
end

% remove from noise too
noiseNames = fieldnames(noiseData);
noiseNames(ismember(noiseNames, noiseData.ignoreNames)) = []; % not the number of frames or sample times
for i = 1:length(noiseNames)
    noiseData.(noiseNames{i})(repmat(toRemove,1,1,size(noiseData.(noiseNames{i}),3))) = NaN;
end

%% shall I normalise them all by deltaC?
% could make it so 0 is equal, +1 is deltaC, -1 is incorrect option
% but then cumsum is always positive
% instead make 0 =deltaC, -1 = -deltaC
% 

if doNorm && noiseData.doNorm==0
    noiseNames = fieldnames(noiseData);
    noiseNames(ismember(noiseNames, noiseData.ignoreNames)) = []; % not the number of frames or sample times
    noiseNames(cellRegexpi(noiseNames, 'Strong')>0) = []; % ignore these too
    for i = 1:length(noiseNames)
        noiseData.(noiseNames{i}) = (noiseData.(noiseNames{i})) ./ noiseData.deltaC;
    end
    noiseData.thresh = 0;
else
    noiseData.thresh = 0;
end

%% recode relative to chosen/unchosen, not correct/incorrect
% i.e. positive=noise favours chosen, negative=favours unchosen

% invert incorrect trials
if recodeToChosen
    toInvert = behData.acc == 0; % make relative to finally chosen option


    recodeNames = {'pre','post','preSum','postSum','respLocked','preSumStrong','postSumStrong','totalSum','totalSumStrong',...
        'respLockedCumul','preSlope','postSlope','totalSlope','preSlopeStrong','postSlopeStrong','totalSlopeStrong',...
        'postEarly','postSumEarly','postSumEarlyStrong','postSlopeEarly','postSlopeEarlyStrong'};
    for i = 1:length(recodeNames)
        
        if regexp(recodeNames{i},'Strong') % should be 0/1 not -1/1
            noiseData.(recodeNames{i})(toInvert) = 1 - noiseData.(recodeNames{i})(toInvert); % invert double(logical)
        else
            noiseData.(recodeNames{i})(toInvert) = -noiseData.(recodeNames{i})(toInvert); % invert continuous
        end
    end
    
    noiseData.recodedToChosen = 1; % marker
    noiseData.ignoreNames{9} = 'recodedToChosen';
end


%% get summed noise within bins

bins = -1050:210:1050;
nBins = length(bins)-1;

noiseBinned = NaN(30,1280,nBins);
for i = 1:nBins
    winInds = isBetween(noiseData.times, bins([i i+1]));
    
    noiseBinned(:,:,i) = nanmean(noiseData.respLockedCumul(:,:,winInds),3);
end

%% use noise slopes?

noiseSlopeBinned = NaN(30,1280,nBins);
for i = 1:nBins
        
    noiseSlopeBinned(:,:,i) = FitCPPSlope(permute(noiseData.respLockedCumul,[1,3,2]), bins([i i+1]), noiseData.times);
end

%% which one to use?
if strcmp(measToUse, 'slope')
    noiseVarBinned = noiseSlopeBinned; % change this between slope and mean
else
    noiseVarBinned = noiseBinned;
end

%% plot

figure();

errorBarPlot(sq(nanmean(noiseVarBinned,2)), 'type','line');
xticks(1:nBins); xticklabels(bins(2:end));
xlabel('time to resp');
ylabel('noise');




%% get CPP slopes in bins

[slopeBinned, amplBinned] = deal(NaN(30,1280,nBins));

for i = 1:nBins
        
    slopeBinned(:,:,i) = FitCPPSlope(cppFilt, bins([i i+1]), eeg.respTimes);
    amplBinned(:,:,i) = sq(nanmean(cppFilt(:,isBetween(eeg.respTimes, bins([i i+1])),:),2));
end

if strcmp(cppVarName,'slope')
    cppVarBinned = slopeBinned;
else
    cppVarBinned = amplBinned;
end
%% plot mean slopes per bin

figure();

h = errorBarPlot(sq(nanmean(cppVarBinned,2)), 'type','line');
xticks(1:nBins); xticklabels(bins(2:end));
xlabel('time to resp');
ylabel('noise');


%% first look at correls

[rho, p] = deal(zeros(nBins));
for i=1:nBins
    for j = 1:nBins
        [rho(i,j), p(i,j)] = corr(col(noiseVarBinned(:,:,i)), col(cppVarBinned(:,:,j)),'rows','pairwise');
    end
end

%% plot

% remove nonsig?
rho(p>.05) = 0;
p(p>.05) = NaN;

figure();

subplot(1,2,1);
imagesc(rho');
colorbar;
xlabel(cppVarName);
ylabel('noise');
title('rho');
xticks(1:nBins);xticklabels(bins(2:end));
yticks(1:nBins);yticklabels(bins(2:end));

subplot(1,2,2);
imagesc(-log10(p)');
colorbar;
xlabel(cppVarName);
ylabel('noise');
title('-log10(p)');
xticks(1:nBins);xticklabels(bins(2:end));
yticks(1:nBins);yticklabels(bins(2:end));

%% do cond plots at each lag

figure();
[r,c] = GetSubPlotShape(nBins);
for i = 1:nBins
    subplot(r,c,i);

    conditionalPlot(noiseVarBinned(:,:,i)', cppVarBinned(:,:,find(bins==0)-1)');
    
    ylabel(cppVarName);
    xlabel(bins(i));
end

makeSubplotScalesEqual(r,c,1:nBins);
%% do regression inside each bin - control for 
% remove trials
% behTab(col(toRemove),:) = array2table(NaN(sum(toRemove,'all'),width(behTab)));


matVars = cat(2, reshape(cppVarBinned,[],nBins), reshape(noiseVarBinned,[],nBins));
regNames2 = strcat(col(repmat({cppVarName,'noise'},nBins,1))', col(repmat(cellfun(@num2str,num2cell(1:nBins),'UniformOutput',0),1,2))');
% matVars = nanzscore(matVars);


regTable = horzcat(regTab, array2table(matVars, 'VariableNames',regNames2));

% add in the total noise in each half - flip when responded right
regTable.noisePre = col(noiseData.preSum);
regTable.noisePost = col(noiseData.postSum);
regTable.noisePreSlope = col(noiseData.preSlope);
regTable.noisePostSlope = col(noiseData.postSlope);

regTabNames = regTable.Properties.VariableNames;

isLogistic = cellRegexpi(regNames2, 'Logistic')>0;
glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression


if strcmp(measToUse,'slope')
    IVs = {'noisePreSlope','noisePostSlope'};
else % sum
    IVs = {'noisePre','noisePost'};
end

% should I set interrupted postSum+slope to zero? allows me to include cond
% as a factor in pre&post
for i = 1:2
    regTable.(IVs{i})(regTable.cond==1) = 0;
end

%% normality

% figure();
% [r,c] = GetSubPlotShape(size(matVars,2));
% for i = 1:size(matVars,2)
%     subplot(r,c,i);
%     hist(regTable.(regNames2{i}),100);
%     title(regNames2{i});
%     xlim(nanstd(regTable.(regNames2{i})) .* [-4 4]);
% end


%% regress each noise bin to see which predicts slope at <0ms

% fitglme(regTable, 'slope5 ~ 1 + cond + acc + conf3 + CoM + (1 | pp)')
formula = '%s ~ 1 + %s (1 | pp)';

pVals = NaN(nBins);
bVals = zeros(nBins);
for i = 1:nBins
    noiseNames = strcat(regNames2((nBins+1):(nBins+i)),'+');
    f = fitglme(regTable, sprintf(formula, regNames2{i},[noiseNames{:}]));
    pVals(1:i,i) = f.Coefficients.pValue(2:end);
    bVals(1:i,i) = f.Coefficients.Estimate(2:end);
end

%% plot

% % remove nonsig
% pVals(pVals>.05) = NaN;
% bVals(isnan(pVals)) = 0;

figure();

subplot(1,2,1);
imagesc(bVals);
colorbar;
xlabel(cppVarName);
ylabel('noise');
title('beta');
xticks(1:nBins);xticklabels(bins(2:end));
yticks(1:nBins);yticklabels(bins(2:end));

subplot(1,2,2);
imagesc(-log10(pVals),[1.3010 5]);
colorbar;
xlabel(cppVarName);
ylabel('noise');
title('-log10(p)');
xticks(1:nBins);xticklabels(bins(2:end));
yticks(1:nBins);yticklabels(bins(2:end));

%% pre/post only

clear pVals2 bVals2
for i = 1:nBins
    f = fitglme(regTable, sprintf('%s ~ 1 + %s*%s*acc  + (1 | pp)', regNames2{i},IVs{1}, IVs{2}));
    pVals2(:,i) = f.Coefficients.pValue(2:end);
    bVals2(:,i) = f.Coefficients.Estimate(2:end);
end

%% plot

figure();

subplot(1,2,1);
imagesc(bVals2);
colorbar;
xlabel(cppVarName);
ylabel('factors');
title('beta');
xticks(1:nBins);xticklabels(bins(2:end));
yticks(1:nBins);yticklabels(f.CoefficientNames(2:end));

subplot(1,2,2);
pVals2(pVals2>.05)=NaN; % make min colour
pTicks = [1 .05 .01 .001 .0001 .00001];
imagesc(-log10(pVals2),-log10(flip(minMax(pTicks))));
c=colorbar;
c.Ticks = -log10(pTicks);
c.TickLabels = pTicks;
xlabel(cppVarName);
ylabel('factors');
title('-log10(p)');
xticks(1:nBins);xticklabels(bins(2:end));
yticks(1:nBins);yticklabels(f.CoefficientNames(2:end));


%% does total post-sum noise affect maybe trials?

cppVars.amplPostLate = sq(nanmean(cppFilt(:,isBetween(eeg.respTimes,[700 1000]),:),2));
cppVars.slopePostLate = FitCPPSlope(cppFilt, [700 1000], eeg.respTimes); % get slopes
% plot postsum noise vs cpp vars for each certainty
cppFactor = 'certainty';
noiseName = 'postSlopeEarly';

figure();
set(gca,'ColorOrder',colours.(cppFactor),'nextplot','replacechildren');
x = groupMeans(noiseData.(noiseName),2,behData.(cppFactor),'dim'); % [pp factor tr] (Extinguished already missing)
y = groupMeans(cppVars.amplPostLate,2,behData.(cppFactor),'dim');% [pp factor tr]

    
[~,~,~,~,h] = conditionalPlot(permute(x,[3,1,2]), permute(y,[3,1,2]));
ylabel('post-CPP mean');
xlabel(noiseName);
    
legend(h(2:2:end), labels.(cppFactor),'Location','Best');

%% use regression?

regTable.amplPostLate = nanzscore(col(cppVars.amplPostLate));
regTable.noisePostStrong = nanzscore(col(noiseData.postSumStrong));
regTable.noisePostSlopeStrong = nanzscore(col(noiseData.postSlopeStrong));
% noiseData.postSlopeSteep = abs(noiseData.postSlope) > nanmedian(abs(noiseData.postSlope(:))); % steep slope?
% regTable.noisePostSlopeSteep = nanzscore(col(noiseData.postSlopeSteep));
regTable.noisePostEarlyStrong = nanzscore(col(noiseData.postSumEarlyStrong));
regTable.noisePostSlopeEarlyStrong = nanzscore(col(noiseData.postSlopeEarlyStrong));



% fitglme(regTable(regTable.cond>=0 & regTable.certainty>0,:),'amplPostLate ~ 1 +noisePostSlopeStrong+ (1 | pp)') % certain
% fitglme(regTable(regTable.cond>=0 & regTable.certainty<0,:),'amplPostLate ~ 1 +noisePostSlopeStrong+ (1 | pp)') % maybe/prob

noiseVar = 'noisePostSlopeEarlyStrong'; % for regTAble
cppFactor = 'certainty';

fitglme(regTable(regTable.cond>=0,:),sprintf('amplPostLate ~ 1 + %s*%s+ (1 | pp)', noiseVar, cppFactor)) % maybe/prob
u=unique(regTable.(cppFactor)(~isnan(regTable.amplPost)));

clear f;
for i =1:length(u)
    f{i} = fitglme(regTable(regTable.cond>=0 & regTable.(cppFactor)==u(i),:), sprintf('amplPostLate ~ 1 + %s + (1 | pp)', noiseVar));
end
StatsTableFromLMECells(f,labels.(cppFactor), 'pValue')

% just in may/prob
fitglme(regTable(regTable.cond>=0 & regTable.certainty<0,:),sprintf('amplPostLate ~ 1 + %s*%s+ (1 | pp)', noiseVar,cppFactor)) % maybe/prob


%% split into pos/neg

noiseName = 'postSlopeEarlyStrong'; % >0

nT = size(cppFilt,2);
cppByNoise = groupMeans(cppFilt,3,repmat(permute(noiseData.(noiseName),[1,3,2]),1,nT),'dim'); %[pp t noise tr]

figure();
h = errorBarPlot(nanmean(cppByNoise,4),'area',1,'xaxisvalues',eeg.respTimes);
xlabel('time from init RT (ms)');
ylabel('CPP');
xlim([-300 1000]);
xline(0);yline(0);
legend([h{:,1}], labels.noise,'Location','Best');

%% and by fac
doStats = 1;
noiseVar = 'noisePostSlopeStrong'; % for regTAble

cppFactor = 'confInR1';
facByNoise = groupMeans(behData.(cppFactor),2,noiseData.(noiseName),'dim'); %[pp noise tr]
cppByNoiseFac = groupMeans(cppByNoise,4,repmat(permute(facByNoise,[1,4,2,3]),1,nT)); %[pp t noise fac]

% split into cert/conf3?
vecs.certainty = [3 2 1 1 2 3]; % cert
vecs.conf3 = [1 1 2 2 3 3];
% vecs.confInR1 = 1:6;
cppFactor2 = 'confInR1';
if ~strcmp(cppFactor2,cppFactor) && strcmp(cppFactor,'confInR1')
    cppByNoiseFac = groupMeans(cppByNoiseFac,4,repmat(permute(vecs.(cppFactor2),[1,3,4,2]),30,nT,2));
end

figure();
t=tiledlayout('flow');
times = -300:100:1000;
stats = NaN(1,length(times)-1,3);

for i = 1:length(labels.(cppFactor2))
    ax(i)=nexttile(t);
    h = errorBarPlot(cppByNoiseFac(:,:,:,i),'area',1,'xaxisvalues',eeg.respTimes);
    xlabel('time from init RT (ms)');
    ylabel('CPP');
    xlim([-300 1000]);
    xline(0);yline(0);

    if 0%doStats
        u=unique(regTable.(cppFactor2)(~isnan(regTable.amplPost)));
        formula = sprintf('amplWin ~ 1 + acc*%s + (1 | pp)', noiseVar); % formula to run, in this cond
        
        regTable.amplWin = nanzscore(col(cppFilt(:,1,:))); % use mean
        fit = fitglme(regTable(regTable.cond>=0 & regTable.(cppFactor2)==u(1),:), formula);
        ind = find(strcmp(fit.CoefficientNames,noiseVar));
        for iT = 1:length(times)-1           
            regTable.amplWin = nanzscore(col(nanmean(cppFilt(:,isBetween(eeg.respTimes, times([iT iT+1])),:),2))); % use mean            
            if sum(~isnan(regTable.amplWin)) > 100 % if at least 100 data
                fit = fitglme(regTable(regTable.cond>=0 & regTable.(cppFactor2)==u(i),:), formula);
                stats(:,iT,i) = fit.Coefficients.pValue(ind); % p-value
            end
        end
        hold on;
        yVal = min(ylim); 
        pbar(col([stats(1,:,i);stats(1,:,i)])', 'xVals',col([times(1:end-1);times(2:end)]), 'yVal',yVal,'plotargs',{'Color','k','LineWidth',3});
    end

    if i==1; legend([h{:,1}], labels.noise,'Location','Best'); end
    title(labels.(cppFactor2){i});
end
linkaxes(ax);

%fitglme(regTable(regTable.cond>=0,:), 'amplPost ~ 1 + certainty*noisePostStrong + (1 | pp)');

%% try using 0:500ms noise rather than whole post-dec time (rel to fin chosen still)

%% do rev-corr with noise rel to chosen

