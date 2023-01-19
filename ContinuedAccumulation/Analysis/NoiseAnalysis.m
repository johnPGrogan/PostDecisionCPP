% NoiseAnalysis
% Load up noise
% sum within each period
% split cpp etc by whether it is pos/neg

% RT vs abs(noise) neg correl is because it is total, so is bigger on slower RTs
% explains why is positive correl with slope - because dividing by larger time
% need to normalise by time if looking at that, or just take first 200ms etc?

%% 
clc; clear; close all;

outFolder = './Saves';

load('./Saves/NoiseData.mat');

load('./Saves/BehDataLoad.mat','behData','labels','fileInfo', 'behTab');

load('./Saves/FlagArtefacts3.mat', 'isFlagged','ppsToExclude');

doNorm = 1; % normalise noise by deltaC?
excludeInterrupted = 0; % remove interrupted to give similar number of pre/post?
recodeToChosen = 1; % make noise relative to finally chosen option
measToUse = 'sum'; % slope or sum

excludeBadPps = 1; % remove pps with <640 good trials?
excludeTooFew = 1; % remove pps with <20 per conf3
excludeByRT = 1; % remove trials outside [100 1500] ms
excludeCoMFromCert = 0; % remove CoM trials from behData.certainty

if strcmp(measToUse,'slope')
    winNames = {'preSlope','postSlope';'preSlopeStrong','postSlopeStrong'};
elseif strcmp(measToUse,'sum')
    winNames = {'preSum','postSum';'preSumStrong','postSumStrong'};
else
    error('measToUse must be sum or slope');
end


% need to remove bad trials/pps
behData.acc = double(behData.acc); % convert
behNames = fieldnames(behData);
nBeh = length(behNames);

% remove flagged trials and bad pps
for i = 1:nBeh
    
    behData.(behNames{i})(repmat(isFlagged',1,1,size(behData.(behNames{i}),3))) = NaN;
    behData.(behNames{i})(ppsToExclude,:,:) = NaN;

    % also remove bad RTs?
    rtLims = [100 1540];
    behData.(behNames{i})(repmat(~isBetween(behData.RT, rtLims), 1,1,size(behData.(behNames{i}),3))) = NaN;

    % remove interrupts
    if excludeInterrupted
        behData.(behNames{i})(repmat(behData.cond~=2,1,1,size(behData.(behNames{i}),3))) = NaN;
    end
    
end


labels.noise = {'weak','strong'};
labels.wins = {'pre','post'};
cols = [1 0 1; 0.4940 0.1840 0.5560];

toRemove = isFlagged'; % flagged trials

if excludeInterrupted
    % also remove interrupts from noiseData
    toRemove(behData.cond==1) = 1;
end

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
%     behData.acc(behData.CoM==1)=NaN;
end

if excludeTooFew
    % only remove those trials, not the entire person
    behData.certainty(toRemove) = NaN;
    toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.cond);

    toRemove(toRemove2) = 1;

end

% actually exclude now
noiseNames = fieldnames(noiseData);
noiseNames(ismember(noiseNames, noiseData.ignoreNames)) = []; % not the number of frames or sample times
for i = 1:length(noiseNames)
    noiseData.(noiseNames{i})(repmat(toRemove,1,1,size(noiseData.(noiseNames{i}),3))) = NaN;
end

varNames = {'acc','RT','confAcc','confRT','certainty','CoM','confInR1','confInCorr'};
nVars = length(varNames);
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


%% invert incorrect trials
if recodeToChosen
    toInvert = behData.acc == 0; % make relative to finally chosen option


    recodeNames = {'pre','post','preSum','postSum','respLocked','preSumStrong','postSumStrong','totalSum','totalSumStrong',...
        'respLockedCumul','preSlope','postSlope','totalSlope','preSlopeStrong','postSlopeStrong','totalSlopeStrong'};
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

%% 


factors = {'acc','certainty','CoM','confInR1','conf3'};
nF = length(factors);
labels.diffNames = {'Error - Correct' , 'Low - High', 'CoM - NoCoM','Low - High','Low - High' };
diffFlips = [-1 -1 1 -1 -1]; % invert some to match order

cols = {[0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2];
        [crameri('roma',6); .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2];};


%% regressions
if excludeCoMFromCert
    behTab.certainty(behTab.CoM>0) = NaN;
end


% use laoded behTab, exclude trials/pps
behTab(col(toRemove),:) = array2table(NaN(sum(toRemove,'all'),width(behTab)));

matNames2 = {'preSumStrong','postSumStrong','preSum','postSum','preSlopeStrong','postSlopeStrong','preSlope','postSlope','postSumEarly','postSumEarlyStrong','postSlopeEarly','postSlopeEarlyStrong'};
noiseNames = fieldnames(noiseData);
matData2 = struct2table(structfun(@col, rmfield(noiseData, noiseNames(~ismember(noiseNames, matNames2))), 'UniformOutput',0));


regTab = horzcat(behTab, matData2);
regNames = regTab.Properties.VariableNames;

% put abs noise in for RT?
regTab.absPreSum = col(abs(noiseData.preSum));
regTab.absPostSum = col(abs(noiseData.postSum));

isLogistic = cellRegexpi(regNames, 'Logistic')>0;
regTab(:,~isLogistic) = varfun(@nanzscore, regTab(:,~isLogistic));

glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression, if behVars are DV

% maybe delete all NaN rows?
regTab(all(isnan(table2array(regTab)),2),:) = [];

dvLogs = {'acc','accLogistic';'confAcc','confAccLogistic';'CoM','CoMLogistic'};

varNamesDV = varNames; % when these are DVs, use logistic ones
varNamesDV(ismember(varNames, dvLogs(:,1))) = dvLogs(:,2);

%% which of sum/sumstrong/slope/slopestrong gives best fits?

% ivOpts = {'preSum','postSum';'preSumStrong','postSumStrong';'preSlope','postSlope';'preSlopeStrong','postSlopeStrong'};
% 
% f1 = '%s ~ 1 + %s*%s*acc + (1 | pp)';
% for i = 1:4
%     fits1{i,1} = fitglme(regTab, sprintf(f1, 'CoM', ivOpts{i,1}, ivOpts{i,2}), 'link','logit','distribution','binomial');
%     bic(i) = fits1{i,1}.ModelCriterion.BIC;
% end
% [m,i] = min(bic)
% % seems that slope < sum << strongs for beh
% % sum < slope << strongs for cpp

% %% do behavioural regs
% 
% formulas = {'%s ~ 1 + preSum*postSum+ (1 | pp)';
%       '%s ~ 1 + preSum*postSum*acc + (1 | pp)';
%       '%s ~ 1 + preSum*postSum+ (1 | pp)';
%       '%s ~ 1 + preSum*postSum+ (1 | pp)'};
% nForms = length(formulas);
% inds = [~isnan(regTab.acc), ~isnan(regTab.acc), regTab.acc<=0, regTab.acc>=0]; % which trials to include
% 
% glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression
% 
% behRegs = cell(nForms, nVars);
% behRegTables = cell(nForms,1);
% for iF = 1:nForms
%     for iV = 1:length(varNames)
%         behRegs{iF,iV} = fitglme(regTab(inds(:,iF),:), sprintf(formulas{iF}, varNamesDV{iV}), glmeArgs{ isLogistic(strcmp(regNames, varNamesDV{iV})) + 1}{:});
%     end
% %     behRegTables{iF} = array2table(sq(permute(nancat(cellfun(@(x) x.Coefficients.pValue(2:end), behRegs(iF,:), 'UniformOutput',0)),[3,1,2]))','VariableNames',varNames,'RowNames',behRegs{iF,1}.CoefficientNames(2:end)); 
%     behRegTables{iF} = StatsTableFromLMECells(behRegs(iF,:), varNamesDV, 'pValue');
%     behRegTablesBeta{iF} = StatsTableFromLMECells(behRegs(iF,:), varNamesDV, 'Estimate');
% end


%% stats from paper supplement

paperFormulae = {'accLogistic ~ 1 + preSum + (1 | pp)';
                 'RTLog ~ 1 + preSum + (1 | pp)';
                 'RTLog ~ 1 + absPreSum + (1 | pp)';
                 'certainty ~ 1 + acc*preSum + (1 | pp)';
                 'confAcc ~ 1 + acc*preSum*postSum + (1 | pp)'; 
                 'certainty ~ 1 + acc*preSum*postSum + (1 | pp)';
                 'CoMLogistic ~ 1 + acc*postSum + (1 | pp)';
                 'confRTLog ~ 1 + absPostSum + acc*RTLog + (1 | pp)';
                 };

fitglmeCell = @(x) fitglme(regTab, x, glmeArgs{ ~isempty(regexp(x,'Logistic'))+ 1}{:});

paperFits = cellfun(fitglmeCell, paperFormulae, 'UniformOutput',0);

%% save


saveOpts = {'Raw','Norm'; '', 'ExclCoMFromCert'};
saveName = sprintf('NoiseAnalysis_%s_%s.mat', saveOpts{1,doNorm+1}, saveOpts{2, excludeCoMFromCert+1});

% save(fullfile(outFolder, saveName));


%% check noise

figure();

names = {'preSum','postSum','totalSum','preSlope','postSlope','totalSlope'};
for i = 1:length(names)
    subplot(2,3,i);
    histogram(col(noiseData.(names{i})),100);
    title(names{i});
end

% makeSubplotScalesEqual(2,3,1:3);
% makeSubplotScalesEqual(2,3,4:6);

%% are pre/post correlated?

figure();
subplot(2,2,1);
scatterRegress(col(noiseData.preSum), col(noiseData.postSum));
xlabel('pre resp noise sum');
ylabel('post resp noise sum');

subplot(2,2,2);
conditionalPlot(noiseData.preSum', noiseData.postSum');
xlabel('pre resp noise sum');
ylabel('post resp noise sum');

subplot(2,2,3);
scatterRegress(col(noiseData.preSlope), col(noiseData.postSlope));
xlabel('pre resp noise slope');
ylabel('post resp noise slope');

subplot(2,2,4);
conditionalPlot(noiseData.preSlope', noiseData.postSlope');
xlabel('pre resp noise slope');
ylabel('post resp noise slope');

%% correlate sum&slope

figure();
subplot(1,2,1);
scatterRegress(col(noiseData.preSum), col(noiseData.preSlope));
xlabel('sum');
ylabel('slope');
title('pre');


subplot(1,2,2);
scatterRegress(col(noiseData.postSum), col(noiseData.postSlope));
xlabel('sum');
ylabel('slope');
title('post');



%% show one person's data

inds = noiseData.times <= inf;
figure();
clear h;
plot(noiseData.times(inds), sq(noiseData.respLockedCumul(30,noiseData.totalSumStrong(30,:)==0,inds)), 'Color', cols{1}(1,:));
hold on;
plot(noiseData.times(inds), sq(noiseData.respLockedCumul(30,noiseData.totalSumStrong(30,:)==1,inds)), 'Color', cols{1}(2,:));
xlim([-1500 1000]);
xlabel('time to response (ms)');
ylabel('cumulative noise');
emptyLegend(2, { {'Color', cols{1}(2,:)}, {'Color',cols{1}(1,:)}}, {'LineWidth',2,'LineStyle','-'},flip(labels.noise), {'Location','NorthWest'});


%% show slope fitting to one trial

i=30;%j=2; % pp tr
js = find(~isnan(noiseData.postSum(i,:)),9,'first');
[r,c] = GetSubPlotShape(length(js));
figure();
for j = 1:length(js)
    subplot(r,c,j);
    plot(noiseData.times(noiseData.times<=0),sq(noiseData.respLockedCumul(i,js(j),noiseData.times<=0)),'x'); % pre
    hold on;
    plot(noiseData.times(noiseData.times>0),sq(noiseData.respLockedCumul(i,js(j),noiseData.times>0)),'x'); % post
    % also show cumsum of just post
    plot(noiseData.times(noiseData.times>0), sq(cumsum(noiseData.post(i,js(j),:),'omitnan')),'x')
    
    xlabel('time to resp (ms)');
    ylabel('cumulative noise');
    
    
    l = lsline;
end
makeSubplotScalesEqual(r,c);
slopes = diff(cat(1, l.YData),[],2) ./ diff(cat(1,l.XData),[],2);
disp(slopes);
disp([noiseData.postSlope(i,js(end)); noiseData.preSlope(i,js(end))]);

% %% how often is the evidence actually negative?
% 
% % 0 = deltaC, so -ve is less than deltaC 
% % <-deltaC is evidence towards the incorrect option
% % <-2*deltaC is suprathreshold evidence towards the incorrect
% c = [flipud(cols); 0 0 0; 1 0 0];
% 
% sumNames = {'preSum','postSum'; winNames{1,:}}; % 1st is what to split 2nd row by
% figure();
% for i = 1:2
%     subplot(2,2,i)
%     hold on;
%     splits = repmat(noiseData.(sumNames{2,i}),1,1,4);
%     splits(~cat(3, noiseData.(sumNames{1,i}) > 0, noiseData.(sumNames{1,i})<0 & noiseData.(sumNames{1,i})>=-thresh, noiseData.(sumNames{1,i}) < -thresh & noiseData.(sumNames{1,i})>=-2*thresh, noiseData.(sumNames{1,i})< -2*thresh)) = NaN;
%     for j = 1:4
%         plot([1:30]' + rand(30, 1280)/3, splits(:,:,j),'x','Color', c(j,:));
%     end
%     xlabel('participant');
%     ylabel('evidence')
%     xlim([0 31]);
%     title(sumNames{2,i});
%     emptyLegend(4, { {'Color', c(1,:)}, {'Color', c(2,:)}, {'Color', c(3,:)}, {'Color', c(4,:)}}, {'Marker','x','LineStyle','None'}, {'strong','weak','negative','supraNegative'},{'Location','Best'});
%     
%     subplot(2,2,i+2);gca
%     nPerThresh(:,:,i) = sq(sum(~isnan(splits),2));
%     nPerThresh(all(nPerThresh(:,:,i)==0,2),:,i) = NaN;
%     set(gca,'ColorOrder',c,'nextplot','replacechildren')
%     h = plot(nPerThresh(:,:,i),'x-');
%     xlabel('particpant');
%     ylabel('number of -ve evidence trials');
% end

%% figure 3.3


ylims = [870 940; .77 .82; 4.97 5.05; .12 .18; 570 630];
[r,c] = GetSubPlotShape(nVars);
anovas = cell(2,1);
for j = 1:2

    splitFactor = winNames{2,j};
    figure();
    for i = 1:length(varNames)
        subplot(r,c,i); 
        
        data = groupMeans(behData.(varNames{i}), 2, noiseData.(splitFactor)); % [pp cond]
        
        h = errorBarPlot(data,'type','bar');
        h.FaceColor = 'flat';
        h.CData = cols{1}(1:2,:);
        xticks(1:2);
        xticklabels({'weak','strong'});
    %     ylim(ylims(i,:));
        ylim(minMax([h.YData],2) .* [.95 1.05]);
        ylabel(varNames{i});
    
        anovas{j}.(varNames{i}) = rmanova(data, {'pp','noise'},'categorical',2,'DummyVarCoding','effects');
        if anovas{j}.(varNames{i}).pValue(2) < .05
            text(0.5, 0.9, '*','Units','Normalized','FontSize',16);
        end
    end
    SuperTitle(splitFactor);
    anovasTabs(j,:) = struct2table(structfun(@(x) x.pValue(2:end), anovas{j},'UniformOutput',0),'rowNames',{splitFactor});
end


%% no do that with acc as factor too
% accByNoise = groupMeans(behData.acc,2,noiseData.postStrong,'dim');

for j = 1:2
figure();
splitFactor = winNames{2,j};
for i = 1:nVars
    subplot(r,c,i); 
    set(gca,'ColorOrder',cols{1},'nextplot','replacechildren');

    % this only includes continued trials
    data = groupMeans(behData.(varNames{i}), 2, noiseData.(splitFactor)*2 + behData.acc); % [pp cond]
    data = reshape(data,[],2,2); %[pp acc noise 
    
    h = errorBarPlot(data,'type','bar');
    xticks(1:2);
    xticklabels(labels.acc);
%     ylim(ylims(i,:));
    ylim(minMax([h.YData],2) .* [.95 1.05]);
    ylabel(varNames{i});
    if i==1; legend(h, labels.noise, 'Location','Best'); end

    anovasAccPost.(varNames{i}) = rmanova(data, {'pp','acc','noise'},'categorical',2:3,'DummyVarCoding','effects');
    if anovasAccPost.(varNames{i}).pValue(end) < .05
        text(0.5, 0.9, '*','Units','Normalized','FontSize',16);
    end
end
SuperTitle(splitFactor);
anovasAccPostTab{j} = struct2table(structfun(@(x) x.pValue(2:end), anovasAccPost,'UniformOutput',0),'rowNames',anovasAccPost.acc.Term(2:end));
end
% %% pre & post as factors - nonexclusive. only includes continued trials
% 
% figure();
% splitFactor = winNames{2,2};
% for i = 1:nVars
%     subplot(r,c,i); 
%     set(gca,'ColorOrder',cols,'nextplot','replacechildren');
% 
%     data = groupMeans(behData.(varNames{i}), 2, noiseData.(winNames{2,1}) + noiseData.(winNames{2,2})*2, 'dim'); % [pp cond]
%     data = reshape(data,30,2,2,[]); %[pp acc noise
%     data = nancat(3, nanmean(data,[4 3]), sq(nanmean(data,[4 2]))); % average over each factor, collapsing across other
%     
%     h = errorBarPlot(data,'type','bar');
%     xticks(1:2);
%     xticklabels(strcat(labels.wins{1}, labels.noise));
% %     ylim(ylims(i,:));
%     ylim(minMax([h.YData],2) .* [.9 1.1]);
%     ylabel(varNames{i});
%     if i==1; legend(h, strcat(labels.wins{2}, labels.noise), 'Location','Best'); end
% 
%     anovasPrePost.(varNames{i}) = rmanova(data, {'pp','pre','post'},'categorical',2:3,'DummyVarCoding','effects');
%     if anovasPrePost.(varNames{i}).pValue(end) < .05
%         text(0.5, 0.9, '*','Units','Normalized','FontSize',16);
%     end
% end
% SuperTitle('post-decision noise');
% anovasPrePostTab = struct2table(structfun(@(x) x.pValue(2:end), anovasPrePost,'UniformOutput',0),'rowNames',anovasPrePost.acc.Term(2:end));

%% do exclusively, split by pre strength, then post strength (w-w, w-s, s-w, s-s)
figure();
splitFactor = winNames{2,2};
for i = 1:nVars
    subplot(r,c,i); 
    set(gca,'ColorOrder',cols{1},'nextplot','replacechildren');

    data = groupMeans(behData.(varNames{i}), 2, noiseData.(winNames{2,1}) + noiseData.(winNames{2,2})*2, 'dim'); % [pp cond]
    data = reshape(data,30,2,2,[]); %[pp acc noise
    data = nanmean(data,4); % only includes continued trials
    
    
    h = errorBarPlot(data,'type','bar');
    xticks(1:2);
    xticklabels(strcat(labels.wins{1}, labels.noise));
%     ylim(ylims(i,:));
    ylim(minMax([h.YData],2) .* [.9 1.1]);
    ylabel(varNames{i});
    if i==1; legend(h, strcat(labels.wins{2}, labels.noise), 'Location','Best'); end

    anovasPrePost.(varNames{i}) = rmanova(data, {'pp','pre','post'},'categorical',2:3,'DummyVarCoding','effects');
    if anovasPrePost.(varNames{i}).pValue(end) < .05
        text(0.5, 0.9, '*','Units','Normalized','FontSize',16);
    end
end
SuperTitle('post-decision noise');
anovasPrePostTab = struct2table(structfun(@(x) x.pValue(2:end), anovasPrePost,'UniformOutput',0),'rowNames',anovasPrePost.acc.Term(2:end));

%% acc+ pre & post as factors - exclusive. only continued

noisesByAcc = groupMeans(noiseData.(winNames{2,1}) + noiseData.(winNames{2,2})*2, 2, behData.acc, 'dim');
figure();
splitFactor = winNames{2,2};
for i = 1:nVars
    subplot(r,c,i); 
    set(gca,'ColorOrder', cols{1}(1:2,:) ,'nextplot','replacechildren');
    
    data = groupMeans(behData.(varNames{i}),2,behData.acc,'dim');
    data = groupMeans(data, 3, noisesByAcc, 'dim'); % [pp acc cond tr]
    data = reshape(data,30,2,2,2,[]); %[pp acc noise
%     data = nancat(3, nanmean(data,[5 4]), sq(nanmean(data,[5 3]))); % average over each factor, collapsing across other
    data = reshape(nanmean(data,5),30,2,4); % exclusivity
    
    h = errorBarPlot(data,'type','bar');
    [h(2:3).FaceAlpha] = deal(.5);
    xticks(1:2);
    xticklabels(labels.acc);
%     ylim(ylims(i,:));
    ylim(minMax([h.YData],2) .* [.9 1.1]);
    ylabel(varNames{i});
    if i==1; legend(h,[strcat(labels.noise{1},labels.noise), strcat(labels.noise{2},labels.noise)], 'Location','Best'); end

    anovasAccPrePost.(varNames{i}) = rmanova(reshape(data,30,2,2,2), {'pp','acc','preNoise','postNoise'},'categorical',2:4,'DummyVarCoding','effects');
    if anovasAccPrePost.(varNames{i}).pValue(end) < .05
        text(0.5, 0.9, '*','Units','Normalized','FontSize',16);
    end
end
SuperTitle('post-decision noise');
anovasAccPrePostTab = struct2table(structfun(@(x) x.pValue(2:end), anovasAccPrePost,'UniformOutput',0),'rowNames',anovasAccPrePost.acc.Term(2:end));

%% plot noise vs RT

figure();

iVs = find(ismember(varNames, {'confAcc','RT','confRT','CoM','confInR1'}));
if length(iVs)==1; iVs = [iVs iVs]; end
n = length(iVs);
x = cat(3, noiseData.(winNames{1,1})', noiseData.(winNames{1,2})');
for i = 1:n
    subplot(n,2,i*2-1);
    y = repmat(behData.(varNames{iVs(i)})',1,1,2);
    [~,~,~,~,h] = conditionalPlot(x, y);
    ylabel(varNames{iVs(i)});
    if i==1; legend(h(2:2:end), labels.wins, 'Location','Best'); end
        

    subplot(n,2,i*2);
    [~,~,~,~,h] = conditionalPlot(reshape(permute(groupMeans(x,1, repmat(behData.acc',1,1,2),'dim'),[4,2,1,3]),[],30,4),...
        reshape(permute(groupMeans(y,1, repmat(behData.acc',1,1,2),'dim'),[4,2,1,3]),[],30,4));
    if i==1; legend(h(2:2:end), col(strcat(repmat(labels.wins,2,1), repmat(labels.acc',1,2))), 'Location','Best'); end
    if i==n; xlabel(measToUse); end
end







%% as post-sum is NaN for interrupted, this is excluding all interrupted from the previous regressions
% run pre & post separately

inds = [1 1 2 2];
formulas = {'%s ~ 1 + %s*cond + (1 | pp)';
      '%s ~ 1 + %s*cond*acc + (1 | pp)';
      '%s ~ 1 + %s + (1 | pp)';
      '%s ~ 1 + %s*acc + (1 | pp)'};
nForms = length(formulas);

glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression

behRegs2 = cell(nForms, nVars);
behRegTables2 = cell(nForms,1);
for iF = 1:nForms
    for iV = 1:length(varNames)
        behRegs2{iF,iV} = fitglme(regTab, sprintf(formulas{iF}, varNamesDV{iV}, winNames{1,inds(iF)}), glmeArgs{ isLogistic(strcmp(regNames, varNamesDV{iV})) + 1}{:});
    end
    behRegTables2{iF} = StatsTableFromLMECells(behRegs2(iF,:), varNamesDV, 'pValue');
    behRegTablesBeta2{iF} = StatsTableFromLMECells(behRegs2(iF,:), varNamesDV, 'Estimate');
end

% names = {'pre','pre','post','post'};
figure();
for i = 1:nForms
    subplot(2,2,i);
    imagep(table2array(behRegTables2{i}),behRegTables2{i}.Properties.RowNames, behRegTables2{i}.Properties.VariableNames);
    title(winNames{1,inds(i)});
end


%% plot Betas for sigs

names = {'all','all','error','correct'};
figure();
for i = 1:nForms
    subplot(2,2,i);
    colormap(crameri('vik'));
    beta = table2array(behRegTablesBeta2{i});
    beta(table2array(behRegTables2{i}) >= .05) = 0;
%     clims =  [-1 1] .* max(abs(beta),[],'all');
%     if all(clims == 0), clims = [-.01 .01]; end
    clims = [-.05 .05];
    imagesc(beta,clims);
    colorbar;
    yticks(1:10); yticklabels(behRegTables2{i}.Properties.RowNames);
    xticks(1:12); xticklabels(behRegTables2{i}.Properties.VariableNames);
    title(winNames{1,inds(i)});
end

%% try running with noise as DV?

formulas = {'%s ~ 1 + %sacc + (1 | pp)';
      '%s ~ 1 + %sacc*certainty + (1 | pp)';
      '%s ~ 1 + %sacc*certainty*CoM + (1 | pp)';};
cond = {'cond*', []};

nForms = length(formulas);

glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression

behRegs3 = cell(nForms, 2);
behRegTables3 = cell(nForms,2);
for iF = 1:nForms
    for iV = 1:2 % pre/post
        behRegs3{iF,iV} = fitglme(regTab, sprintf(formulas{iF}, winNames{1,iV}, cond{iV}));
    end
end

names = {'pre','pre','post','post'};
figure();
for iF = 1:nForms
    for iV = 1:2
        subplot(nForms,2,(iF-1)*2+iV);
        imagep(behRegs3{iF,iV}.Coefficients.pValue,behRegs3{iF,iV}.CoefficientNames);
        if iF==1; title(winNames{1,iV}); end
    end
end


%% do pre in int/cond sep
formulas = {'%s ~ 1 + acc + (1 | pp)';
      '%s ~ 1 + acc*certainty + (1 | pp)';
      '%s ~ 1 + acc*certainty*CoM + (1 | pp)';};

names = {'all','int','cont'};
inds = [~isnan(regTab.preSum), regTab.cond<=0, regTab.cond>=0];

nForms = length(formulas);

glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression

behRegs4 = cell(nForms, 3);
behRegTables4 = cell(nForms,3);
for iF = 1:nForms
    for iV = 1:3 % all/int/cont
        behRegs4{iF,iV} = fitglme(regTab(inds(:,iV),:), sprintf(formulas{iF}, winNames{1,1}));
    end
    behRegTables4{iF} = StatsTableFromLMECells(behRegs4(iF,:), names, 'pValue');
    behRegTablesBeta4{iF} = StatsTableFromLMECells(behRegs4(iF,:), names, 'Estimate');

end


figure();
for i = 1:nForms
    subplot(2,2,i)
    imagep(table2array(behRegTables4{i}),behRegTables4{i}.Properties.RowNames, behRegTables4{i}.Properties.VariableNames);
end
SuperTitle(winNames{1,1});


%% plot average cumulative noise
% NoiseAveragesByFactors