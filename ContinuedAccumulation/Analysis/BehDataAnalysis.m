% BehDataAnalysis


%% load up data
clear;clc;

excludeByRT = 1;
excludeBadPps = 1;
excludeTooFew = 1;

load('./Saves/BehDataLoad.mat');

load('./Saves/FlagArtefacts3.mat', 'isFlagged','ppsToExclude');

%% filter + exclude!!!


toRemove = isFlagged';

if excludeBadPps
    toRemove(ppsToExclude,:) = 1;     
end

if excludeByRT
    rtLims = [100 1540];
    toRemove(~isBetween(behData.RT, rtLims)) = 1;
end

if excludeTooFew
    
    % only remove those trials, not the entire person
    behData.certainty(toRemove) = NaN;
    toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.cond);

    toRemove(toRemove2) = 1;
end

% exclude from each field
behNames = fieldnames(behData);
nBeh = length(behNames);

% remove flagged trials and bad pps
for i = 1:nBeh
    
    behData.(behNames{i})(repmat(toRemove,1,1,size(behData.(behNames{i}),3))) = NaN;
    
end


cols = {[0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410];
        crameri('roma',6);
        [0    0.4470    0.7410; 0.8500    0.3250    0.0980]};
labels.bothAcc = {'initial','second'};
%% histograms etc

% histogram confidence responses

figure();
names = {'acc','RT','confAcc','confRT','confResp','confInR1','conf3','CoM','scoretemp','certainty','confInCorr'};
[r,c] = GetSubPlotShape(length(names));
for i = 1:length(names)
    subplot(r,c,i);
    histogram(col(behData.(names{i})));
    xlabel(names{i});
end

%% hist per pp

figure();
[r,c] = GetSubPlotShape(length(names));

for i = 1:length(names)
    subplot(r,c,i);

    
    % lines for continuous, bars otherwise
    u = unique(behData.(names{i}));
    u(isnan(u)) = [];
    if length(u) > 7
        x = prctile(col(behData.(names{i})), 0:1:100)';
        y = [];
        for j = 1:30
            if any(~isnan(behData.(names{i})(j,:)))
                [y(:,j)] = ksdensity(behData.(names{i})(j,:),x);
            end
        end
        y(:,all(y==0,1)) = NaN;
        h = plot(x,y);
        for j = 1:30; h(j).Color(4) = .5; end

        hold on;
        plot(x, nanmean(y,2), '-k');
        if i==4; xlim([0 1000]); end
        
    else
        y = [];
        for j = 1:30
            if any(~isnan(behData.(names{i})(j,:)))
                y(:,j) = histcounts(behData.(names{i})(j,:),[-1;u] + .5);
            end
        end
        
        y(:,all(y==0,1)) = NaN;
        bar(u, nanmean(y,2));
        hold on;
        rands = rand(size(y)) .* .2 - .1;
        plot(u + rands, y, '-x');
    end
    xlabel(names{i});
end

%% count when splitting

% split by cond, acc, confAcc, conf3

% count how many trials of each conf3 for each
% per person

conf6ByCond = groupMeans(behData.confInR1,2,behData.cond,'dim');

conf6ByAcc =  groupMeans(behData.confInR1,2,behData.acc,'dim');

conf6ByCondAcc = groupMeans(conf6ByCond,3, groupMeans(behData.acc,2,behData.cond,'dim'),'dim');

% count # confs per cell (1st dim is conf level)
nConf6 = sum(permute(behData.confInR1,[3,1,2])==[1:6]',3);
nConf6ByCond = sum(permute(conf6ByCond,[4,1,2,3]) == [1:6]',4);
nConf6ByAcc = sum(permute(conf6ByAcc,[4,1,2,3]) == [1:6]',4);
nConf6ByCondAcc = sum(permute(conf6ByCondAcc,[5,1,2,3,4]) == [1:6]',5);

nConf6(:,ppsToExclude) = NaN;
nConf6ByCondAcc(:,ppsToExclude,:,:) = NaN;

% sum into conf3 + CoM
nConf3ByCondAcc = [ sum(nConf6ByCondAcc([1 2],:,:,:)); sum(nConf6ByCondAcc([3 4],:,:,:)); sum(nConf6ByCondAcc([5 6],:,:,:))];
nCoMByCondAcc = sq(sum(reshape(nConf6ByCondAcc,3,2,30,2,2),1));


%% remove pps with < 20 trials per conf3?

tooFewTrials = sq(any(sum(reshape(nConf6,2,3,[]))<20,2));


%% av freqs per pp, across cond/acc

% log them?
% nConf6=log(nConf6);
% nConf6(isinf(nConf6)) = NaN;

figure();
errorBarPlot(nConf6','type','line','plotargs',{'LineWidth',2,'Color','k'});
hold on;
h = plot(nConf6, '-');
for ii = 1:30
    h(ii).Color(4) = .5;
end
set(gca,'XTick', 1:6, 'XTickLabel', labels.confInR1);
xlim([.5 6.5]);
ylabel('# trials included');
% yticklabels(round(exp(yticks))); % if log



%% 

figure();
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        errorBarPlot(nConf6ByCondAcc(:,:,i,j)','type','line','plotargs',{'LineWidth',2,'Color','k'});
        hold on;
        h = plot(nConf6ByCondAcc(:,:,i,j), '-');
        for ii = 1:30
            h(ii).Color(4) = .5;
        end
        xlim([.5 6.5]);    
        title([labels.cond{i}, ', ', labels.acc{j}]);
        if j==1; ylabel('# trials'); end
        if i==2; xlabel('conf in 1st decision'); end
    end
end

%% conf3

figure();
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        errorBarPlot(nConf3ByCondAcc(:,:,i,j)','type','line','plotargs',{'LineWidth',2,'Color','k'});
        hold on;
        h = plot(nConf3ByCondAcc(:,:,i,j), '-');
        for ii = 1:30
            h(ii).Color(4) = .5;
        end
        xticks(1:3);
        xlim([.5 3.5]);
        title([labels.cond{i}, ', ', labels.acc{j}]);
        if j==1; ylabel('# trials'); end
        if i==2; xlabel('conf in 1st decision'); end
    end
end

%% CoM


figure();
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        errorBarPlot(nCoMByCondAcc(:,:,i,j)','type','line','plotargs',{'LineWidth',2,'Color','k'});
        hold on;
        h = plot(nCoMByCondAcc(:,:,i,j), '-');
        for ii = 1:30
            h(ii).Color(4) = .5;
        end
        xticks(1:2); xticklabels({'No Change', 'Change'});
        xlim([.5 2.5]);
        title([labels.cond{i}, ', ', labels.acc{j}]);
        if j==1; ylabel('# trials'); end
        if i==2; xlabel('change of mind?'); end
    end
end



%% make the exact same figs as wouter

figure();
ylims = [.77 .82; 4.7 5.1; 82 90; 0 .65; 0 .25; 200 700; 800 1050; 1 7];



subplot(2,5,1); % acc by cond
data = groupMeans(behData.acc, 2, behData.cond);
h = errorBarPlot(data, 'type','bar');
h.FaceColor = 'flat';
h.CData = cols{1};
ylim(ylims(1,:));
xticks(1:2); xticklabels(labels.cond);
ylabel('accuracy');

subplot(2,5,2); % conf by cond
data = groupMeans(behData.scoretemp, 2, behData.cond); % scoretemp
h = errorBarPlot(data, 'type','bar');
h.FaceColor = 'flat';
h.CData = cols{1};
ylim(ylims(2,:));
xticks(1:2); xticklabels(labels.cond);
ylabel('confidence');

subplot(2,5,3); % score by cond
data = groupMeans(dataMat(:,:,14), 2, behData.cond); % score
h = errorBarPlot(data, 'type','bar');
h.FaceColor = 'flat';
h.CData = cols{1};
ylim(ylims(3,:));
xticks(1:2); xticklabels(labels.cond);
ylabel('score');

subplot(2,5,4:5); % CoM by acc & conf3
data = groupMeans(behData.CoM, 2, behData.acc,'dim');
data = groupMeans(data, 3, groupMeans(behData.certainty,2,behData.acc,'dim')); % certainty
h = errorBarPlot(data, 'type','bar');
for i =1:3
    h(i).FaceColor = cols{2}(i,:);
end
ylim(ylims(4,:));
xticks(1:2); xticklabels(labels.acc);
legend(labels.conf3);
ylabel('CoM');

subplot(2,5,6); % CoM by cond
data = groupMeans(behData.CoM, 2, behData.cond);
h = errorBarPlot(data, 'type','bar');
h.FaceColor = 'flat';
h.CData = cols{1};
ylim(ylims(5,:));
xticks(1:2); xticklabels(labels.cond);
ylabel('CoM');

subplot(2,5,7); % confRT by cond
data = groupMeans(behData.confRT, 2, behData.cond);
h = errorBarPlot(data, 'type','bar');
h.FaceColor = 'flat';
h.CData = cols{1};
ylim(ylims(6,:));
xticks(1:2); xticklabels(labels.cond);
ylabel('confRT (ms)');

subplot(2,5,8); % RT by conf3
data = groupMeans(behData.RT, 2, behData.certainty); % certainty
h = errorBarPlot(data, 'type','bar');
h.FaceColor = 'flat';
h.CData = cols{2};
ylim(ylims(7,:));
xticks(1:2); xticklabels(labels.conf3);
ylabel('RT (ms)');

subplot(2,5,9:10); % conf by acc & cond
data = groupMeans(behData.scoretemp, 2, behData.acc,'dim');
data = groupMeans(data, 3, groupMeans(behData.cond,2,behData.acc,'dim')); % Wouter used 'behData.cond' here it seems
h = errorBarPlot(data, 'type','bar');
for i =1:2
    h(i).FaceColor = cols{1}(i,:);
end
ylim(ylims(8,:));
xticks(1:2); xticklabels(labels.acc);
legend(labels.cond);
ylabel('confidence');

%% speed-acc

saNames = {'acc','RT'; 'confAcc','confRT'};
toSplit = {'cond', 'certainty','conf3','confInR1'};
n = length(toSplit);

for j = 1:2
    figure();
    subplot(2,1,1);
    conditionalPlot(behData.(saNames{j,2})', behData.(saNames{j,1})');
    xlabel((saNames{j,2}));
    ylabel((saNames{j,1}));
    
    % split by factors
    for i = 1:n
        subplot(2,n,i+n);
        if strcmp(toSplit{i},'confInR1'); set(gca,'ColorOrder',crameri('roma',6),'nextplot','replacechildren');end
        [~,~,~,~,h] = conditionalPlot(permute(groupMeans(behData.(saNames{j,2}),2,behData.(toSplit{i}),'dim'),[3,1,2]), permute(groupMeans(behData.(saNames{j,1}),2,behData.(toSplit{i}),'dim'),[3,1,2]));
        xlabel((saNames{j,2}));
        ylabel((saNames{j,1}));
        legend(flip(h(2:2:end)), flip(labels.(toSplit{i})), 'Location','Best');
        title(toSplit{i});
        axis([500 1300 0.5*(j-1) 1]);
    end

end

%% split by cond + factors
saNames = {'acc','RT'; 'confAcc','confRT'};
toSplit = {'certainty','conf3','confInR1'};
n = length(toSplit);

for j = 2 % no point doing 1st decision
    figure();
%     subplot(2,1,1);
    x = permute(groupMeans(behData.(saNames{j,2}),2,behData.cond,'dim'),[3,1,2]);
    y = permute(groupMeans(behData.(saNames{j,1}),2,behData.cond,'dim'),[3,1,2]);
%     [~,~,~,~,h] = conditionalPlot(x,y);
%     legend(flip(h(2:2:end)), flip(labels.cond), 'Location','Best');
%     xlabel((saNames{j,2}));
%     ylabel((saNames{j,1}));
%     
    % split by factors
    for i = 1:n
        fac = permute(groupMeans(behData.(toSplit{i}),2,behData.cond,'dim'),[3,1,2]);
        x1 = permute(groupMeans(x,1,fac,'dim'),[4,2,1,3]);
        y1 = permute(groupMeans(y,1,fac,'dim'),[4,2,1,3]);

        for k = 1:2
            subplot(2,n,i+n*(k-1));
            if strcmp(toSplit{i},'confInR1'); set(gca,'ColorOrder',crameri('roma',6),'nextplot','replacechildren');end

            conditionalPlot(x1(:,:,:,k), y1(:,:,:,k));
            h = findobj(gca,'Type','Patch');
            xlabel((saNames{j,2}));
            ylabel((saNames{j,1}));
            legend((h), flip(labels.(toSplit{i})), 'Location','Best');
            title([labels.cond{k}, ' ' toSplit{i}]);
            axis([500 1300 0.5*(j-1) 1]);
        end
    end

end

%% is confAcc higher than acc, and does continued ev increase this?

figure();
set(gca,'ColorOrder',cols{4},'nextplot','replacechildren');

behDataByCond.bothAcc = nancat(3, groupMeans(behData.acc,2,behData.cond), groupMeans(behData.confAcc,2,behData.cond)); %[pp cond init/conf]
h = errorBarPlot(permute(behDataByCond.bothAcc,[1,3,2]), 'type','line','plotargs',{'LineWidth',2}); % permute to make cond the legend
set(gca,'XTick',1:2,'XTickLabel',labels.bothAcc);
ylabel('accuracy');
xlabel('response');
ylim([.7 .9]);
xlim([ .75 2.25]);
legend(h, labels.cond,'Location','Best');

condAnovas.bothAcc = rmanova(behDataByCond.bothAcc, {'pp','cond','respTime'},'categorical',2:3,'DummyVarCoding','effects');

% try lme
% use cond, not cond. pos = continued effect
regTableBoth = [ [col(behData.acc); col(behData.confAcc)], repmat(1:30,1,1280*2)', [ zeros(30*1280,1); ones(30*1280,1)], repmat(col(behData.cond),2,1)];
regTableBoth(:,2:end) = nanzscore(regTableBoth(:,2:end)); % leave acc as 1/0 for logistic reg
regTableBoth = array2table(regTableBoth, 'VariableNames',{'bothAcc','pp','respTime','cond'});

% logistic
fitsBothAcc = fitglme(regTableBoth, 'bothAcc ~ 1 + cond*respTime + (1 | pp)', 'link','logit','distribution','binomial');
% 


%% draw figs

% vars
varNames = {'acc','RT','confAcc','confRT','confInR1','certainty','CoM','confInCorr'};
nVars = length(varNames); % just these

% factors
factorNames = {'cond'};

% split by cond
figure();

[r,c] = GetSubPlotShape(nVars);
for i = 1:nVars
    subplot(r,c,i);
    behDataByCond.(varNames{i}) = groupMeans(behData.(varNames{i}),2,behData.(factorNames{1}),'dim');
    
    h = errorBarPlot(nanmean(behDataByCond.(varNames{i}),3), 'type','bar');
    h.FaceColor = 'flat';
    h.CData = cols{4};
    set(gca,'XTickLabel',labels.(factorNames{1}));
    
    % need to adjust ylims for some reason
    yl = ylim;
    if abs(diff(h.YData)) < (abs(diff(yl))/10)
        ylim(mean(h.YData) .* [.95 1.05]);
    end
    
    ylabel(varNames{i});
    condAnovas.(varNames{i}) = rmanova(nanmean(behDataByCond.(varNames{i}),3), {'pp',factorNames{1}},'categorical',2);

    if condAnovas.(varNames{i}).pValue(2) < .05
        text(0.5, 0.95, '*','Units','Normalized','FontSize',16);
    end
    
end

%% what is happening on CoM trials?
% hist confinr1 by acc/cond
figure();
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        histogram(behData.confInR1(behData.cond==i & behData.acc==(j-1) ));
        title([labels.cond{i} ', ' labels.acc{j}]);
    end
end
makeSubplotScalesEqual(2,2);

%% hist RT by factors
behDataByCond = structfun(@(x) groupMeans(x,2,behData.cond,'dim'), behData,'UniformOutput',0);

figure(); hold on;
for i = 1:2
    ksdensity(col(behDataByCond.confRT(:,i,:)),0:10:3000);
end
xlim([0 3000]);
legend(labels.cond)


%% and acc

behDataByCondAcc = structfun(@(x) groupMeans(x,3,behDataByCond.acc,'dim'), behDataByCond,'UniformOutput',0);

figure();
for iC = 1:2
    subplot(1,2,iC);
    set(gca,'ColorOrder', cols{1}, 'nextplot','replacechildren');
    hold on;
    for iA = 1:2
        ksdensity(col(behDataByCondAcc.confRT(:,iC,iA,:)),0:10:3000);
    end
    xlim([0 3000]);
    title(labels.cond(iC));
end
legend(labels.acc)

%% and cert

factor = 'confInR1';
iF=3; % for colour
confRTByCondAccFac = groupMeans(behDataByCondAcc.confRT,4,behDataByCondAcc.(factor),'dim');

figure();
for iC = 1:2
    for iA = 1:2
        subplot(2,2,(iA-1)*2+iC);
        set(gca,'ColorOrder', cols{iF}, 'nextplot','replacechildren');
        
        % get ksdensities and plot
        [f,x] = ksdensities(reshape(permute(confRTByCondAccFac(:,iC,iA,:,:),[1,2,3,5,4]),[],length(labels.(factor))), 0:10:2000);
        plot(x', f', 'LineWidth',2);

%         for j = 1:3
%             ksdensity(col(confRTByCondAccCert(:,iC,iA,j,:)),0:10:3000);
%         end
        xlim([0 2000]);
        title([labels.cond{iC} ', ' labels.acc{iA}]);
    end
end
legend(labels.(factor))

%% and CoM

confRTByCondAccCert = groupMeans(behDataByCondAcc.confRT,4,behDataByCondAcc.certainty,'dim');

CoMByCondAccCert = groupMeans(behDataByCondAcc.CoM,4,behDataByCondAcc.certainty,'dim');
confRTByCondAccCertCoM = groupMeans(confRTByCondAccCert,5,CoMByCondAccCert,'dim');
    
figure();

for iC = 1:2
    for iA = 1:2
        for iCoM = 1:2
            subplot(4,2,((iA-1)*2+iCoM-1)*2+iC)
            set(gca,'ColorOrder', cols{2}, 'nextplot','replacechildren');
            hold on;
            % get ksdensities and plot
            [f,x] = ksdensities(reshape(permute(confRTByCondAccCertCoM(:,iC,iA,:,iCoM,:),[1,2,3,5,6,4]),[],length(labels.(factor))), 0:10:2000);
            plot(x', f', 'LineWidth',2);
            xlim([0 2000]);
            title([labels.cond{iC} ', ' labels.acc{iA}, ', ' labels.CoM{iCoM}]);
        end
    end
end
legend(labels.certainty)
         

%% repeated measures anovas

% behTab loaded up
regTab = behTab; clear behTab;
% remove all NaN rows
regTab(all(isnan(table2array(regTab)),2),:) = [];
regNames = regTab.Properties.VariableNames;

%% regress each variable by cond?

lmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression

formulas = '%s ~ 1 + cond + (1 | pp)';

regTables = struct();
fits = cell(nVars,1);
for iV = 1:nVars
    fits{iV} = fitglme(regTab, sprintf(formulas, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:});
end
regTables.p = array2table(sq(nancat(cellfun(@(x) x.Coefficients.pValue(2:end), fits, 'UniformOutput',0)))','VariableNames',varNames,'RowNames',fits{1}.CoefficientNames(2:end)); 
regTables.beta = array2table(sq(nancat(cellfun(@(x) x.Coefficients.Estimate(2:end), fits, 'UniformOutput',0)))','VariableNames',varNames,'RowNames',fits{1}.CoefficientNames(2:end)); 


%%  interaction of cond*acc

formulas2 = '%s ~ 1 + cond*acc + (1 | pp)';

regTables2 = struct();
fits2 = cell(nVars,1);
for iV = 1:nVars
    fits2{iV} = fitglme(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:});
end
regTables2.p = array2table(sq(nancat(cellfun(@(x) x.Coefficients.pValue(2:end), fits2, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits2{1}.CoefficientNames(2:end)); 
regTables2.beta = array2table(sq(nancat(cellfun(@(x) x.Coefficients.Estimate(2:end), fits2, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits2{1}.CoefficientNames(2:end)); 

%% do these plots

% split condData by accuracy
[r,c] = GetSubPlotShape(nVars);

figure();
for i = 1:nVars
    subplot(r,c,i);
    set(gca,'ColorOrder',cols{4},'nextplot','replacechildren');
    condDataByAcc.(varNames{i}) = groupMeans(behDataByCond.(varNames{i}),3,behDataByCond.acc,'dim'); %[ pp cond acc tr]
    
    h = errorBarPlot(permute(nanmean(condDataByAcc.(varNames{i}),4),[1,3,2]), 'type','bar');
%     h.FaceColor = 'flat';
%     h.CData = cols{1};
    set(gca,'XTickLabel',labels.acc);
    
    % need to adjust ylims for some reason
    ylim( minMax([h.YData],2) .* [.9 1.1]);
    ylabel(varNames{i});
    
    if i==1; legend(h, labels.cond, 'Location','Best'); end
    
    
%     condAnovasByAcc.(varNames{i}) = rmanova(nanmean(condDataByAcc.(varNames{i}),4), {'pp','acc','cond'},'categorical',2:3);    
    if regTables2.p.(varNames{i})(3) < .05 % interaction
        text(0.5, 0.95, '*','Units','Normalized','FontSize',16);
    end
end

%% do above with cond on x-axis


% split condData by accuracy
figure();
for i = 1:nVars
    subplot(r,c,i);
    set(gca,'ColorOrder',cols{1},'nextplot','replacechildren');
    condDataByAcc.(varNames{i}) = groupMeans(behDataByCond.(varNames{i}),3,behDataByCond.acc,'dim'); %[ pp cond acc tr]
    
    h = errorBarPlot(permute(nanmean(condDataByAcc.(varNames{i}),4),[1,2,3]), 'type','bar');
%     h.FaceColor = 'flat';
%     h.CData = cols{1};
    set(gca,'XTickLabel',labels.cond);
    
    % need to adjust ylims for some reason
    ylim( minMax([h.YData],2) .* [.9 1.1]);
    ylabel(varNames{i});
    
    if i==1; legend(h, labels.acc, 'Location','Best'); end
    
    
%     condAnovasByAcc.(varNames{i}) = rmanova(nanmean(condDataByAcc.(varNames{i}),4), {'pp','acc','cond'},'categorical',2:3);    
    if regTables2.p.(varNames{i})(3) < .05 % interaction
        text(0.5, 0.95, '*','Units','Normalized','FontSize',16);
    end
end

%% regressions
% 
% bNames = [varNames factorNames 'corrLR'];
% regTab = rmfield(behData, behNames(~ismember(behNames, bNames))); % get without splits2
% regTab.pp = repmat([1:fileInfo.nPP]',1,fileInfo.maxTr);
% 
% regTab = structfun(@col, regTab, 'UniformOutput',0); % make into columns
% 
% 
% regTab = structfun(@nanzscore, regTab, 'UniformOutput',0); % zscore across all pps and conds
% 
% regTab = struct2table(regTab);
% 
% %% check normality?
% 
% figure();
% [r,c] = GetSubPlotShape(nVars);
% for iV = 1:nVars
%     subplot(r,c,iV);
%     hist(regTab.(varNames{iV}), 100);
%     title(varNames{iV});
% end
% 
% 
% 
% 
% %% do regressions
% 
% f = {'%s ~ 1 + acc*cond*conf3 + (1 | pp)';
%      '%s ~ 1 + acc*cond*CoM + (1 | pp)';
%      '%s ~ 1 + acc*cond*conf6 + (1 | pp)';
%      };
%  
% glmeOpts = { {'Distribution','binomial','Link','logit'};
%              {}; {};  
%              };
% 
% for iF = 1%:length(f)
%     for iV = 1:nVars
%         formula = sprintf(f{iF}, varNames{iV});
%         fits{iV, iF} = fitglme(regTab, formula);
%     end
% end
%      
% %% 
% 
% bic = cellfun(@(x) x.ModelCriterion.BIC, fits); % is which conf measure explains it best?
% [m,i] = min(bic,[],2);


%% paper figures

logHists = 1; % use log counts?

ylabels = { 'proportion correct', '# trials included', '# trials includede', ...
            'mean confidence (6pt scale)', 'mean confidence (3pt scale)', 'proportion of changes-of-mind'};

% counts per conf level
nConf6 = sum(permute(behData.confInR1,[3,1,2])==[1:6]',3);
conf6ByCond = groupMeans(behData.confInR1,2,behData.cond,'dim');
nConf6ByCond = sum(permute(conf6ByCond,[4,1,2,3]) == [1:6]',4);
nConf6(:,ppsToExclude) = NaN;
nConf6ByCond(:,ppsToExclude,:) = NaN;

behDataByCond = structfun(@(x) groupMeans(x,2,behData.cond,'dim'), behData, 'UniformOutput',0);

% combine init and conf acc
behDataByCond.bothAcc = nancat(3, groupMeans(behData.acc,2,behData.cond), groupMeans(behData.confAcc,2,behData.cond)); %[pp cond init/conf]

figure();

% 1 = acc: time*cond
subplot(2,3,1);

labels.bothAcc = {'initial','second'};
set(gca,'ColorOrder',cols{4},'nextplot','replacechildren');

h = errorBarPlot(permute(behDataByCond.bothAcc,[1,3,2]), 'type','line','plotargs',{'LineWidth',2}); % permute to make cond the legend
set(gca,'XTick',1:2,'XTickLabel',labels.bothAcc);
ylabel(ylabels(1));
xlabel('response');
ylim([.7 .9]);
xlim([ .75 2.25]);
legend(h, labels.cond,'Location','Best');
% title('accuracy');


% 2-3: histograms of conf6 by condiitons

if logHists % log them?
    nConf6ByCond2 =log(nConf6ByCond);
    nConf6ByCond2(isinf(nConf6ByCond2)) = 0; % remove -infs
else
    nConf6ByCond2 = nConf6ByCond;
end

for i = 1:2
    subplot(2,3,1+i);

    h = errorBarPlot(nConf6ByCond2(:,:,i)','type','bar');
    h.FaceAlpha = .5;
    h.FaceColor = 'flat';
    h.CData = [flipud(cols{2}); cols{2}];
    hold on;
    for ii = 1:6
        iCol = round(abs(ii-3.5));
        h1(ii) = plot(ii + rand(fileInfo.nPP,1)*.3-.15, nConf6ByCond2(ii,:,i), '.', 'MarkerSize', 8, 'Color', cols{2}(iCol,:));
    end

%     set(gca,'XTick', 1:6, 'XTickLabel', [fliplr(labels.certainty) labels.certainty], 'XTickLabelRotation', 45);
    set(gca,'XTick', 1:6, 'XTickLabel', []);
    xlabel(sprintf('%s       %s', labels.CoM{2}, labels.CoM{1}));
    xlim([0 7]);
    ylabel(ylabels(1+i));
    if logHists
        y = [0 5 10 20 50 100 200 300 400 500 600];
        set(gca,'YTick',[0 log(y(2:end))], 'YTickLabel', y);
    end% if log
    title(labels.cond{i});

    if i==1; legend(h1(4:6), labels.certainty,'Location','NorthWest'); end
end
makeSubplotScalesEqual(2,3,2:3);



varNames2 = {'confInR1', 'certainty', 'CoM'};
% 4-6: certainty: acc*cond, and Com + confInr1
for i = 1:3
    subplot(2,3,3+i);
    set(gca,'ColorOrder',cols{4},'nextplot','replacechildren');

    condDataByAcc.(varNames2{i}) = groupMeans(behDataByCond.(varNames2{i}),3,behDataByCond.acc,'dim'); %[ pp cond acc tr]
    h = errorBarPlot(permute(nanmean(condDataByAcc.(varNames2{i}),4),[1,3,2]), 'type','bar');
    
    set(gca,'XTickLabel',labels.acc);
    % need to adjust ylims for some reason
    ylim( minMax([h.YData],2) .* [.9 1.1]);
    ylabel(ylabels(3+i));
    xlabel('initial accuracy');
    
%     if i==1; legend(h, labels.cond, 'Location','Best'); end
%     title(varNames2{i});

end

% 4: certainty: acc*cond*CoM?


%% plot accuracy per conf response

varNames = {'acc','confAcc'};
names = {'initial','final'};
figure();
factor = 'confInR1';
for i = 1:2
    subplot(1,2,i);
    set(gca,'ColorOrder',cols{1},'nextplot','replacechildren'); % cond
    data = groupMeans(behDataByCond.(varNames{i}), 3,behDataByCond.confInR1);
    h = errorBarPlot(permute(data,[1 3 2]), 'type','line','plotargs',{'LineWidth',2});
    xticks(1:6);xticklabels(labels.confInR1);
    xlabel('confidnece in initial choice');
    ylabel(sprintf('%s accuracy',names{i}));
    
    title(sprintf('%s response',names{i}));
    
    [~,p(i,:)] = ttest(sq(data(:,1,:)), sq(data(:,2,:))); % paired t-tests
    hold on;
    pbar(p(i,:), 'yVal',1, 'alpha',.05, 'plotargs', {'Color','k','LineStyle','none','Marker','*','MarkerSize',10});

    if i==1; legend(h, labels.cond,'Location','Best','AutoUpdate','off'); end
    
end

%% ksdensity it
colours.cond = cols{1};
varNames = {'acc','confAcc'};
names = {'initial','final'};
figure();
factor = 'confInR1';
for i = 1:2
    data = groupMeans(behDataByCond.confInR1, 3,behDataByCond.(varNames{i}),'dim');
    c = CountUnique(data, 4); %[pp cond acc conf]
    c = c ./ sum(c,2); % normalise per pp
    
    subplot(1,2,i);
    set(gca,'ColorOrder',colours.cond(1:2,:),'nextplot','replacechildren', 'LineStyleOrder',{'--','-'});

    h = errorBarPlot(reshape(permute(c,[1,4,2,3]),30,6,4),'type','line','plotargs',{'LineWidth',2});
    xticks(1:6);xticklabels(labels.confInR1);
    xlabel('confidnece in initial choice');
    ylabel('% frequency');
    
    title(sprintf('%s response',names{i}));
    
%     [~,p(i,:)] = ttest(sq(data(:,1,:)), sq(data(:,2,:))); % paired t-tests
%     hold on;
%     pbar(p(i,:), 'yVal',1, 'alpha',.05, 'plotargs', {'Color','k','LineStyle','none','Marker','*','MarkerSize',10});
% 
    if i==1
     hold on; 
     emptyLegend(4, { {'Color', colours.cond(2,:),'LineStyle','-'},{'Color', colours.cond(1,:),'LineStyle','-'}, {'k','LineStyle','-'}, {'k','LineStyle','--'} }, {'LineWidth',2}, [flip(labels.cond) {'Initial Correct'},{'Initial Error'}], {'Location','Best'});

    end
end
