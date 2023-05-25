function PlotNewFig4()
% Plot pre-resp vs pre-ev baselined Pe
% Requires:
% ContinuedAccumulation/Analysis/Saves/CPPAnalysis_CSD_.mat
gf = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',14);

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
    'behData', 'cppFilt', 'factors','cols','nF', 'labels',...
    'diffFlips','eeg','amplWindows',...
    'cppTopoAmplWindow','wins','cppMeanTopo','chans5','cmap');

data.diffFlips(2) = 1; % High - Low
data.labels.diffNames(2) = {'Certain - Maybe'};
data.labels.confInR1 = {'certain CoM', 'probably CoM', 'maybe CoM', 'maybe no-CoM', 'probably no-CoM', 'certain no-CoM'};
% % invert CoM so that change comes first
% data.behData.CoM = 2 - data.behData.CoM;
% data.labels.CoM = flip(data.labels.CoM);
data.labels.conf3 = {'certain/probably CoM', 'maybe CoM/no-CoM', 'probably/certain no-CoM'};
data.labels.certainty = {'maybe CoM/no-CoM', 'probably CoM/no-CoM', 'certain CoM/no-CoM'};
data.labels.baseline = {'pre-evidence baseline'};

%% check things match

dataNames = fieldnames(data);
data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));

if ~equals(opts, data1)
    warning('loaded data and options do not match');
    keyboard;
end

%%
data.amplWindows = [200 350; 700 1000]; % use B&Y2015 Pe window + our window
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
data.names.certainty = 'final-confidence';
data.names.conf3 = 'confidence-in-initial-choice: binned';
data.names.CoM = 'change-of-mind';
data.names.acc = 'initial accuracy';


%% split by ev-cond

data.behDataByCond = structfun(@(x) groupMeans(x, 2, data.behData.cond,'dim'),data.behData,'UniformOutput',0); %[pp cond tr]


%% split trace too (filtered)

data.cppFiltCond = groupMeans(data.cppFilt,3,repmat(permute(data.behData.cond,[1,3,2]),1,size(data.cppFilt,2)),'dim'); %[pp t cond tr]

data.nF = 3; % incl CoM

%% do pre-response baselined Pe too

respBaseline = [-100 0];

dataResp = data; % copy
baselineInds = isBetween(data.eeg.respTimes, respBaseline);
dataResp.cppFiltCond = data.cppFiltCond - nanmean(data.cppFiltCond(:,baselineInds,:,:),2);
dataResp.cppFilt = data.cppFilt - nanmean(data.cppFilt(:,baselineInds,:),2);
dataResp.cppMeanTopo = data.cppMeanTopo - nanmean(data.cppMeanTopo(:,:,:,isBetween(data.wins, respBaseline+1)),4);
% data2.cpp = data.cpp - nanmean(data.cpp(:,baselineInds,:),2);
yl2 = [-30 10]; % for means

dataResp.labels.baseline = {'pre-response baseline'};
% data.labels.baselines = {'pre-response'; 'pre-evidence'};
%% combine the means



%% plot

plotContinued = 0;

% % get un-filtered for means
% d2 = load(fullfile(opts.outFolder, opts.saveName), optsNames{:}, ...
%     'cpp');
% data.cpp = d2.cpp;

for iC = 1:(1+plotContinued)
    figure();
    dataResp.amplWindows = data.amplWindows; % reset for 2nd iter
    plotBy1CondFac(211, dataResp, 4, iC); % plot resp-BS traces, with both windows shown
%     plotBy1CondFac(222, data, 4, iC); % pre=ev

    % plot means sep for diff windows
    for iW = 1:2
        dataResp.amplWindows = data.amplWindows(iW,:);
        plotMeanBy1CondFac(222+iW, dataResp, 4, iC, yl2);
    end
    %     plotMeanBy1CondFac(224, data, 4, iC, yl2); % late window

end

%% topos

% split by cond
for iC = 1:(1+plotContinued)
    figure(); 
    for iW = 1:2
        dataResp.amplWindows = data.amplWindows(iW,:);
        plotToposCond(120+iW, dataResp, iC);
    end
%     plotToposCond(122, data, iC);
end

% plot CoM effect 
for iC = 1:(1+plotContinued)
    figure(); 
    for iW=1:2
        dataResp.amplWindows = data.amplWindows(iW,:);
        plotToposCondCoM(120+iW, dataResp, iC)
    end
%     plotToposCondCoM(122, data, iC)
end

% %% split by init acc too
% 
% for iC = 1:(1+plotContinued)
%     figure();
%     plotMeanBy1CondAccFac(121, dataResp, 4, iC, yl2)
%     plotMeanBy1CondAccFac(122, data, 4, iC, yl2)
% end


%%  
set(0,'DefaultAxesFontSize',gf); % reset

end



function plotBy1CondFac(subplotInds, data, iFs, iC)
struct2workspace(data);
if ~exist('iFs','var') || isempty(iFs)
    iFs = 2:3; % cert & CoM
end
doStats=1;
for i = 1:length(iFs)
    iF = iFs(i);

%     resplocked = groupMeans(cppFiltCond,4,repmat(permute(behDataByCond.(factors{iF}),[1,4,2,3]),1,size(cppFilt,2)),'dim'); %[pp t cond fac]
%     resplocked = reshape(permute(resplocked,[1 5 2 3 4]),[],size(cppFilt,2),2,size(resplocked,4)); % fold trials+pps into one, average across that - more similar to lme
    resplocked = groupMeans(cppFiltCond,4,repmat(permute(behDataByCond.(factors{iF}),[1,4,2,3]),1,size(cppFilt,2))); %[pp t cond fac]
%     resplocked(:,:,:,end+1) = diff(resplocked(:,:,:,[1 end]),[],4) * diffFlips(iF); % get diff

    if any(strcmpi(factors{iF},{'certainty','conf3','CoM'}))
        % average over confInR1 first
        resplocked = groupMeans(cppFiltCond,4,repmat(permute(behDataByCond.confInR1,[1,4,2,3]),1,size(cppFilt,2))); %[pp t cond confInR1]
        facByConfInR1 = groupMeans(behDataByCond.(factors{iF}),3,behDataByCond.confInR1); %[pp cond confInR1]
        resplocked = groupMeans(resplocked,4,repmat(permute(facByConfInR1,[1,4,2,3]),1,size(cppFilt,2))); %[pp t cond fac]
    end

%     for iC=1
        subplot(subplotInds(i));
        set(gca,'ColorOrder',colours.(factors{iF}),'nextplot','replacechildren');

        h = errorBarPlot(sq(resplocked(:,rInds,iC,:)),'area',1,'xaxisvalues',eeg.respTimes(rInds));
        if ~showErrBars
            for j=1:size(h,1) 
                h{j,2}.Visible='off';
                h{j,1}.LineWidth = 2;
            end % remove error bars
        else % make more transparent
            for j=1:size(h,1)
                h{j,2}.FaceAlpha = alphas.(factors{iF});
                h{j,1}.LineWidth = 2;
            end 
        end
        for j=1:size(h,1)
            h{j,1}.LineWidth = data.lines.widths.(factors{iF})(j);
            h{j,1}.LineStyle = data.lines.styles.(factors{iF}){j};
        end

        % shade the mean window
        hold on;
        for iB = 1:size(amplWindows,1)
            fill(col(repmat(amplWindows(iB,:),2,1)), [yl(1,:) fliplr(yl(1,:))], 'k','FaceAlpha',faceAlpha,'LineStyle','none');
        end


        xticks(-1000:500:1000); xtickangle(0);
        xline(0);yline(0);
        xlim(respWin); %ylim(ylims(2,:));
        xlabel('time from initial RT (ms)');

        if iC==1; ylabel('Pe \muV/m^2'); end
        ylim(yl(1,:));

        if i==1; title([labels.cond(iC); labels.baseline]);
        else; title(names.(factors{iF})); end
        

        if doStats % regress factor within this cond, at every 100ms bin

            if ~exist('regTab','var')
                regTab = table;
                regTab.pp = nanzscore(col(behData.pp)); % participant for RE
                regTab.cond = col(behData.cond);
            end
            regTab.(factors{iF}) = nanzscore(col(behData.(factors{iF}))); % add this factor
            times = 0:100:1000;
%             x = movmean(times,2,2,'endpoints','discard');
            formula = sprintf('amplWin ~ 1 + %s + (1 | pp)', factors{iF}); % formula to run, in this cond
            
            % run one reg to get number of terms
            regTab.amplWin = nanzscore(col(cppFilt(:,1,:))); % use mean
            fit = fitglme(regTab(regTab.cond==iC,:), formula);
            stats = NaN(length(fit.CoefficientNames)-1,length(times)-1,2);

            for iT = 1:length(times)-1
            
                % take mean within window
                regTab.amplWin = nanzscore(col(nanmean(cppFilt(:,isBetween(eeg.respTimes, times([iT iT+1])),:),2))); % use mean
            
            
                if sum(~isnan(regTab.amplWin)) > 100 % if at least 100 data
                    fit = fitglme(regTab(regTab.cond==iC,:), formula);
%                     stats(:,iT,1) = fit.Coefficients.Estimate(2:end); % beta
                    stats(:,iT,2) = fit.Coefficients.pValue(2:end); % p-value
                end
            end

            % plot sig timepoints
%             plot(, stats(:,:,2) < alpha, '.k', 'MarkerSize',20)
            hold on;
            yVal = min(ylim); 
            % duplicate to get it from start:end of a single window if
            % needed
            pbar(col([stats(1,:,2);stats(1,:,2)])', 'xVals',col([times(1:end-1);times(2:end)]), 'yVal',yVal,'plotargs',{'Color','k','LineWidth',3});
%             h1.LineStyle = 'none'; h1.Marker = '.'; h1.MarkerSize = 20;

        end

        if mod(subplotInds(i),2)
            legend([h{:,1}], [labels.(factors{iF}) labels.diffNames{iF}], 'Location','Best', 'FontSize',12,'Box','off');
        end

%     end



end
end

function plotMeanBy1CondFac(subplotInds, data, iFs, iC, ylims)
struct2workspace(data);
if ~exist('iFs','var') || isempty(iFs)
    iFs = 2:3; % cert & CoM
end

isEvBS = strcmp(labels.baseline, 'pre-evidence baseline');

doStats=1;
for i = 1:length(iFs)
    iF = iFs(i);

    amplByCond = sq(nanmean(cppFiltCond(:,isBetween(eeg.respTimes,amplWindows),:,:),2)); %[pp cond tr]
    amplByCondFac = groupMeans(amplByCond,3,behDataByCond.(factors{iF})); %[pp cond fac]

    if any(strcmpi(factors{iF},{'certainty','conf3','CoM'}))
        % average over confInR1 first
        amplByCondFac = groupMeans(amplByCond,3,behDataByCond.confInR1); %[pp cond confInR1]
        facByConfInR1 = groupMeans(behDataByCond.(factors{iF}),3,behDataByCond.confInR1); %[pp cond confInR1]
        amplByCondFac = groupMeans(amplByCondFac, 3, facByConfInR1); %[pp cond fac]
    end


    subplot(subplotInds(i));
    set(gca,'ColorOrder',colours.(factors{iF}),'nextplot','replacechildren');

    h = errorBarPlot(sq(amplByCondFac(:,iC,:)),'type','bar');
    h.FaceColor = 'flat';
    h.CData = colours.(factors{iF})(1:size(amplByCondFac,3),:);
        
    xticks(3.5); xticklabels(names.(factors{iF})); xtickangle(0);
%     xlabel('Baseline period');
%     xticks(1:6); xticklabels(labels.(factors{iF})); xtickangle(90);
%     xlabel(names.(factors{iF}));
%     title(sprintf('%s baseline', labels.baseline{1}));

%     ylabel(sprintf('mean amplitude %d:%d ms', amplWindows(2,1), amplWindows(2,2)));
    ylabel('mean amplitude');
    title(sprintf('%d: %d ms', amplWindows(1), amplWindows(2)));
    ylim(ylims); 

    if doStats
        if ~exist('regTab','var')
            regTab = table;
            regTab.pp = nanzscore(col(behDataByCond.pp(:,iC,:))); % participant for RE
        end
        regTab.(factors{iF}) = nanzscore(col(behDataByCond.(factors{iF})(:,iC,:))); % add this factor

        
            
            newAmpl = sq(nanmean(cppFiltCond(:, isBetween(eeg.respTimes, amplWindows),iC,:),2)); %[pp tr]
            regTab.f = nanzscore(col(newAmpl));
            f = fitglme(regTab, sprintf('f ~ 1 + %s + (1 | pp)', factors{iF}));
            if f.Coefficients.pValue(2)<.05
                hold on;
                text(3.5, max(ylim)*.95, '*','FontSize',20);
            end
        
    end
    
    if mod(subplotInds(i),2)
        for j = 1:length(labels.(factors{iF}))
            c{1,j} = {'FaceColor', colours.(factors{iF})(j,:)};
        end
        hold on;
        emptyLegend(j, c, {}, labels.(factors{iF}), {'Location','Best', 'Box','off'}, @bar)
    end     

end
end

function plotMeanBy1CondAccFac(subplotInds, data, iFs, iC, ylims)
struct2workspace(data);
if ~exist('iFs','var') || isempty(iFs)
    iFs = 2:3; % cert & CoM
end

isEvBS = strcmp(labels.baseline, 'pre-evidence baseline');

doStats=1;
for i = 1:length(iFs)
    iF = iFs(i);

    amplByCond = sq(nanmean(cppFiltCond(:,isBetween(eeg.respTimes,amplWindows(2,:)),:,:),2)); %[pp cond tr]
    amplByCondFac = groupMeans(amplByCond,3,behDataByCond.(factors{iF}) + behDataByCond.acc/10); %[pp cond acfac]
    amplByCondAccFac = reshape(amplByCondFac,30,2,2,length(labels.(factors{iF}))); %[pp cond acc fac]

    if any(strcmpi(factors{iF},{'certainty','conf3','CoM'}))
        % average over confInR1 first
        amplByCondFac = reshape(groupMeans(amplByCond,3,behDataByCond.confInR1 + behDataByCond.acc/10),30,2,2,length(labels.(factors{iF}))); %[pp cond acc+confInR1]
        facByConfInR1 = reshape(groupMeans(behDataByCond.(factors{iF}),3,behDataByCond.confInR1 + behDataByCond/10), 30, 2,2,length(labels.(factors{iF}))); %[pp cond acc+confInR1]
        amplByCondAccFac = groupMeans(amplByCondFac, 3, facByConfInR1); %[pp cond acc fac]
    end


    subplot(subplotInds(i));
    set(gca,'ColorOrder',colours.(factors{iF}),'nextplot','replacechildren');

    h = errorBarPlot(sq(amplByCondAccFac(:,iC,:,:)),'type','bar');
%     h.FaceColor = 'flat';
%     h.CData = colours.(factors{iF})(1:size(amplByCondFac,3),:);
        
    xticks(1:2); xticklabels(labels.acc);xtickangle(0);
    xlabel('initial accuracy');
    ylabel('mean amplitude');
    ylim(ylims); 

    if doStats
        if ~exist('regTab','var')
            regTab = table;
            regTab.pp = nanzscore(col(behDataByCond.pp(:,iC,:))); % participant for RE
            regTab.acc = (col(behDataByCond.acc(:,iC,:))); % participant for RE
        end
        regTab.(factors{iF}) = nanzscore(col(behDataByCond.(factors{iF})(:,iC,:))); % add this factor

        
        newAmpl = sq(nanmean(cppFiltCond(:, isBetween(eeg.respTimes, amplWindows(2,:)),iC,:),2)); %[pp tr]
        regTab.f = nanzscore(col(newAmpl));
        for iA = 1:2
            f = fitglme(regTab(regTab.acc==(iA-1),:), sprintf('f ~ 1 + %s + (1 | pp)', factors{iF}));
            if f.Coefficients.pValue(2)<.05
                hold on;
                text(iA, max(ylim)*.95, '*','FontSize',20);
            end
        end
        % interaction too?
        f = fitglme(regTab, sprintf('f ~ 1 + acc*%s + (1 | pp)', factors{iF}));
        if f.Coefficients.pValue(end)<.05
            hold on;
            text(1.5, max(ylim)*.95, 'r', '*','FontSize',20);
        end
    end
    
    if mod(subplotInds(i),2)
%         for j = 1:length(labels.(factors{iF}))
%             c{1,j} = {'FaceColor', colours.(factors{iF})(j,:)};
%         end
%         hold on;
%         emptyLegend(j, c, {}, labels.(factors{iF}), {'Location','Best', 'Box','off'}, @bar)
        legend(h, labels.(factors{iF}), 'Location','Best', 'FontSize',12,'Box','off');

    end     

end
end

function plotToposCond(subplotInds, data, iC)
struct2workspace(data);

mapLims = [-15 15];

cppMeanTopoByCond = groupMeans(cppMeanTopo,2,repmat(behData.cond,1,1,eeg.nChans,length(wins)-1)); %[pp cond chan time]

% put chans5 electrodes on?
[~, chans5Inds] = ismember(chans5, eeg.chanNames);

subplot(subplotInds);
winInds = isBetween(wins, amplWindows+1); % add 1, as bin is from >min to >=max
topoplot(sq(nanmean(nanmean(cppMeanTopoByCond(:,iC,:,winInds(2:end)),4),1)),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'mapLimits',mapLims,...
    'emarker2',{chans5Inds,'.','k',10,1});
title([labels.cond(iC); labels.baseline; {sprintf('%d:%d ms', amplWindows(1), amplWindows(2))}]);

c = colorbar;

end


function plotToposCondCoM(subplotInds, data, iC)
struct2workspace(data);

mapLims = [-10 10];

% set other cond to NaN
cppMeanTopo(repmat(behData.cond~=iC,1,1,eeg.nChans,length(wins)-1)) = NaN;
% split by CoM
cppMeanTopoByCoM = groupMeans(cppMeanTopo,2,repmat(behData.CoM,1,1,eeg.nChans,length(wins)-1)); %[pp CoM chan time]
cppMeanTopoByCoM = diff(cppMeanTopoByCoM,[],2);

% put chans5 electrodes on?
[~, chans5Inds] = ismember(chans5, eeg.chanNames);

subplot(subplotInds);
winInds = isBetween(wins, amplWindows+1); % add 1, as bin is from >min to >=max
topoplot(sq(nanmean(nanmean(cppMeanTopoByCoM(:,:,:,winInds(2:end)),4),1)),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'mapLimits',mapLims,...
    'emarker2',{chans5Inds,'.','k',10,1});
title([labels.cond(iC); labels.baseline; {sprintf('%d:%d ms', amplWindows(1), amplWindows(2))}; {'Effect of CoM'}]);

c = colorbar;

end