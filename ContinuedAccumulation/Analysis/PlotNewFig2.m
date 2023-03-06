function PlotNewFig2(doSuppl)
% PlotNewFig2(doSuppl)
% Plot Figure 2 (behavioural data plus first-order CPP traces +
% topographies)
% if doSuppl==1, plot supplementary figure (just Task 2 CPP traces + topog,
% without excluding CoM trials from confidence)
% 
% Requires: 
% ../../RDMManualConf/Analysis/Saves/CPPAnalysis_CSD.mat
% ../../RDMManualConf/Analysis/Saves/BehDataLoad.mat
% ../../RDMManualConf/Analysis/Saves/FlagArtefacts.mat
% ../../ContinuedAccumulation/Analysis/Saves/CPPAnalysis_CSD_ExclCoMFromCert.mat
% ../../ContinuedAccumulation/Analysis/Saves/CPPAnalysis_CSD_.mat
% ../../ContinuedAccumulation/Analysis/Saves/BehDataLoad.mat
% ../../ContinuedAccumulation/Analysis/Saves/FlagArtefacts3.mat
gf = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',14);

    respWin = [-1000 100];
    showErrBars = 1; % error bars on indiv conditions?

    colours.certainty = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];
%     colours.certainty = [0.3731    0.8864    0.7382; 0.6115    0.5872    0.1107; 0.5505    0.0068    0.4520; .2 .2 .2];
    colours.cond = [0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
    colours.CoM = [0 0 1; 1 0 0; .2 .2 .2]; %[0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2]; % CoM
    colours.conf3 = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];%colours.certainty;
    % colours.conf3 = [0.0937    0.3636    0.6709; 0.6796    0.8444    0.6806; 0.5829    0.2780    0.0740; .2 .2 .2]; % av within confinr1 pairs
    colours.stimDur = tail(crameri('lajolla',5),4);
    colours.acc = [1 0 1; 0 .8 0; .2 .2 .2];
    colours.confInR1 =  [flipud(crameri('roma',6)); .2 .2 .2];

    colours.cmap = crameri('vik');

    lines.widths.certainty = [2 2 2 2];
    lines.styles.certainty = {'-','-','-','--'};
    lines.widths.confInR1 = [2 2 2 2 2 2 2];
    lines.styles.confInR1 = {'-','-','-','-','-','-', '--'};
    lines.alphas = table2struct(array2table([.2 .2 .1 .2],'VariableNames',{'certainty','CoM','confInR1','conf3'}));
   
    if exist('doSuppl','var') && doSuppl % do supplementary figure
        % plot supplementary figure with excludeCoMFromCert=0
        figure();
        PlotTask2Traces([221; 222], respWin, {'Experiment 2: Extinguished', 'Experiment 2: Continued'}, showErrBars, colours, 0, lines); % don't excludeCoMFromCert
    
    %     figure();
        if ~exist('topoplot','file'); eeglab nogui; end
        
        mapLims = [-10 10];
        PlotTask2Topos([2 1 2], {'Task 2'}, colours.cmap, mapLims, 0); % do excludeCoMFromCert

    else % do main figures
        figure();
        PlotTask1Traces([3 2 1:2], respWin, 'Experiment 1', showErrBars, colours, lines)
        PlotBehFigsTask1([3 2 3; 3 2 4; 3 2 5; 3 2 6], {''}, colours, [.5 .9; .5 .9; 1.8 2.8; 350 650]);

%         % with CoM excluded
        figure();
        PlotTask2Traces([321; 322], respWin, {'Experiment 2: Extinguished', 'Experiment 2: Continued'}, showErrBars, colours, 0, lines); % do excludeCoMFromCert
        PlotBehFigsTask2([3 2 3; 3 2 4; 3 2 5; 3 2 6], {''}, colours, [.7 .9; 300 1100; 1.8 2.8; 0 .7]);
    
        figure();
        if ~exist('topoplot','file'); eeglab nogui; end
        mapLims = [-20 20];
        PlotTask1Topo(121, 'Experiment 1', colours.cmap, mapLims);
        PlotTask2Topos([1 2 2], {'Experiment 2'}, colours.cmap, mapLims, 0); % do excludeCoMFromCert

    end

    set(0,'DefaultAxesFontSize',gf); % reset
end

function PlotTask1Traces(subplotInds, respWin, titles, showErrBars, colours, lines)
doStats=1;
    opts.useCSD = 1;
    opts.excludeBadPps = 1; % remove pps with <640 good trials?
    opts.excludeTooFew = 1; % remove pps with <20 per conf3
    opts.excludeByRT = 1; % remove trials outside [100 1500] ms
    opts.doFilt = 1; % whether to plot with loPass filtered data
    
    opts.outFolder = './Saves';
    
    loadFolder = '../../RDMManualConf/Analysis/Saves';
    
    %% load
    
    opts.saveOpts = {'Volt','CSD'};
    opts.saveName = sprintf('CPPAnalysis_%s.mat', opts.saveOpts{1,opts.useCSD+1});
    
    optsNames = fieldnames(opts);
    data = load(fullfile(loadFolder, opts.saveName), optsNames{:}, ...
        'behData', 'cppFilt', 'factors','nF', 'labels',...
        'diffFlips','eeg','amplWindows');
    
    
    %% check things match
    
    dataNames = fieldnames(data);
    data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));
    
    if ~equals(opts, data1)
        warning('loaded data and options do not match');
        keyboard;
    end
    
    %% move into workspace
    struct2workspace(data)

    diffFlips(2) = 1; % High - Low, for difference wave
    labels.diffNames(2) = {'Certain - Maybe'};

    %% plot trace, mean, topo
    
    rInds = isBetween(eeg.respTimes, respWin);
    
    ylims = [-0 42];
    
    for iF = 2
    
        subplot(subplotInds(1,1), subplotInds(1,2), subplotInds(1,3:end));
        set(gca,'ColorOrder',colours.(factors{iF}),'nextplot','replacechildren');
    
        resplocked = groupMeans(cppFilt,3,repmat(permute(behData.(factors{iF}),[1,3,2]),1,size(cppFilt,2))); %[pp t stimDur]
%         resplocked(:,:,end+1) = diff(resplocked(:,:,[1 end]),[],3) * diffFlips(iF); % get diff wave

        h = errorBarPlot(resplocked(:,rInds,:),'area',1,'xaxisvalues',eeg.respTimes(rInds));
    
        if showErrBars==0
            for j=1:size(h,1) % still show diff wave
                h{j,2}.Visible='off'; 
                h{j,1}.LineWidth = 2;
            end % remove error bars
        else % make more transparent
            for j=1:size(h,1)
                h{j,2}.FaceAlpha = lines.alphas.(factors{iF});
                h{j,1}.LineWidth = 2;
            end
        end
        for j=1:size(h,1)
            h{j,1}.LineWidth = lines.widths.(factors{iF})(j);
            h{j,1}.LineStyle = lines.styles.(factors{iF}){j};
        end


        % shade the mean window
        hold on;
        fill(col(repmat(amplWindows,2,1)), [ylims(1,:) fliplr(ylims(1,:))], 'k','FaceAlpha',.2,'LineStyle','none')
    
    
        xticks(-1000:250:500); xtickangle(0);
        xline(0);yline(0);
        xlim(respWin); %ylim(ylims(2,:));
        xlabel('time to intial RT (ms)');
        ylabel('CPP \muV/m^2');
        ylim(ylims(1,:));

        if doStats % regress factor within this cond, at every 100ms bin

            if ~exist('regTab','var')
                regTab = table;
                regTab.pp = nanzscore(col(behData.pp)); % participant for RE
                regTab.stimDur = col(behData.stimDur);
            end
            regTab.(factors{iF}) = nanzscore(col(behData.(factors{iF}))); % add this factor
            times = -1000:100:100;
%             x = movmean(times,2,2,'endpoints','discard');
            formula = sprintf('amplWin ~ 1 + %s*stimDur + (1 | pp)', factors{iF}); % formula to run, in this cond
            
            % run one reg to get number of terms
            regTab.amplWin = nanzscore(col(cppFilt(:,1,:))); % use mean
            fit = fitglme(regTab, formula);
            stats = NaN(length(fit.CoefficientNames)-1,length(times)-1,2);
            ind = strcmp(fit.CoefficientNames,factors{iF});

            for iT = 1:length(times)-1
            
                % take mean within window
                regTab.amplWin = nanzscore(col(nanmean(cppFilt(:,isBetween(eeg.respTimes, times([iT iT+1])),:),2))); % use mean
            
            
                if sum(~isnan(regTab.amplWin)) > 100 % if at least 100 data
                    fit = fitglme(regTab, formula);
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
            pbar(col([stats(find(ind)-1,:,2);stats(find(ind)-1,:,2)])', 'xVals',col([times(1:end-1);times(2:end)]), 'yVal',yVal,'plotargs',{'Color','k','LineWidth',3});
%             h1.LineStyle = 'none'; h1.Marker = '.'; h1.MarkerSize = 20;

        end
    
    
        legend(flip([h{:,1}]), flip(labels.(factors{iF})) , 'Location','Best','Box','off');
    
        title(titles)
    end
end

function PlotTask2Traces(subplotInds, respWin, titles, showErrBars, colours, excludeCoMFromCert, lines)
doStats=1;
    %% set options
    
    opts.useCSD = 1;
    opts.excludeBadPps = 1; % remove pps with <640 good trials?
    opts.excludeTooFew = 1; % remove pps with <20 per conf3
    opts.excludeByRT = 1; % remove trials outside [100 1500] ms
    opts.doFilt = 1; % whether to plot with loPass filtered data
    opts.excludeCoMFromCert = excludeCoMFromCert; % remove CoM trials from behData.conf3
    
    opts.outFolder = './Saves';
        
    %% load
    
    opts.saveOpts = {'Volt','CSD'; '', 'ExclCoMFromCert'};
    opts.saveName = sprintf('CPPAnalysis_%s_%s.mat', opts.saveOpts{1,opts.useCSD+1}, opts.saveOpts{2, opts.excludeCoMFromCert+1});
    
    optsNames = fieldnames(opts);
    data = load(fullfile(opts.outFolder, opts.saveName), optsNames{:}, ...
        'behData', 'cppFilt','cppVarNames','factors','nF', 'labels',...
        'diffFlips', 'eeg','amplWindows');
    
    
    %% check things match
    
    dataNames = fieldnames(data);
    data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));
    
    if ~equals(opts, data1)
        warning('loaded data and options do not match');
        keyboard;
    end
    
    %% move into workspace
    struct2workspace(data)
    
    diffFlips(2) = 1; % High - Low
    labels.diffNames(2) = {'Certain - Maybe'};
    if excludeCoMFromCert
        labels.certainty = {'maybe no-CoM', 'probably no-CoM', 'certain no-CoM'}; % no-com
    else
        labels.certainty = {'maybe CoM/no-CoM', 'probably CoM/no-CoM', 'certain CoM/no-CoM'};
    end
    %% final figure - prep
    
    % for now, strip unused vars out
    
    behDataByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),behData,'UniformOutput',0); %[pp cond tr]
    cppFiltCond = groupMeans(cppFilt,3,repmat(permute(behData.cond,[1,3,2]),1,size(cppFilt,2)),'dim'); %[pp t cond tr]
    
    faceAlpha = .2;
    rInds = isBetween(eeg.respTimes, respWin);
    
    %% plot
    ylims = [0 30;];
    
    for iF = 2 % certainty
        for iC = 1:2
            subplot(subplotInds(iC,:));
            set(gca,'ColorOrder',colours.(factors{iF}),'nextplot','replacechildren');
    
            resplocked = groupMeans(cppFiltCond,4,repmat(permute(behDataByCond.(factors{iF}),[1,4,2,3]),1,size(cppFilt,2))); %[pp t cond fac]
%             resplocked(:,:,:,end+1) = diff(resplocked(:,:,:,[1 end]),[],4) * diffFlips(iF); % get diff

            h = errorBarPlot(sq(resplocked(:,rInds,iC,:)),'area',1,'xaxisvalues',eeg.respTimes(rInds));
            
            if showErrBars==0
                for j=1:size(h,1)
                    h{j,2}.Visible='off'; 
                    h{j,1}.LineWidth = 2;
                end % remove error bars
            else % make more transparent
                for j=1:size(h,1)
                    h{j,2}.FaceAlpha = lines.alphas.(factors{iF});
                    h{j,1}.LineWidth = 2;
                end
            end
            for j=1:size(h,1) % width + style
                h{j,1}.LineWidth = lines.widths.(factors{iF})(j);
                h{j,1}.LineStyle = lines.styles.(factors{iF}){j};
            end

            % shade the mean window
            hold on;
            fill(col(repmat(amplWindows(1,:),2,1)), [ylims(1,:) fliplr(ylims(1,:))], 'k','FaceAlpha',faceAlpha,'LineStyle','none')
    
    
            xticks(-1000:250:500); xtickangle(0);
            xline(0);yline(0);
            xlim(respWin); %ylim(ylims(2,:));
            xlabel('time to initial RT (ms)');
            if iC==1; ylabel('CPP \muV/m^2'); end
            ylim(ylims(1,:));
    
            title(titles(iC))

            if doStats % regress factor within this cond, at every 100ms bin

                if ~exist('regTab','var')
                    regTab = table;
                    regTab.pp = nanzscore(col(behData.pp)); % participant for RE
                    regTab.cond = col(behData.cond);
                end
                regTab.(factors{iF}) = nanzscore(col(behData.(factors{iF}))); % add this factor
                times = -1000:100:100;
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
%                 pbar(stats(1,:,2), 'xVals',movmean(times,2,2,'endpoints','discard'), 'yVal',yVal,'plotargs',{'Color','k','LineWidth',3,'Marker','x'});
    %             h1.LineStyle = 'none'; h1.Marker = '.'; h1.MarkerSize = 20;
    
            end

            if iC==1
                legend(flip([h{:,1}]), flip(labels.(factors{iF})), 'Location','Best','Box','off');
            end
        end
    
    end

end

function PlotBehFigsTask1(subplotInds, titles, colours, ylims)
    % load up data, do some splitting, and then plot the figure
            
    excludeByRT = 1;
    excludeBadPps = 1;
    excludeTooFew = 1;

    loadFolder = '../../RDMManualConf/Analysis/Saves';
    
    load(fullfile(loadFolder, './BehDataLoad.mat'),'behData','labels');
    behData = rmfield(behData, {'dataMatNames', 'dataMat'});
    
    load(fullfile(loadFolder, './FlagArtefacts.mat'), 'isFlagged','ppsToExclude');
    
    %% filter + exclude!!!

    toRemove = isFlagged';
    
    if excludeBadPps
        toRemove(ppsToExclude,:) = 1;     
    end
    
    if excludeByRT
        rtLims = [100 1500];
        toRemove(~isBetween(behData.RT, rtLims)) = 1;
    end
    
    if excludeTooFew
        
        % only remove those trials, not the entire person
        behData.certainty(toRemove) = NaN;
        toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.stimDur);
    
        toRemove(toRemove2) = 1;
    end
    
    % exclude from each field
    behNames = fieldnames(behData);
    nBeh = length(behNames);
    
    % remove flagged trials and bad pps
    for i = 1:nBeh
        
        behData.(behNames{i})(repmat(toRemove,1,1,size(behData.(behNames{i}),3))) = NaN;
        
    end
        
    
    %% paper figures
    
       
    behDataByDur= structfun(@(x) groupMeans(x,2,behData.stimDur,'dim'), behData, 'UniformOutput',0);
    
    ylabels = { 'p(correct)', 'p(correct)', 'mean confidence', 'mean RT (ms)'};
    varNames2 = {'acc','accByCert','certainty','rtByCert'};
    % 4-6: certainty: acc*cond, and Com + confInr1
    for i = 1:size(subplotInds,1)
        subplot(subplotInds(i,1), subplotInds(i,2), subplotInds(i,3));
        set(gca,'ColorOrder',colours.stimDur,'nextplot','replacechildren');

        if strcmp(varNames2{i}, 'acc')
    
            h = errorBarPlot(nanmean(behDataByDur.(varNames2{i}),3), 'type','bar');
            h.FaceColor = 'flat';
            h.CData = colours.stimDur(1:3,:);
            set(gca,'XTickLabel',labels.stimDur);

            ylabel(ylabels{i});
            xlabel('stimulus duration');
            ylim(ylims(i,:));
            xlim([ .25 3.75]);

%             % make legend
%             hold on; emptyLegend(3, { {'FaceColor', colours.stimDur(1,:)},{'FaceColor', colours.stimDur(2,:)},{'FaceColor', colours.stimDur(3,:)} }, {}, labels.stimDur, {'Location','Best'}, @bar);
%             title(titles{1});
        
        elseif strcmp(varNames2{i},'accByCert')
            accByCert = groupMeans(behData.acc,2,behData.certainty); %[pp dur cert]
            h = errorBarPlot(accByCert, 'type','bar');
            h.FaceColor = 'flat';
            h.CData = colours.certainty(1:3,:);
            set(gca,'XTickLabel',labels.certainty);

            ylim(ylims(i,:));
            xlim([ .25 3.75]);

            ylabel(ylabels(i));
            xlabel('confidence rating');

        elseif strcmp(varNames2{i},'rtByCert')
            % RT by Cert/acc/dur
            xFactor = 'stimDur';
            legendFactor = 'certainty';
            if strcmp(legendFactor,'stimDur'), inds = [1,4,2,3,5];
            else, inds = [1,2,4,3,5];
            end
            durDataByAcc.RT = groupMeans(behDataByDur.RT,3,behDataByDur.acc,'dim'); %[ pp dur acc tr]
            durDataByAcc.certainty = groupMeans(behDataByDur.certainty,3,behDataByDur.acc,'dim'); %[ pp dur acc tr]
            rtByDurAccCert = groupMeans(durDataByAcc.RT,4,durDataByAcc.certainty,'dim');
            % colour is dur, xlabel is cert
            set(gca,'ColorOrder',colours.(legendFactor)(1:3,:),'nextplot','replacechildren');
            hold on;
            h = errorBarPlot(permute(nanmean(rtByDurAccCert(:,:,1,:,:),5),inds),'type','line','plotargs',{'LineStyle','--','LineWidth',2}); % errors
            h1 = errorBarPlot(permute(nanmean(rtByDurAccCert(:,:,2,:,:),5),inds),'type','line','plotargs',{'LineStyle','-','LineWidth',2});
            set(gca,'XTick',1:3,'XTickLabel',labels.(xFactor));
            ylabel(ylabels(4));
            xlabel('stimulus duration');
            xlim([0.5 3.5]);
            emptyLegend(5, { {'Color', colours.(legendFactor)(1,:)},{'Color', colours.(legendFactor)(2,:)},{'Color', colours.(legendFactor)(3,:)}, {'k','LineStyle','-'}, {'k','LineStyle','--'} }, {'LineWidth',2}, [labels.(legendFactor) {'Correct'},{'Error'}], {'Location','Best'});
        
        else

            durDataByAcc.(varNames2{i}) = groupMeans(behDataByDur.(varNames2{i}),3,behDataByDur.acc,'dim'); %[ pp cond acc tr]
            h = errorBarPlot(permute(nanmean(durDataByAcc.(varNames2{i}),4),[1,3,2]), 'type','bar');
            
            set(gca,'XTickLabel',labels.acc);
            ylim(ylims(i,:));

            ylabel(ylabels(i));
            xlabel('accuracy');

            legend(labels.stimDur, 'Location','Best');

        end
        if i==1
            t = title(titles); 
            t.Position = t.Position + [1.3 0 0];
        end
    end


end

function PlotBehFigsTask2(subplotInds, titles, colours, ylims)
    % load up data, do some splitting, and then plot the figure
            
    excludeByRT = 1;
    excludeBadPps = 1;
    excludeTooFew = 1;
    
    load('./Saves/BehDataLoad.mat','behData','labels');
    
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
        
    labels.bothAcc = {'initial','final'};
    labels.bothRT = labels.bothAcc;
    
    %%   
       
    behDataByCond = structfun(@(x) groupMeans(x,2,behData.cond,'dim'), behData, 'UniformOutput',0);
    
    % combine init and conf acc
    behDataByCond.bothAcc = nancat(3, groupMeans(behData.acc,2,behData.cond), groupMeans(behData.confAcc,2,behData.cond)); %[pp cond init/conf]
    
    behDataByCond.bothRT = nancat(3, groupMeans(behData.RT,2,behData.cond), groupMeans(behData.confRT,2,behData.cond)); %[pp cond init/conf]

    % 1 = acc: time*cond
    varNames1 = {'bothAcc','bothRT'};
    ylabels1 = {'p(correct)', 'mean initial RT (ms)'};
    for i = 1
        subplot(subplotInds(i,1), subplotInds(i,2), subplotInds(i,3));
        
        set(gca,'ColorOrder',colours.cond,'nextplot','replacechildren');
        
        h = errorBarPlot(permute(behDataByCond.(varNames1{i}),[1,3,2]), 'type','bar');%','plotargs',{'LineWidth',2}); % permute to make cond the legend
        set(gca,'XTick',1:2,'XTickLabel',labels.(varNames1{i}));
        ylabel(ylabels1{i});
        xlabel('response');
        ylim(ylims(i,:));

        if i==1; legend(h, labels.cond,'Location','Best'); end
        
        if i==1
            t = title(titles); 
            t.Position = t.Position + [1.3 0 0];
        end
    end

    % resp time * acc * cond. use confAcc for 2nd order RT
%     rt = cat(3, groupMeans(behDataByCond.RT, 3, behDataByCond.acc), groupMeans(behDataByCond.confRT, 3, behDataByCond.confAcc));
    subplot(subplotInds(2,1), subplotInds(2,2), subplotInds(2,3));

%     set(gca,'ColorOrder',colours.cond,'nextplot','replacechildren');
% 
%     h = errorBarPlot(permute(rt,[1,3,2]), 'type','bar');%','plotargs',{'LineWidth',2}); % permute to make cond the legend
%     set(gca,'XTick',1:4,'XTickLabel',labels.acc);
%     ylabel(ylabels1{2});
%     xlabel({[labels.bothRT{1}, '            ' labels.bothRT{2}]; 'response'});
%     ylim([300 1000]);
%     ylim(ylims(2,:));
% 
    
    % initial RT first
    % RT into [pp cond acc cert tr]
    rtByCondAcc = groupMeans(behDataByCond.RT,3,behDataByCond.acc,'dim'); %[pp cond acc tr]
    certByCondAcc = groupMeans(behDataByCond.certainty,3,behDataByCond.acc,'dim');
    rtByCondAccCert = groupMeans(rtByCondAcc,4,certByCondAcc,'dim'); %[pp cond acc cert tr]

    % also get confRT
    rtConfByCondAcc = groupMeans(behDataByCond.confRT,3,behDataByCond.acc,'dim'); %[pp cond acc tr]
    rtConfByCondAccCert = groupMeans(rtConfByCondAcc,4,certByCondAcc,'dim'); %[pp cond acc cert tr]
    

%     lineStyles = {'--','-'};
    set(gca,'ColorOrder',colours.cond(1:2,:),'nextplot','replacechildren', 'LineStyleOrder',{'--','-'});
    hold on;
    h = errorBarPlot(reshape(permute(nanmean(rtByCondAccCert,5),[1,4,2,3]),30,[],4), 'type','line','plotargs',{'LineWidth',2});
    
    set(gca,'XTick',1:6,'XTickLabel',labels.certainty);
    ylabel(ylabels1(2));
%     ylim(ylims(2,:));
    yt = yticks;
    yticks(yt(mod(yt,100)==0)); % only 100-markers
    xlabel('confidence report');
    xlim([.5 3.5]);
    hold on; emptyLegend(4, { {'Color', colours.cond(2,:),'LineStyle','-'},{'Color', colours.cond(1,:),'LineStyle','-'}, {'k','LineStyle','-'}, {'k','LineStyle','--'} }, {'LineWidth',2}, [flip(labels.cond) {'Initial Correct'},{'Initial Error'}], {'Location','Best'});

%     % plot conf RT on right yaxis
%     yyaxis right
%     set(gca,'ColorOrder',colours.cond(1:2,:),'nextplot','replacechildren', 'LineStyleOrder',{'--','-'});
%     hold on;
%     % now plot conf RT
%     h = errorBarPlot(reshape(permute(nanmean(rtConfByCondAccCert,5),[1,4,2,3]),30,[],4), 'type','line','plotargs',{'LineWidth',2},'xaxisvalues',[7:12]');
% 
%     ax = gca;
%     ax.YAxis(2).Color = 'k';
%     ylabel('confidence RT (ms)');
%     yticks(100:200:1200);
% 
%     xlabel([labels.bothRT{1} ' response', '                       ', labels.bothRT{2} ' response']);
%     xlim([.5 6.5]);
%     hold on; emptyLegend(4, { {'Color', colours.cond(2,:),'LineStyle','-'},{'Color', colours.cond(1,:),'LineStyle','-'}, {'k','LineStyle','-'}, {'k','LineStyle','--'} }, {'LineWidth',2}, [flip(labels.cond) {'Initial Correct'},{'Initial Error'}], {'Location','Best'});


    ylabels = { 'mean confidence-in-final-choice', 'p(change-of-mind)'};
    varNames2 = {'certainty', 'CoM'};
    % 4-6: conf3: acc*cond, and Com + confInr1
    for i = 1:length(varNames2)
        subplot(subplotInds(2+i,1), subplotInds(2+i,2), subplotInds(2+i,3));
        set(gca,'ColorOrder',colours.cond,'nextplot','replacechildren');
    
        condDataByAcc.(varNames2{i}) = groupMeans(behDataByCond.(varNames2{i}),3,behDataByCond.acc,'dim'); %[ pp cond acc tr]
        h = errorBarPlot(permute(nanmean(condDataByAcc.(varNames2{i}),4),[1,3,2]), 'type','bar');
        set(gca,'XTickLabel',labels.acc);

%         h = violin(reshape(permute(nanmean(condDataByAcc.(varNames2{i}),4),[1,2,3]),[],4),'x',[.7 1.3 2.7 3.3], 'medc',[],'facecolor',repmat(colours.cond(1:2,:),2,1),'facealpha',1,'plotlegend',0);
%         set(gca,'XTick',[1 3], 'XTickLabel',labels.acc);

        ylim(ylims(2+i,:));
        ylabel(ylabels(i));
        xlabel('initial accuracy');

    end

end


function PlotTask1Topo(subplotInds, titles, cmap, mapLims)

    opts.useCSD = 1;
    opts.excludeBadPps = 1; % remove pps with <640 good trials?
    opts.excludeTooFew = 1; % remove pps with <20 per conf3
    opts.excludeByRT = 1; % remove trials outside [100 1500] ms
    opts.doFilt = 1; % whether to plot with loPass filtered data
    
    opts.outFolder = './Saves';
    
    loadFolder = '../../RDMManualConf/Analysis/Saves';
    
    %% load
    
    opts.saveOpts = {'Volt','CSD'};
    opts.saveName = sprintf('CPPAnalysis_%s.mat', opts.saveOpts{1,opts.useCSD+1});
    
    optsNames = fieldnames(opts);
    data = load(fullfile(loadFolder, opts.saveName), optsNames{:}, ...
        'behData', 'factors','nF', 'labels',...
        'diffFlips', 'eeg','cppTopoAmplWindow','chans5','cmap');

    
    
    %% check things match
    
    dataNames = fieldnames(data);
    data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));
    
    if ~equals(opts, data1)
        warning('loaded data and options do not match');
        keyboard;
    end
    
    %% move into workspace
    struct2workspace(data)


    %% plot
    [~, chans5Inds] = ismember(chans5, eeg.chanNames);
    
    for iF = 2

        subplot(subplotInds);
        topoplot(sq(nanmean(nanmean(cppTopoAmplWindow(:,:,:,1),2),1)),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', mapLims,...
            'emarker2',{chans5Inds,'.','k',20,1});

        colorbar;
        title(titles);
    end

end


function PlotTask2Topos(subplotInds, titles, cmap, mapLims, excludeCoMFromCert)

    %% set options
    
    opts.useCSD = 1;
    opts.excludeBadPps = 1; % remove pps with <640 good trials?
    opts.excludeTooFew = 1; % remove pps with <20 per conf3
    opts.excludeByRT = 1; % remove trials outside [100 1500] ms
    opts.doFilt = 1; % whether to plot with loPass filtered data
    opts.excludeCoMFromCert = excludeCoMFromCert; % remove CoM trials from behData.conf3
    
    opts.outFolder = './Saves';
    

    %% load
    
    opts.saveOpts = {'Volt','CSD'; '', 'ExclCoMFromCert'};
    opts.saveName = sprintf('CPPAnalysis_%s_%s.mat', opts.saveOpts{1,opts.useCSD+1}, opts.saveOpts{2, opts.excludeCoMFromCert+1});
    
    optsNames = fieldnames(opts);
    data = load(fullfile(opts.outFolder, opts.saveName), optsNames{:}, ...
        'behData','factors','nF', 'labels',...
        'diffFlips', 'eeg','cppTopoAmplWindow','chans5','cmap');
    
    
    %% check things match
    
    dataNames = fieldnames(data);
    data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));
    
    if ~equals(opts, data1)
        warning('loaded data and options do not match');
        keyboard;
    end
    
    %% move into workspace
    struct2workspace(data)

    %%  plot
    [~, chans5Inds] = ismember(chans5, eeg.chanNames);
    
    for iF = 2 % cert

        for iC=1 % is pre, so combine
            
            subplot(subplotInds(iC,1), subplotInds(iC,2), subplotInds(iC,3));
            topoplot(sq(nanmean(nanmean(cppTopoAmplWindow(:,:,:,1),2),1)),...
                eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', mapLims,...
                'emarker2',{chans5Inds,'.','k',20,1});

            colorbar;
            title(titles(iC));
        end        
    end

end

