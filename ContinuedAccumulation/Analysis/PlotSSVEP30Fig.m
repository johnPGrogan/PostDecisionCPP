function PlotSSVEP30Fig()
% Plot 30Hz SSVEP 
% Requires:
% ContinuedAccumulation/Analysis/Saves/SSVEPAnalysis_CSD_.mat
% ContinuedAccumulation/Analysis/Saves/GetSSVEPTopo_CSD.mat
% ContinuedAccumulation/Analysis/Saves/ExtractEpochs.mat
% ContinuedAccumulation/Analysis/Saves/BehDataLoad.mat
% ContinuedAccumulation/Analysis/Saves/FlagArtefacts3.mat


respWin = [-100 1000; -500 1000];
showErrBars = 1; % show error bars in phase figures?
faceAlpha = .1;

colours.certainty = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];
% colours.certainty = [0.3731    0.8864    0.7382; 0.6115    0.5872    0.1107; 0.5505    0.0068    0.4520; .2 .2 .2];
colours.cond = [0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
colours.CoM = [0 0 1; 1 0 0; .2 .2 .2]; %[0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2]; % CoM
colours.conf3 = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];%colours.certainty;
% colours.conf3 = [0.0937    0.3636    0.6709; 0.6796    0.8444    0.6806; 0.5829    0.2780    0.0740; .2 .2 .2]; % av within confinr1 pairs
colours.confInR1 =  [flipud(crameri('roma',6)); .2 .2 .2];
colours.cmap = crameri('vik');
% colours.cmap = colours.cmap(129:end,:); % if using baselined

widths.certainty = [2 2 2 2];
styles.certainty = {'-','-','-','--'};
widths.conf3 = widths.certainty;
styles.conf3 = styles.certainty;
widths.CoM = [2 2 2];
styles.CoM = {'-','-','--'};
widths.acc = [1.5 3 2];
styles.acc = {'-','-','--'};
widths.confInR1 = [2 2 2 2 2 2];
styles.confInR1 = {'-','-','-','-','-','-'};
lines.widths = widths; 
lines.styles = styles;
lines.alphas = table2struct(array2table([.2 .2 .1 .2],'VariableNames',{'certainty','CoM','confInR1','conf3'}));


% % plot phase err + corr, by conf + com
figure();
PlotSSVEPFigs(respWin(1,:), showErrBars, faceAlpha, colours, lines, [4 5 2]);

%% plot topographies - average over all continued good trials
% plots the beta coeff from regression of confInR1 in each acc separately

figure();
PlotSSVEPTopo([1 1 1], [700:100:1000], colours, [1]); % just mean across all conds/trials

% figure();
% topoWin = [700:100:1000]; %ms,  end of 100ms times to plot mean of , i.e.[400 500] = 300:500ms
% mapLim = [.06]; % from zero to this
% subplotInds = [321 322; 323 324; 325 326;]; % just pass 1st row to only do grand mean, others do CoM effect in err/corr
% PlotSSVEPTopoBeta(subplotInds, topoWin, colours, mapLim, [4 5 2]);% plot per factors




end


function PlotSSVEPFigs(respWin, showErrBars, faceAlpha, colours, lines, iFs)

    opts.useCSD = 1;
    opts.excludeByRT = 1; % remove trials outside [100 1500] ms
    opts.excludeBadPps = 1;
    opts.excludeTooFew = 1;
    opts.useCPPSSVEP = 0;
    undoBaseline = 0;
    
    opts.outFolder = './Saves';
    
    
    %% load
    
    opts.saveOpts = {'Volt','CSD'; '', 'CPPChans'};
    opts.saveName = sprintf('SSVEPAnalysis_%s_%s.mat', opts.saveOpts{1,opts.useCSD+1}, opts.saveOpts{2, opts.useCPPSSVEP+1});
    
    optsNames = fieldnames(opts);
    data = load(fullfile(opts.outFolder, opts.saveName), optsNames{:}, ...
        'ssvep','lockNames','windows','wins','nTimes',...
        'behData','factors','nF', 'labels','baseline','doBaseline');
    
    
    %% check things match
    
    dataNames = fieldnames(data);
    data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));
    
    if ~equals(opts, data1)
        warning('loaded data and options do not match');
        keyboard;
    end

    %% move into workspace
    struct2workspace(data);

    %% undo baseline

    if undoBaseline==1 && doBaseline==1
        for i = 1:2
            ssvep.(lockNames{i}) = exp(ssvep.(lockNames{i})) .* baseline;
        end
    end

    %%

    times = data.windows;

    wins.resp = [700 1000];
    labels.confInR1 = {'certain CoM', 'probably CoM', 'maybe CoM', 'maybe no-CoM', 'probably no-CoM', 'certain no-CoM'};
    labels.conf3 = {'certain/probably CoM', 'maybe CoM/no-CoM', 'probably/certain no-CoM'};
    labels.certainty = {'maybe CoM/no-CoM', 'probably CoM/no-CoM', 'certain CoM/no-CoM'};
    names.confInR1 = 'confidence-in-initial-choice';
    names.certainty = 'final-certainty';
    names.conf3 = 'confidence-in-initial-choice: binned';
    names.CoM = 'change-of-mind';

    % invert CoM so that change comes first
    behData.CoM = 2 - behData.CoM;
    labels.CoM = flip(labels.CoM);

    %% split things
    
    % behData ones
    behDataByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), behData, 'UniformOutput',0);
    
    % phase data ones
    ssvepByAcc = structfun(@(x) groupMeans(x, 3, repmat(permute(behData.acc,[1,3,2]),1,size(x,2)), 'dim'), ssvep, 'UniformOutput',0); %[pp t acc tr]
    
    
    %% final figure - acc*fac
    
    % just post-response? trying to answer the question, do people stop paying
    % attention when they are highly confident?
    
    ylims = [-1 .2;];
    iPh = find(strcmp(lockNames, 'resp')); % just post-resp
    nF=2;

    %% just fac

%     for iF = 3 % +1
%         phaseByFac = groupMeans(phase.(lockNames{iPh}), 3, repmat(permute(behData.(factors{iF+1}),[1,3,2]),1,nTimes.(lockNames{iPh}))); %[pp t fac]   
%         
%         subplot(2,3,1);%(iF-1)*3+1);
% 
%         set(gca,'ColorOrder',colours.(factors{iF+1}),'nextplot','replacechildren');
%         h = errorBarPlot(phaseByFac, 'area',1, 'xaxisvalues', times.(lockNames{iPh}));
%         if ~showErrBars
%             for i=1:size(h,1), h{i,2}.Visible = 'off'; end
%         end
%         for j=1:size(h,1)
%             h{j,1}.LineWidth = lines.widths.(factors{iF+1})(j);
%             h{j,1}.LineStyle = lines.styles.(factors{iF+1}){j};
%         end
%         hold on;
%         fill(col(repmat(times.resp,2,1)), [ylims(1,:) fliplr(ylims(1,:))], 'k','FaceAlpha',faceAlpha,'LineStyle','none')
%    
%         xlabel('time from initial RT (ms)');
%         xline(0);
%         yline(0);
%         ylim(ylims(1,:));
%         xlim(respWin);
%         
%         if iF==2; legend([h{:,1}], labels.(factors{iF+1}), 'Location','Best'); end
%         if iF; ylabel('30Hz SSVEP'); end
%         title('Continued Evidence');
%     
%     end

    %% acc*fac
    nF = length(iFs);
    for i = 1:nF % +1
        iF=iFs(i);
        data = groupMeans(ssvepByAcc.(lockNames{iPh}), 4, repmat(permute(behDataByAcc.(factors{iF}),[1,4,2,3]),1,nTimes.(lockNames{iPh}))); %[pp t ac fac]
    
        if any(strcmpi(factors{iF},{'certainty','conf3','CoM'}))
            % average over confInR1 first
            data = groupMeans(ssvepByAcc.(lockNames{iPh}), 4, repmat(permute(behDataByAcc.confInR1,[1,4,2,3]),1,nTimes.(lockNames{iPh}))); %[pp t ac confinr1]
            facByConfInR1 = groupMeans(behDataByAcc.(factors{iF}),3,behDataByAcc.confInR1); %[pp acc confInR1]
            data = groupMeans(data,4,repmat(permute(facByConfInR1,[1,4,2,3]),1,nTimes.(lockNames{iPh}))); %[pp t cond fac]
        end
        for iA = 1:2
            subplot(nF,2,(i-1)*2+iA);
    
            set(gca,'ColorOrder',colours.(factors{iF}),'nextplot','replacechildren');
            h = errorBarPlot(sq(data(:,:,iA,:)), 'area',1, 'xaxisvalues', times.(lockNames{iPh}));
            if ~showErrBars
                for i=1:size(h,1), h{i,2}.Visible = 'off'; end
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
            hold on;
            fill(col(repmat(wins.resp,2,1)), [ylims(1,:) fliplr(ylims(1,:))], 'k','FaceAlpha',faceAlpha,'LineStyle','none')
       
            xlabel('time from initial RT (ms)');
            xline(0);
            yline(0);
            ylim(ylims(1,:));
            xlim(respWin);

            % do stats
            if ~exist('regTab','var')
                regTab = table;
                regTab.pp = nanzscore(col(behData.pp)); % participant for RE
                regTab.acc = col(behData.acc);
            end
            regTab.(factors{iF}) = nanzscore(col(behData.(factors{iF}))); % add this factor
            t = -100:100:1000;
%             x = movmean(times,2,2,'endpoints','discard');
            formula = sprintf('amplWin ~ 1 + %s + (1 | pp)', factors{iF}); % formula to run, in this cond
            
            % run one reg to get number of terms
            regTab.amplWin = nanzscore(col(ssvep.resp(:,1,:))); % use mean
            fit = fitglme(regTab(regTab.acc==(iA-1),:), formula);
            stats = NaN(length(fit.CoefficientNames)-1,length(times)-1,2);

            for iT = 1:length(t)-1
            
                % take mean within window
                regTab.amplWin = nanzscore(col(nanmean(ssvep.resp(:,isBetween(times.resp, t([iT iT+1])),:),2))); % use mean
            
            
                if sum(~isnan(regTab.amplWin)) > 100 % if at least 100 data
                    fit = fitglme(regTab(regTab.acc==(iA-1),:), formula);
                    stats(:,iT,1) = fit.Coefficients.Estimate(2:end); % beta
                    stats(:,iT,2) = fit.Coefficients.pValue(2:end); % p-value
                end
            end

            % plot sig timepoints
%             plot(, stats(:,:,2) < alpha, '.k', 'MarkerSize',20)
            hold on;
            yVal = min(ylim); 
            % duplicate to get it from start:end of a single window if
            % needed
            pbar(col([stats(1,:,2);stats(1,:,2)])', 'xVals',col([t(1:end-1);t(2:end)]), 'yVal',yVal,'plotargs',{'Color','k','LineWidth',3});
%             h1.LineStyle = 'none'; h1.Marker = '.'; h1.MarkerSize = 20;




            if iA==1
                legend([h{:,1}], labels.(factors{iF}), 'Location','Best')
%                 if iF; ylabel('30Hz SSVEP'); end
            end
            if i==1; title([ {['Initial ' labels.acc{iA}]}; names.(factors{iF})]); 
            else; title(names.(factors{iF})); end
            if iA==1; ylabel('30Hz SSVEP'); end

        end
    
    end

%     %% acc*fac means
%     meanPhase.(lockNames{iPh}) = sq(nanmean(phase.(lockNames{iPh})(:,isBetween(times.(lockNames{iPh}), times.(lockNames{iPh})),:),2));
% 
%     meanssvepByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), meanPhase, 'UniformOutput',0);
%     for i = 1:length(iFs)
%         iF=iFs(i);
%         data = groupMeans(meanssvepByAcc.(lockNames{iPh}), 3, behDataByAcc.(factors{iF})); %[pp ac fac]
%     
% %         for iA = 1:2
%             subplot(nF,3,(i)*3);
%     
%             set(gca,'ColorOrder',colours.(factors{iF}),'nextplot','replacechildren');
%             h = errorBarPlot(data,'type','bar');
% %             h = violin(reshape(permute(data,[1 3 2]),[],2*size(data,3)),'medc',[],'facecolor',repmat(colours.(factors{iF})(1:end-1,:),2,1),'facealpha',1,'plotlegend',0);
%             % violin looks awful
%             % super dots instead
%             
% %             % don't for now - maybe suppl fig
% %             hold on
% %             for ii = 1:size(data,3)
% %                 plot(h(ii).XEndPoints+(rand(30,2)-.5)/20,sq(data(:,:,ii)),'.','Color','k');%h(ii).FaceColor);
% %             end
% 
% 
% %             ylim(minMax([h.YData],'all') .* [.8 1.2]');
%             xticks(1:2); xticklabels(labels.acc)
%             xlim([.5 2.5]);
%             ylabel('mean 30Hz SSVEP');
%             if i==1; title('mean within grey window'); end
% %             legend(labels.(factors{iF}),'Location','Best');
%             
% %         end
%     
%     end

%     SuperTitle('Continued Post-Decision Evidence');
end




function PlotSSVEPTopo(subplotInds, topoWins, colours, mapLim)

    useCSD = 1;
    excludeBadPps = 1; % these are already excluded, as are isFlagged, during GetPhaseTagging.m
    excludeTooFew = 1; % remove if <10 trials per cert*cond
    excludeByRT = 1; % bad pps and flagged trials already gone
    % bad pps are already excluded
    undoBaseline = 1;
    
    outFolder = './Saves';
    load('./Saves/ExtractEpochs.mat');
    load('./Saves/BehDataLoad.mat','behData', 'labels','behTab');
    load('./Saves/FlagArtefacts3.mat','isFlagged','ppsToExclude');
    
    if useCSD
        topoName = fullfile(outFolder, 'GetSSVEPTopoCSD.mat');
    else
        topoName = fullfile(outFolder, 'GetSSVEPTopo.mat');
    end
    
    load(topoName, 'ssvepTopoResp', 'respWindows'); %[pp chan t tr]
    
    if undoBaseline==1
        d = load(topoName,'doBaseline','baseline');
        if d.doBaseline==1 % un-baseline it
            ssvepTopoResp = exp(ssvepTopoResp) .* d.baseline;
        end
    end
           
    %% other exclusions?
    toRemove = isFlagged'; % flagged trials or cpp outliers
    
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
        toRemove2 = GetCellsWithLowN(10, x);
    
        toRemove(toRemove2) = 1;
    
    end
    
    ssvepTopoResp(repmat(permute(toRemove,[1,3,4,2]), 1, eeg.nChans, size(ssvepTopoResp,3))) = NaN; % remove
    %% plot
    
    mapLims = [-1 1] .* mapLim;
    iT = ismember(respWindows(2:end), topoWins); % all times to average over
    meanTopo = nanmean(ssvepTopoResp(:,:,iT,:),3); % av over times


    visChanNames = {'A23'}; % midline only

    [~, visChanInds] = ismember(visChanNames, eeg.chanNames);
    
    subplot(subplotInds(1,1),subplotInds(1,2),subplotInds(1,3))
    topoplot(sq(nanmean(nanmean(meanTopo,4),1)),...
        eeg.chanlocs, 'colormap',colours.cmap, ...
         'electrodes','off', 'emarker2', {visChanInds, '.','k',20,1});
    c = colorbar;%('Limits', mapLims); % there are no negatives
    title('grand mean')

%     if size(subplotInds,1) > 1
%         % split by init acc + CoM, topoplot CoM effects
% 
%         for i = 1:(size(subplotInds,1)-1)
%             
%             topoByCoMAcc = sq(groupMeans(meanTopo,4,repmat(permute(behData.CoM + behData.acc*10, [1 3 4 2]),1,1,128,1)));
%             %[pp ch conds(4)]
%             comEffect = topoByCoMAcc(:,:,[2 4]) - topoByCoMAcc(:,:,[1 3]); % CoM - no change
% 
%             subplot(subplotInds(i+1,1),subplotInds(i+1,2),subplotInds(i+1,3))
%             topoplot(sq(nanmean(comEffect(:,:,i))),...
%                 eeg.chanlocs, 'colormap',colours.cmap, 'mapLimits', mapLims, ...
%                  'electrodes','off', 'emarker2', {visChanInds, '.','k',20,1});
%             c = colorbar;%('Limits', [-mapLim mapLim]); % there are no negatives
%             title(['CoM effect: ', labels.acc{i}]);
%         end
%     end


end


function PlotSSVEPTopoBeta(subplotInds, topoWins, colours, mapLim, iFs)
% do regressions at each channel, plot beta

    factors = {'acc','certainty','CoM','confInR1','conf3'};
    

    useCSD = 1;
    excludeBadPps = 1; % these are already excluded, as are isFlagged, during GetPhaseTagging.m
    excludeTooFew = 1; % remove if <10 trials per cert*cond
    excludeByRT = 1; % bad pps and flagged trials already gone
    % bad pps are already excluded
    undoBaseline = 0;
    
    outFolder = './Saves';
    load('./Saves/ExtractEpochs.mat');
    load('./Saves/BehDataLoad.mat','behData', 'labels','behTab');
    load('./Saves/FlagArtefacts3.mat','isFlagged','ppsToExclude');
    
    if useCSD
        topoName = fullfile(outFolder, 'GetSSVEPTopoCSD.mat');
    else
        topoName = fullfile(outFolder, 'GetSSVEPTopo.mat');
    end
    
    load(topoName, 'ssvepTopoResp', 'respWindows');
    
    if undoBaseline==1
        d = load(topoName,'doBaseline','baseline');
        if d.doBaseline==1 % un-baseline it
            ssvepTopoResp = exp(ssvepTopoResp) .* d.baseline;
        end
    end
               
    %% other exclusions?
    toRemove = isFlagged'; % flagged trials or cpp outliers
    
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
        toRemove2 = GetCellsWithLowN(10, x);
    
        toRemove(toRemove2) = 1;
    
    end
    
    ssvepTopoResp(repmat(permute(toRemove,[1,3,4,2]), 1, eeg.nChans, size(ssvepTopoResp,3))) = NaN; % remove
    %% plot
    
    mapLims = [-1 1] .* mapLim;
    iT = ismember(respWindows(2:end), topoWins); % all times to average over
    meanTopo = nanmean(ssvepTopoResp(:,:,iT,:),3); % av over times

    %% do regressions on each channel

    if ~exist('regTab','var')
        regTab = table;
        regTab.pp = nanzscore(col(behData.pp)); % participant for RE
        regTab.cond = col(behData.cond);
        regTab.acc = col(behData.acc);
    end
    iF=4; % placeholder
    regTab.(factors{iF}) = nanzscore(col(behData.(factors{iF}))); % add this factor
    formula = sprintf('amplWin ~ 1 + %s + (1 | pp)', factors{iF}); % formula to run, in this cond

    % run one reg to get number of terms
    regTab.amplWin = nanzscore(col(meanTopo(:,1,1,:))); % use first chan
    fit = fitglme(regTab(regTab.cond==2 & regTab.acc==1,:), formula);


    visChanNames = {'A23'}; % midline only

    [~, visChanInds] = ismember(visChanNames, eeg.chanNames);
    nF = length(iFs);
    for i = 1:nF
        iF = iFs(i);
        regTab.(factors{iF}) = nanzscore(col(behData.(factors{iF}))); % add this factor
        formula = sprintf('amplWin ~ 1 + %s + (1 | pp)', factors{iF}); % update with factor name

        stats = NaN(length(fit.CoefficientNames)-1,eeg.nChans,2,2);
    
        for iCh = 1:eeg.nChans
    
            % take mean within window
            regTab.amplWin = nanzscore(col(meanTopo(:,iCh,1,:))); % each channel, within post-dec window
    
    
            if sum(~isnan(regTab.amplWin)) > 100 % if at least 100 data
                for iA = 1:2
                    fit = fitglme(regTab(regTab.cond==2 & regTab.acc==iA-1,:), formula);
                    stats(:,iCh,iA,1) = fit.Coefficients.Estimate(2:end); % beta
                    stats(:,iCh,iA,2) = fit.Coefficients.pValue(2:end); % p-value
                end
            end
        end
    
        
        for iA=1:2 % acc

        % highlight sig chans?
%         visChanInds = find(stats(1,:,iA,2)<.0001);
        
            subplot(subplotInds(i, iA))
            topoplot(stats(:,:,iA,1),...
                eeg.chanlocs, 'colormap',colours.cmap, 'mapLimits',mapLims,...
                 'electrodes','off', 'emarker2', {visChanInds, '.','k',20,1});
            c = colorbar;%('Limits', [-mapLim mapLim]); % there are no negatives
            title([labels.acc{iA} ': \beta coefficient for ', factors{iF}]);
            
        end
    end

end