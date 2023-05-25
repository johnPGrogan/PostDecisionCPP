function PlotNewFigS6()
% PlotNewFigS6
% Plot Figure S6 - CPP traces + topoplots on voltage (non-CSD) data
% expt 1 pre-response, expt 2 pre-response + post-choice CPPs
% 
% Requires: 
% ../../RDMManualConf/Analysis/Saves/CPPAnalysis_Volt.mat
% ../../ContinuedAccumulation/Analysis/Saves/CPPAnalysis_Volt_ExclCoMFromCert.mat
% ../../ContinuedAccumulation/Analysis/Saves/CPPAnalysis_Volt_.mat
gf = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',14);

    respWin = [-1000 100];
    showErrBars = 1; % error bars on indiv conditions?

    colours.certainty = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];
    colours.confInR1 =  [flipud(crameri('roma',6)); .2 .2 .2];

    colours.cmap = crameri('vik');

    lines.widths.certainty = [2 2 2 2];
    lines.styles.certainty = {'-','-','-','--'};
    lines.widths.confInR1 = [2 2 2 2 2 2 2];
    lines.styles.confInR1 = {'-','-','-','-','-','-', '--'};
    lines.alphas = table2struct(array2table([.2 .2 .1 .2],'VariableNames',{'certainty','CoM','confInR1','conf3'}));
   
   
    figure();
    % expt 1 pre-choice CPP by confidence
    PlotTask1Traces([3 2 1:2], respWin, 'Experiment 1', showErrBars, colours, lines)


    % expt 2 pre-choice CPP by confidence (no-CoM only)
    PlotTask2Traces([323; 324], respWin, {'Experiment 2: Extinguished', 'Experiment 2: Continued'}, showErrBars, colours, 1, lines); % do excludeCoMFromCert
    
    % expt 2 post-choice CPP by confidence-in-initial-choice
    plotByCondFac([325:326], 4, showErrBars, colours, lines, 0);

    
    
    
    %% topoplots too - plot separately and combine in inkplot
    figure();
    if ~exist('topoplot','file'); eeglab nogui; end
    mapLims = [-3 3];
    PlotTask1Topo(131, 'Experiment 1', colours.cmap, mapLims);
    PlotTask2Topos([1 3 2], {'Experiment 2 Pre-choice'}, colours.cmap, mapLims, 1); % do excludeCoMFromCert
    plotTopoContLow(133, {'Experiment 2 Post-choice'; 'Continued: Maybe'}, colours.cmap, mapLims, 0);

    set(0,'DefaultAxesFontSize',gf); % reset


end





function PlotTask1Traces(subplotInds, respWin, titles, showErrBars, colours, lines)
doStats=1;
    opts.useCSD = 0;
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
    
    ylims = [-0 5];
    
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
    
    opts.useCSD = 0;
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
        labels.certainty = {'maybe (no-CoM)', 'probably (no-CoM)', 'certain (no-CoM)'}; % no-com
    else
        labels.certainty = {'maybe (CoM/no-CoM)', 'probably (CoM/no-CoM)', 'certain (CoM/no-CoM)'};
    end
    %% final figure - prep
    
    % for now, strip unused vars out
    
    behDataByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),behData,'UniformOutput',0); %[pp cond tr]
    cppFiltCond = groupMeans(cppFilt,3,repmat(permute(behData.cond,[1,3,2]),1,size(cppFilt,2)),'dim'); %[pp t cond tr]
    
    faceAlpha = .2;
    rInds = isBetween(eeg.respTimes, respWin);
    
    %% plot
    ylims = [0 5;];
    
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


function PlotTask1Topo(subplotInds, titles, cmap, mapLims)

    opts.useCSD = 0;
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
        'diffFlips', 'eeg','cppTopoAmplWindow','chans5');

    
    
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
    
    opts.useCSD = 0;
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
        'diffFlips', 'eeg','cppTopoAmplWindow','chans5');
    
    
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


function plotByCondFac(subplotInds, iFs, showErrBars, colours, lines, excludeCoMFromCert)

opts.useCSD = 0;
opts.excludeBadPps = 1; % remove pps with <640 good trials?
opts.excludeTooFew = 1; % remove pps with <20 per conf3
opts.excludeByRT = 1; % remove trials outside [100 1500] ms
opts.doFilt = 1; % whether to plot with loPass filtered data
opts.excludeCoMFromCert = excludeCoMFromCert; % remove CoM trials from behData.certainty

opts.outFolder = './Saves';


%% load

opts.saveOpts = {'Volt','CSD'; '', 'ExclCoMFromCert'};
opts.saveName = sprintf('CPPAnalysis_%s_%s.mat', opts.saveOpts{1,opts.useCSD+1}, opts.saveOpts{2, opts.excludeCoMFromCert+1});

optsNames = fieldnames(opts);
data = load(fullfile(opts.outFolder, opts.saveName), optsNames{:}, ...
    'behData', 'cppFilt', 'factors', 'labels',...
    'eeg','amplWindows');



%% check things match

dataNames = fieldnames(data);
data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));

if ~equals(opts, data1)
    warning('loaded data and options do not match');
    keyboard;
end

struct2workspace(data); % unpack

%%
labels.confInR1 = {'certain CoM', 'probably CoM', 'maybe CoM', 'maybe no-CoM', 'probably no-CoM', 'certain no-CoM'};

respWin = [0 1000];
faceAlpha = .1;
yl = [0 5];
rInds = isBetween(eeg.respTimes, respWin);

names.confInR1 = 'confidence-in-initial-choice';
names.certainty = 'final-confidence';
names.conf3 = 'confidence-in-initial-choice: binned';
names.CoM = 'change-of-mind';
names.acc = 'initial accuracy';

%% split by ev-cond

behDataByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),behData,'UniformOutput',0); %[pp cond tr]


%% split trace too (filtered)

cppFiltCond = groupMeans(cppFilt,3,repmat(permute(behData.cond,[1,3,2]),1,size(cppFilt,2)),'dim'); %[pp t cond tr]

nF = 3; % incl CoM

if ~exist('iFs','var') || isempty(iFs)
    iFs = 2:3; % cert & CoM
end
doStats=1;
for i = 1:length(iFs)
    iF = iFs(i);

   resplocked = groupMeans(cppFiltCond,4,repmat(permute(behDataByCond.(factors{iF}),[1,4,2,3]),1,size(cppFilt,2))); %[pp t cond fac]

    for iC=1:2
        subplot(subplotInds((i-1)*2+iC));
        set(gca,'ColorOrder',colours.(factors{iF}),'nextplot','replacechildren');

        h = errorBarPlot(sq(resplocked(:,rInds,iC,:)),'area',1,'xaxisvalues',eeg.respTimes(rInds));
        if ~showErrBars
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
        for j=1:size(h,1)
            h{j,1}.LineWidth = lines.widths.(factors{iF})(j);
            h{j,1}.LineStyle = lines.styles.(factors{iF}){j};
        end

        % shade the mean window
        hold on;
        fill(col(repmat(amplWindows(2,:),2,1)), [yl(1,:) fliplr(yl(1,:))], 'k','FaceAlpha',faceAlpha,'LineStyle','none')


        xticks(-1000:500:1000); xtickangle(0);
        xline(0);yline(0);
        xlim(respWin); %ylim(ylims(2,:));
        xlabel('time from initial RT (ms)');

        if iC==1; ylabel('CPP \muV/m^2'); end
        ylim(yl(1,:));

        if i==1; title([labels.cond(iC); names.(factors{iF})]);
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

        if iC==1
            legend([h{:,1}], [labels.(factors{iF}) labels.diffNames{iF}], 'Location','Best', 'FontSize',10,'Box','off');
        end

    end


end

end

function plotTopoContLow(subplotInds, titles, cmap, mapLims, excludeCoMFromCert)

opts.useCSD = 0;
opts.excludeBadPps = 1; % remove pps with <640 good trials?
opts.excludeTooFew = 1; % remove pps with <20 per conf3
opts.excludeByRT = 1; % remove trials outside [100 1500] ms
opts.doFilt = 1; % whether to plot with loPass filtered data
opts.excludeCoMFromCert = excludeCoMFromCert; % remove CoM trials from behData.certainty

opts.outFolder = './Saves';


%% load

opts.saveOpts = {'Volt','CSD'; '', 'ExclCoMFromCert'};
opts.saveName = sprintf('CPPAnalysis_%s_%s.mat', opts.saveOpts{1,opts.useCSD+1}, opts.saveOpts{2, opts.excludeCoMFromCert+1});

optsNames = fieldnames(opts);
data = load(fullfile(opts.outFolder, opts.saveName), optsNames{:}, ...
    'behData', 'factors', 'labels',...
    'eeg','cppTopoAmplWindow','chans5');



%% check things match

dataNames = fieldnames(data);
data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));

if ~equals(opts, data1)
    warning('loaded data and options do not match');
    keyboard;
end

struct2workspace(data); % unpack

%%
labels.confInR1 = {'certain CoM', 'probably CoM', 'maybe CoM', 'maybe no-CoM', 'probably no-CoM', 'certain no-CoM'};

% just continued-low-certainty
[~, chans5Inds] = ismember(chans5, eeg.chanNames);

subplot(subplotInds);

iF = 2;
iC=2; % just Continued
    %%%%% topo that
    
behDataByCond.(factors{iF}) = groupMeans(behData.(factors{iF}), 2, behData.cond,'dim'); %[pp cond tr]

respTopoCond = groupMeans(cppTopoAmplWindow(:,:,:,2),2,repmat(behData.cond,1,1,eeg.nChans),'dim'); %[pp cond ch tr]

respTopoCondFac = groupMeans(respTopoCond, 4, repmat(permute(behDataByCond.(factors{iF}), [1,2,4,3]),1,1,eeg.nChans));
%     respTopoCondFac = diff(respTopoCondFac(:,:,:,[1 end]),[],4) .* diffFlips(iF);


% just plot low-certainty
topoplot(sq(nanmean(nanmean(respTopoCondFac(:,iC,:,1),4),1)),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', mapLims,...
    'emarker2',{chans5Inds,'.','k',20,1});
t = title(titles);
colorbar;
    
end