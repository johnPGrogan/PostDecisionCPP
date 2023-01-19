function PlotNewFig3(doSuppl)
% Plot second-order CPP traces, split by post-dec evidence, confidence, CoM
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

%% check things match

dataNames = fieldnames(data);
data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));

if ~equals(opts, data1)
    warning('loaded data and options do not match');
    keyboard;
end

%%

data.respWin = [0 1000];
data.showErrBars = 1; % error bars on traces?
data.faceAlpha = .1;
data.yl = [-20 40];
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
data.names.certainty = 'confidence-in-final-choice';
data.names.conf3 = 'confidence-in-initial-choice: binned';
data.names.CoM = 'change-of-mind';
data.names.acc = 'initial accuracy';

%% split by ev-cond

data.behDataByCond = structfun(@(x) groupMeans(x, 2, data.behData.cond,'dim'),data.behData,'UniformOutput',0); %[pp cond tr]


%% split trace too (filtered)

data.cppFiltCond = groupMeans(data.cppFilt,3,repmat(permute(data.behData.cond,[1,3,2]),1,size(data.cppFilt,2)),'dim'); %[pp t cond tr]

data.nF = 3; % incl CoM

%% plot
if exist('doSuppl','var') && doSuppl % do supplementary figure

    % Continued vs Extinguished
    figure();
    plotByCond([111],data)
    
    % supplementary figure - split by initial accuracy
    figure();
    iFs = 1;
    plotByCondFac([121,122], data, iFs);

else

    
    figure();
    plotByCondFac([321:326], data, [4 5 2]);



    %% do topoplots, together or separately
    
    if ~exist('topoplot','file'); eeglab nogui; end
    
    figure();
    plotTopoContLow(data);

end
%%  
set(0,'DefaultAxesFontSize',gf); % reset

end

function plotByCond(subplotInds, data)
%% plot total cpp split by cond, with topoplots below

struct2workspace(data)

subplot(subplotInds)
set(gca,'ColorOrder', colours.cond, 'nextplot','replacechildren');

resplocked = groupMeans(cppFilt,3, repmat(permute(behData.cond,[1,3,2]),1,length(eeg.respTimes)));
resplocked(:,:,3) = diff(resplocked,[],3); % diff wave

h = errorBarPlot(resplocked(:,rInds,:),'area',1,'xaxisvalues',eeg.respTimes(rInds));
if ~showErrBars
    for j=1:size(h,1)-1
        h{j,2}.Visible='off';
        h{j,1}.LineWidth = 2;
    end % remove error bars
end
for j=1:size(h,1)
    h{j,1}.LineWidth = data.lines.widths.cond(j);
    h{j,1}.LineStyle = data.lines.styles.cond{j};
end
hold on;
xline(0,':k');
yline(0, ':k');
yl = ylim;
fill(col(repmat(amplWindows(2,:),2,1)), [yl fliplr(yl)], 'k','FaceAlpha',faceAlpha,'LineStyle','none')
ylim(yl);
ylabel('CPP \muV/m^2');
xlabel('time from initial response (ms)');
box off;
legend([h{:,1}], [labels.cond, 'difference wave'],'Location','NorthEast');
xlim(respWin);
title('grand mean');

end

function plotByCondFac(subplotInds, data, iFs)
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


function plotToposCond(subplotInds, data)
struct2workspace(data);

mapLims2 = [-20 20; -20 20; -10 10];


cppMeanTopoByCond = groupMeans(cppMeanTopo,2,repmat(behData.cond,1,1,eeg.nChans,length(wins)-1)); %[pp cond chan time]
cppMeanTopoByCond(:,3,:,:) = diff(cppMeanTopoByCond(:,1:2,:,:),[],2); % diff



% topoplots within some windows - win 100 is 0:100, so below will go from
% -100 before first window up to second. I will add 1 to each window limit
% to make it match what is actually written here
windows = [100 300; 400 500; 600 800; 900 1000];
nWin = size(windows,1);

% put chans5 electrodes on?
[~, chans5Inds] = ismember(chans5, eeg.chanNames);
for i = 1:nWin
    for j = 3
        subplot(subplotInds(i));
        winInds = isBetween(wins(2:end), windows(i,:) + [1 1]); % add 
        topoplot(sq(nanmean(nanmean(cppMeanTopoByCond(:,j,:,winInds),4),1)),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap,'mapLimits',mapLims2(j,:),...
            'emarker2',{chans5Inds,'.','k',10,1});
        if j==3
            t = {'Cont-Int'};
        else
            t = labels.cond(j);
        end
        if j
            t{2} = sprintf('%d:%dms', windows(i,1), windows(i,2));
        end
        title(flip(t));

        if i==2 
            c = colorbar('Location','North');%,'Position',[0.4767 0.15+(j==1)*.2 0.0786 0.0320]);
            c.Position = c.Position + [0 -.3 0 0];
        end

    end
end

end

function plotToposCondFac(subplotInds, data, iFs)
struct2workspace(data);
[~, chans5Inds] = ismember(chans5, eeg.chanNames);
if ~exist('iFs','var') || isempty(iFs)
    iFs = 2:3; % cert & CoM
end
for i = 1:length(iFs)
    iF=iFs(i);
    %%%%% topo that
    
    respTopoCond = groupMeans(cppTopoAmplWindow(:,:,:,2),2,repmat(behData.cond,1,1,eeg.nChans),'dim'); %[pp cond ch tr]

    respTopoCondFac = groupMeans(respTopoCond, 4, repmat(permute(behDataByCond.(factors{iF}), [1,2,4,3]),1,1,eeg.nChans));
    respTopoCondFac = diff(respTopoCondFac(:,:,:,[1 end]),[],4) .* diffFlips(iF); % diff wave

    for iC=1:2
        
        subplot(subplotInds((i-1)*2+iC));
        topoplot(sq(nanmean(respTopoCondFac(:,iC,:),1)),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', [-10 10],...
            'emarker2',{chans5Inds,'.','k',20,1});
        if iC==2
            t = title(labels.diffNames{iF});
            t.Position = [-.9 .5 .95];
        end
        
    end
    
end
end

function plotTopoContLow(data)
% just continued-low-certainty
struct2workspace(data);
[~, chans5Inds] = ismember(chans5, eeg.chanNames);


for iF = 2
    %%%%% topo that
    
    respTopoCond = groupMeans(cppTopoAmplWindow(:,:,:,2),2,repmat(behData.cond,1,1,eeg.nChans),'dim'); %[pp cond ch tr]

    respTopoCondFac = groupMeans(respTopoCond, 4, repmat(permute(behDataByCond.(factors{iF}), [1,2,4,3]),1,1,eeg.nChans));
%     respTopoCondFac = diff(respTopoCondFac(:,:,:,[1 end]),[],4) .* diffFlips(iF);

    for iC=2
        
        topoplot(sq(nanmean(nanmean(respTopoCondFac(:,iC,:,1),4),1)),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', [-20 20],...
            'emarker2',{chans5Inds,'.','k',20,1});
        if iC==2
            t = title([labels.cond{iC} ', ' labels.(factors{iF}){1}]);
%             t.Position = [-.9 .5 .95];
        end
        colorbar;
    end
    
end
end

function plotTopoDelayEffect(data)
% just extinguished - difference between 1000ms and 0ms
struct2workspace(data);
[~, chans5Inds] = ismember(chans5, eeg.chanNames);

cppMeanTopoByCond = groupMeans(cppMeanTopo,2,repmat(behData.cond,1,1,128,20)); %[pp cond chan time]

inds = any(wins(2:end) == [100;1000]); % 0-100ms and 900-1000ms
cppMeanTopoByCond = diff(cppMeanTopoByCond(:,:,:,inds),[],4); % 1000ms - 100ms

figure();
for i = 1:2
    subplot(1,2,i);
    topoplot(sq(nanmean(cppMeanTopoByCond(:,i,:),1)),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap,'mapLimits',[-15 15],...
            'emarker',{'.','k',10,1}, 'emarker2', {chans5Inds, '.','k',20,1});
    title(labels.cond{i});
end
colorbar
SuperTitle('change in CSD across delay period');

end

function plotToposCondFacRegress(subplotInds, data, iFs)
% regress each channel, topoplot the betas for the factor
struct2workspace(data);
[~, chans5Inds] = ismember(chans5, eeg.chanNames);
if ~exist('iFs','var') || isempty(iFs)
    iFs = 2:3; % cert & CoM
end
for i = 1:length(iFs)
    iF=iFs(i);
    %%%%% topo that
    
%     respTopoCond = groupMeans(cppTopoAmplWindow(:,:,:,2),2,repmat(behData.cond,1,1,eeg.nChans),'dim'); %[pp cond ch tr]
% 
%     respTopoCondFac = groupMeans(respTopoCond, 4, repmat(permute(behDataByCond.(factors{iF}), [1,2,4,3]),1,1,eeg.nChans));
%     respTopoCondFac = diff(respTopoCondFac(:,:,:,[1 end]),[],4) .* diffFlips(iF); % diff wave

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
    regTab.amplWin = nanzscore(col(cppTopoAmplWindow(:,:,1,2))); % use first
    fit = fitglme(regTab(regTab.cond==1,:), formula);

    for iC=1:2

        
        stats = NaN(length(fit.CoefficientNames)-1,eeg.nChans,1);

        for iT = 1:eeg.nChans

            % take mean within window
            regTab.amplWin = nanzscore(col(cppTopoAmplWindow(:,:,iT,2))); % each channel, within post-dec window


            if sum(~isnan(regTab.amplWin)) > 100 % if at least 100 data
                fit = fitglme(regTab(regTab.cond==iC,:), formula);
                stats(:,iT,1) = fit.Coefficients.Estimate(2:end); % beta
%                 stats(:,iT,2) = fit.Coefficients.pValue(2:end); % p-value
            end
        end

        subplot(subplotInds((i-1)*2+iC));
        topoplot(stats(:,:,1),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap,'mapLimits',[-.05 .05],...
            'emarker2',{chans5Inds,'.','k',20,1});
        if iC==2
            t = title(['\beta coefficient for ' factors{iF}]);
            t.Position = [-.9 .5 .95];
        end
        
    end
    
end
end
