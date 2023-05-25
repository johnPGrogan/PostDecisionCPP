function PlotNewFig5()
% actually Fig 4 in paper
% Plot dSSVEP and beta power figures
% Requires:
% ContinuedAccumulation/Analysis/Saves/PhaseAnalysis_CSD_ExclRT.mat
% ContinuedAccumulation/Analysis/Saves/MuBetaAnalysis_CSD__.mat
% ContinuedAccumulation/Analysis/Saves/GetPhaseTaggingTopo_CSD.mat
% ContinuedAccumulation/Analysis/Saves/ExtractEpochs.mat
% ContinuedAccumulation/Analysis/Saves/BehDataLoad.mat
% ContinuedAccumulation/Analysis/Saves/FlagArtefacts3.mat
% ContinuedAccumulation/Analysis/Saves/DoFFTCSD.mat

gf = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',14);

respWin = [-100 1000; -500 1000];
showErrBars = 1; % show error bars in phase figures?
faceAlpha = .1;

colours.certainty = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];
% colours.certainty = [0.3731    0.8864    0.7382; 0.6115    0.5872    0.1107; 0.5505    0.0068    0.4520; .2 .2 .2];
colours.cond = [1 0 0; 0 0 1; .2 .2 .2];
colours.CoM = [0 0 1; 1 0 0; .2 .2 .2]; %[0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2]; % CoM
colours.conf3 = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];%colours.certainty;
% colours.conf3 = [0.0937    0.3636    0.6709; 0.6796    0.8444    0.6806; 0.5829    0.2780    0.0740; .2 .2 .2]; % av within confinr1 pairs
% colours.conf3 =  [flipud(crameri('roma',3)); .2 .2 .2];
colours.confInR1 =  [flipud(crameri('roma',6)); .2 .2 .2];
colours.cmap = crameri('vik');

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
PlotPhaseFigs(respWin(1,:), showErrBars, faceAlpha, colours, lines, [4 5 2]);


% plot mubeta by cond
figure();
PlotBetaFigs(respWin(2,:), showErrBars, faceAlpha, colours, lines, [4 5 2]);

%% plot topographies - phase tagging - average over all continued good trials
% plots the beta coeff from regression of confInR1 in each acc separately

figure();
PlotPhaseTopo([1 1 1], [800 900 1000], colours, 800); % just mean across all conds/trials
% topoWin = [700 1000]; %ms,  end of 100ms windows to plot mean of , i.e.[400 500] = 300:500ms
% mapLim = 0.3; % from zero to this
% subplotInds = [321 322; 323 324; 325 326;]; % just pass 1st row to only do grand mean, others do CoM effect in err/corr
% PlotPhaseTopoBeta(subplotInds, topoWin, colours, mapLim, [4 5 2]);% plot per factors


%% beta topographies - 2nd resp, grand-mean

topoWin = [800 900 1000]; % endpt (ms) of windows to average over
figure();
% subplot(2,2,4);% figure();
PlotBetaTopo(topoWin, colours); % mean across all conds/trials
% PlotBetaTopoRegs(topoWin, colours, [321 322; 323 324; 325 326;], [4 5 2]); % plot per factors

    set(0,'DefaultAxesFontSize',gf); % reset

end


function PlotPhaseFigs(respWin, showErrBars, faceAlpha, colours, lines, iFs)

    opts.useCSD = 1;
    opts.excludeByRT = 1; % remove trials outside [100 1500] ms
    opts.excludeBadPps = 1;
    opts.excludeTooFew = 1;
    
    opts.outFolder = './Saves';
    
    
    %% load
    
    opts.saveOpts = {'Volt','CSD'; '', 'ExclRT'};
    opts.saveName = sprintf('PhaseAnalysis_%s_%s.mat', opts.saveOpts{1,opts.useCSD+1}, opts.saveOpts{2, opts.excludeByRT+1});
    
    optsNames = fieldnames(opts);
    data = load(fullfile(opts.outFolder, opts.saveName), optsNames{:}, ...
        'phase','phaseNames','phaseTimes','nPhase', 'windows','nTimes',...
        'behData','factors','nF', 'labels', 'diffFlips');
    
    
    %% check things match
    
    dataNames = fieldnames(data);
    data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));
    
    if ~equals(opts, data1)
        warning('loaded data and options do not match');
        keyboard;
    end

    %% move into workspace
    struct2workspace(data)

    windows.resp = [700 1000];
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
    phaseByAcc = structfun(@(x) groupMeans(x, 3, repmat(permute(behData.acc,[1,3,2]),1,size(x,2)), 'dim'), phase, 'UniformOutput',0); %[pp t acc tr]
    
    
    %% final figure - acc*fac
    
    % just post-response? trying to answer the question, do people stop paying
    % attention when they are highly confident?
    
    ylims = [-200 1100; 500 1000];
    iPh = find(strcmp(phaseNames, 'resp')); % just post-resp
    nF=2;

    %% just fac

%     for iF = 3 % +1
%         phaseByFac = groupMeans(phase.(phaseNames{iPh}), 3, repmat(permute(behData.(factors{iF+1}),[1,3,2]),1,nTimes.(phaseNames{iPh}))); %[pp t fac]   
%         
%         subplot(2,3,1);%(iF-1)*3+1);
% 
%         set(gca,'ColorOrder',colours.(factors{iF+1}),'nextplot','replacechildren');
%         h = errorBarPlot(phaseByFac, 'area',1, 'xaxisvalues', phaseTimes.(phaseNames{iPh}));
%         if ~showErrBars
%             for i=1:size(h,1), h{i,2}.Visible = 'off'; end
%         end
%         for j=1:size(h,1)
%             h{j,1}.LineWidth = lines.widths.(factors{iF+1})(j);
%             h{j,1}.LineStyle = lines.styles.(factors{iF+1}){j};
%         end
%         hold on;
%         fill(col(repmat(windows.resp,2,1)), [ylims(1,:) fliplr(ylims(1,:))], 'k','FaceAlpha',faceAlpha,'LineStyle','none')
%    
%         xlabel('time from initial RT (ms)');
%         xline(0);
%         yline(0);
%         ylim(ylims(1,:));
%         xlim(respWin);
%         
%         if iF==2; legend([h{:,1}], labels.(factors{iF+1}), 'Location','Best'); end
%         if iF; ylabel('differential SSVEP'); end
%         title('Continued Evidence');
%     
%     end

    %% acc*fac
    nF = length(iFs);
    for i = 1:nF % +1
        iF=iFs(i);
        data = groupMeans(phaseByAcc.(phaseNames{iPh}), 4, repmat(permute(behDataByAcc.(factors{iF}),[1,4,2,3]),1,nTimes.(phaseNames{iPh}))); %[pp t ac fac]
    
        if any(strcmpi(factors{iF},{'certainty','conf3','CoM'}))
            % average over confInR1 first
            data = groupMeans(phaseByAcc.(phaseNames{iPh}), 4, repmat(permute(behDataByAcc.confInR1,[1,4,2,3]),1,nTimes.(phaseNames{iPh}))); %[pp t ac confinr1]
            facByConfInR1 = groupMeans(behDataByAcc.(factors{iF}),3,behDataByAcc.confInR1); %[pp acc confInR1]
            data = groupMeans(data,4,repmat(permute(facByConfInR1,[1,4,2,3]),1,nTimes.(phaseNames{iPh}))); %[pp t cond fac]
        end
        for iA = 1:2
            subplot(nF,2,(i-1)*2+iA);
    
            set(gca,'ColorOrder',colours.(factors{iF}),'nextplot','replacechildren');
            h = errorBarPlot(sq(data(:,:,iA,:)), 'area',1, 'xaxisvalues', phaseTimes.(phaseNames{iPh}));
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
            fill(col(repmat(windows.resp,2,1)), [ylims(1,:) fliplr(ylims(1,:))], 'k','FaceAlpha',faceAlpha,'LineStyle','none')
       
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
            times = -100:100:1000;
%             x = movmean(times,2,2,'endpoints','discard');
            formula = sprintf('amplWin ~ 1 + %s + (1 | pp)', factors{iF}); % formula to run, in this cond
            
            % run one reg to get number of terms
            regTab.amplWin = nanzscore(col(phase.resp(:,1,:))); % use mean
            fit = fitglme(regTab(regTab.acc==(iA-1),:), formula);
            stats = NaN(length(fit.CoefficientNames)-1,length(times)-1,2);

            for iT = 1:length(times)-1
            
                % take mean within window
                regTab.amplWin = nanzscore(col(nanmean(phase.resp(:,isBetween(phaseTimes.resp, times([iT iT+1])),:),2))); % use mean
            
            
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
            pbar(col([stats(1,:,2);stats(1,:,2)])', 'xVals',col([times(1:end-1);times(2:end)]), 'yVal',yVal,'plotargs',{'Color','k','LineWidth',3});
%             h1.LineStyle = 'none'; h1.Marker = '.'; h1.MarkerSize = 20;




            if iA==1
                legend([h{:,1}], labels.(factors{iF}), 'Location','Best')
%                 if iF; ylabel('differential SSVEP'); end
            end
            if i==1; title([ {['Initial ' labels.acc{iA}]}; names.(factors{iF})]); 
            else; title(names.(factors{iF})); end
            if iA==1; ylabel({'differential SSVEP';'(a.u.)'}); end

        end
    
    end

%     %% acc*fac means
%     meanPhase.(phaseNames{iPh}) = sq(nanmean(phase.(phaseNames{iPh})(:,isBetween(phaseTimes.(phaseNames{iPh}), windows.(phaseNames{iPh})),:),2));
% 
%     meanPhaseByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), meanPhase, 'UniformOutput',0);
%     for i = 1:length(iFs)
%         iF=iFs(i);
%         data = groupMeans(meanPhaseByAcc.(phaseNames{iPh}), 3, behDataByAcc.(factors{iF})); %[pp ac fac]
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
%             ylabel('mean differential SSVEP');
%             if i==1; title('mean within grey window'); end
% %             legend(labels.(factors{iF}),'Location','Best');
%             
% %         end
%     
%     end

%     SuperTitle('Continued Post-Decision Evidence');
end

function PlotBetaFigs(respWin,showErrBars, faceAlpha,colours, lines, iFs)

    
    %% set options
    
    opts.useCSD = 1;
    opts.excludeBadPps = 1; % remove pps with <640 good trials?
    opts.excludeTooFew = 1; % remove pps with <20 per conf3
    opts.excludeByRT = 1; % remove trials outside [100 1500] ms
    opts.excludeCoMFromCert = 0; % remove CoM trials from behData.certainty
    opts.doBaseline = 0; % baseline to lsat window before evonset, log
    
    opts.outFolder = './Saves';
       
    %% load
    
    opts.saveOpts = {'Volt','CSD'; '', 'ExclCoMFromCert'; '', 'BS'};
    opts.saveName = sprintf('MuBetaAnalysis_%s_%s_%s.mat', opts.saveOpts{1,opts.useCSD+1}, opts.saveOpts{2, opts.excludeCoMFromCert+1}, opts.saveOpts{3, opts.doBaseline+1});
    
    optsNames = fieldnames(opts);
    data = load(fullfile(opts.outFolder, opts.saveName), optsNames{:}, ...
        'behData', 'betas', 'betaVars','betaVarNames','factors','nF', 'labels',...
        'f','amplWindows');
    
    
    %% check things match
    
    dataNames = fieldnames(data);
    data1 = rmfield(data, dataNames(~ismember(dataNames, optsNames)));
    
    if ~equals(opts, data1)
        warning('loaded data and options do not match');
        keyboard;
    end
    
    %% move into workspace
    struct2workspace(data)
    
    factors = {'acc','certainty','CoM','confInR1','conf3'};
    %% post
    labels.hemispheres = {'ipsi','contra','contra-ipsi'};
    labels.confInR1 = {'certain CoM', 'probably CoM', 'maybe CoM', 'maybe no-CoM', 'probably no-CoM', 'certain no-CoM'};
    labels.conf3 = {'certain/probably CoM', 'maybe CoM/no-CoM', 'probably/certain no-CoM'};
    labels.certainty = {'maybe CoM/no-CoM', 'probably CoM/no-CoM', 'certain CoM/no-CoM'};
    names.confInR1 = 'confidence-in-initial-choice';
    names.certainty = 'final-certainty';
    names.conf3 = 'confidence-in-initial-choice: binned';
    names.CoM = 'change-of-mind';
%     % invert CoM so that change comes first
%     behData.CoM = 2 - behData.CoM;
%     labels.CoM = flip(labels.CoM);

    % invert CoM trials so that it reflects final response
    betas1 = betas;
    betas1 = reshape(betas,[],3);
    c = col(repmat(permute(behData.CoM==1,[1,3,2]),1,length(f.respWindows),1)); % in same shape
    betas1(c,[1 2]) = betas1(c,[2 1]); % swap hemis
    betas1 = reshape(betas1, size(betas));
    betas1(:,:,:,3) = diff(betas1(:,:,:,1:2),[],4); % contra-ipsi
    betas = betas1;
    clear betas1;

    rInds = isBetween(f.respWindows,respWin);
    nT = sum(rInds);
% %     b=behData;
    behDataByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),behData,'UniformOutput',0); %[pp cond tr]
    nF = 3;
    for i = 1:3
        iF = iFs(i);
    % split by this factor
% %     behData = b;
% %     behData.(factors{iF})(behData.CoM==(3-i)) = NaN; % only CoM for now
    facByCond = groupMeans(behData.(factors{iF}), 2, behData.cond,'dim');
    facByCond = repmat(permute(facByCond,[1,4,2,5,3]),1,nT,1,3,1);
   
    resplocked = groupMeans(betas(:,rInds,:,:),3,repmat(permute(behData.cond,[1,3,2]),1,nT,1,3),'dim'); %[pp t cond 3 tr]
    resplocked1 = groupMeans(resplocked,5,facByCond); %[pp t cond 3 fac]

    if any(strcmpi(factors{iF},{'certainty','conf3','CoM'}))
        % split into confInR1, then average those into pairs/triples
        facByCond = repmat(permute(behDataByCond.confInR1,[1,4,2,5,3]),1,nT,1,3,1);
        resplocked1 = groupMeans(resplocked,5,facByCond); %[pp t cond 3 confInR1]
        facByConfInR1 = groupMeans(behDataByCond.(factors{iF}),3,behDataByCond.confInR1); %[pp cond confInR1]
        resplocked1 = groupMeans(resplocked1, 5, repmat(permute(facByConfInR1,[1 4 2 5 3]),1,nT,1,3,1)); %[pp t cond 3 fac acc]
    end
    
    % is already flipped in MuBetaAnalysis now, this will flip back to
    % initial response hand
    % flip CoM
%     resplocked1(:,:,:,:,1:3) = -resplocked1(:,:,:,:,1:3);
    % flip ipsi & contra rather than inverting both
%     resplocked1(:,:,:,[1 2],1:3) = resplocked1(:,:,:,[2 1],1:3);
%     resplocked1(:,:,:,3,1:3) = diff(resplocked1(:,:,:,1:2,1:3),[],4);
    
    ylims = [-1.5 1.5; -1.4 0];
    ylims = repmat([-1.5 1.5],3,1);
    hemi = 3;%4-i; % 1=ipsi, 2=contr, 3=lateralisation
%     ylims = [-1.5 1.5; 5 7.5; 5 7.5];
    
%     c = colours.certainty([3 2 1 1 2 3],:);
%     w = lines.widths.certainty([3 2 1 1 2 3]);
    c = colours.(factors{iF});
    w = lines.widths.(factors{iF});

    %     figure();
    for iC = 1:2
%         subplot(subplotInds(i, iC))
        subplot(nF,2,(i-1)*2+iC);
%         for hemi=1:3
%         subplot(2,3,(iC-1)*3+hemi)
        set(gca,'ColorOrder',c,'nextplot','replacechildren','LineStyleOrder',{'-'});
        % contra-ipsi
        h = errorBarPlot(sq(resplocked1(:,:,iC,hemi,:)),'area',1,'xaxisvalues',f.respWindows(rInds),'plotargs',{'LineWidth',2});

        % turn off error bars entirely?
        if ~showErrBars
            for j=1:size(h,1), h{j,2}.Visible = 'off'; end
        else % make more transparent
            for j=1:size(h,1)
                h{j,2}.FaceAlpha = lines.alphas.(factors{iF});
            end
        end
        for j=1:size(h,1)
            h{j,1}.LineWidth = w(j);
%             h{j,1}.LineStyle = lines.styles.(factors{iV}){j};
        end
        hold on;
        fill(col(repmat(amplWindows(2,:),2,1)), [ylims(1,:) fliplr(ylims(1,:))], 'k','FaceAlpha',faceAlpha,'LineStyle','none')

        xline(0);yline(0);
%         if i==1; title(labels.cond{iC});
%         else title([labels.cond{iC} ', ' labels.CoM{i-1}]);%labels.hemispheres{hemi}]);
%         end
        if i==1; title([labels.cond(iC); names.(factors{iF})]);
        else; title(names.(factors{iF})); end
        ylim(ylims(i,:));
        xlabel('time from initial RT (ms)');
        xlim(respWin); %ylim(ylims(3,:));

        % do stats
        if ~exist('regTab','var')
            regTab = table;
            regTab.pp = nanzscore(col(behData.pp)); % participant for RE
            regTab.cond = col(behData.cond);
            regTab.CoM = col(behData.CoM);
        end
        regTab.(factors{iF}) = nanzscore(col(behData.(factors{iF}))); % add this factor
        times = -400:200:1000;
%             x = movmean(times,2,2,'endpoints','discard');
        formula = sprintf('amplWin ~ 1 + %s + (1 | pp)', factors{iF}); % formula to run, in this cond
        
        % run one reg to get number of terms
        regTab.amplWin = nanzscore(col(betas(:,1,:,3))); % use mean
        fit = fitglme(regTab(regTab.cond==iC,:), formula);
        stats = NaN(length(fit.CoefficientNames)-1,length(times)-1,2);

        for iT = 1:length(times)-1
        
            % take mean within window
%             betaMean = nanmean(betas(:,isBetween(f.respWindows, times([iT iT+1])),:,3),2); % mean within window
            betaMean = FitCPPSlope(betas(:,:,:,3), times([iT iT+1]), f.respWindows); % slope within window
            regTab.amplWin = nanzscore(col(betaMean)); % use mean
        
        
            if sum(~isnan(regTab.amplWin)) > 100 % if at least 100 data
                fit = fitglme(regTab(regTab.cond==iC,:), formula);
                stats(:,iT,1) = fit.Coefficients.Estimate(2:end); % beta
                stats(:,iT,2) = fit.Coefficients.pValue(2:end); % p-value

                % also run in CoM + not
                if strcmp(factors{iF},'confInR1')
                    for iCoM = 1:2
                        fit = fitglme(regTab(regTab.cond==iC & regTab.CoM==(iCoM-1),:), formula);
                        stats(:,iT,1,iCoM+1) = fit.Coefficients.Estimate(2:end); % beta
                        stats(:,iT,2,iCoM+1) = fit.Coefficients.pValue(2:end); % p-value
                    end
                end
            end
        end

        % plot sig timepoints
%             plot(, stats(:,:,2) < alpha, '.k', 'MarkerSize',20)
        hold on;
        yVal = min(ylim); 
        % duplicate to get it from start:end of a single window if
        % needed
        pbar(col([stats(1,:,2,1);stats(1,:,2,1)])', 'xVals',col([times(1:end-1);times(2:end)]), 'yVal',yVal,'plotargs',{'Color','k','LineWidth',3});
%             h1.LineStyle = 'none'; h1.Marker = '.'; h1.MarkerSize = 20;

        if strcmp(factors{iF},'confInR1') % CoM separately
            for iCoM = 1:2
                pbar(col([stats(1,:,2,1+iCoM);stats(1,:,2,1+iCoM)])', 'xVals',col([times(1:end-1);times(2:end)]), 'yVal',yVal*(1-iCoM/20),'plotargs',{'Color',colours.(factors{iF})((iCoM-1)*5+1,:),'LineWidth',3});
            end
        end
        if iC==1
            ylabel('\beta lateralisation (\muV/m^2)');
            hold on;
%             emptyLegend(2, { {'-'}, {':'}}, {'Color','k','LineWidth',2}, labels.CoM, {'Location','Best'});
            legend([h{:,1}],labels.(factors{iF}), 'Location','Best');
        end
    
        end
        
%     end
    end


end


function PlotPhaseTopo(subplotInds, topoWins, colours, mapLim)

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
        topoName = fullfile(outFolder, 'GetPhaseTaggingTopo_CSD.mat');
    else
        topoName = fullfile(outFolder, 'GetPhaseTaggingTopo.mat');
    end
    
    load(topoName, 'projWaveRespMeanWindows', 'respWins');
    projWaveRespMeanWindows = real(projWaveRespMeanWindows); % just the real component
    
           
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
    
    projWaveRespMeanWindows(repmat(permute(toRemove,[1,3,4,2]), 1, size(projWaveRespMeanWindows,2), eeg.nChans)) = NaN; % remove
    %% plot
    
    mapLims = [-1 1] .* mapLim;
    iT = ismember(respWins(2:end), topoWins); % all windows to average over
    meanTopo = nanmean(projWaveRespMeanWindows(:,iT,:,:),2); % av over windows


    visChanNames = {'A22','A23','A24','A25'}; % midline only

    [~, visChanInds] = ismember(visChanNames, eeg.chanNames);
    
    subplot(subplotInds(1,1),subplotInds(1,2),subplotInds(1,3))
    topoplot(sq(nanmean(nanmean(meanTopo,4),1)),...
        eeg.chanlocs, 'colormap',colours.cmap, 'maplimits',mapLims,...
         'electrodes','off', 'emarker2', {visChanInds, '.','k',20,1});
    c = colorbar('Limits', [-mapLim mapLim]); % there are no negatives
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

function PlotBetaTopo(topoWins, colours)

    % set these
    useCSD = 1;
    excludeBadPps = 1; % <640 trials
    excludeTooFew = 1; % remove if <20 trials per conf3
    excludeByRT = 1; % remove outside [100 1500]
    excludeCoMFromCert = 0; % remove CoM trials from behData.certainty
    doBaseline = 0; % baseline to lsat window before evonset, log
    
    
    % load
    outFolder = './Saves';
    load(fullfile(outFolder, 'ExtractEpochs.mat'));
    load(fullfile(outFolder, 'FlagArtefacts3.mat'),'isFlagged','ppsToExclude');
    load(fullfile(outFolder, 'BehDataLoad.mat'),'behData','labels','behTab');
    
    if useCSD
        f = load('./Saves/DoFFTCSD.mat');
        load(fullfile(outFolder, 'GetMuBetaTopoCSD.mat'),'wins','betaTopoWins');
    else
        f = load('./Saves/DoFFT.mat');
        load(fullfile(outFolder, 'GetMuBetaTopo.mat'),'wins','betaTopoWins');
    end
    
    labels.hemispheres = {'ipsi','contra','contra-ipsi'};
    
    
    %% exlude bad?
    
    toRemove = sq(any(isnan(betaTopoWins(:,:,3,10)),2)); % flagged trials or cpp outliers
    
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
    end
    
    if excludeTooFew
    
        % only remove those trials, not the entire person
        behData.certainty(toRemove) = NaN;
        toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.cond);
    
        toRemove(toRemove2) = 1;
    
    end
    
    % actually exclude now
    betaTopoWins(repmat(toRemove,1,1,eeg.nChans,size(betaTopoWins,4))) = NaN;

    %% calculate lat index

    iT = ismember(wins(2:end), topoWins); % all windows to average over
    betaTopoMean = nanmean(betaTopoWins(:,:,:,iT),4); % av over time %[pp tr ch]

    % get contra, ipsi
    betaTopoHemis = groupMeans(betaTopoMean, 2, repmat(behData.respLR2,1,1,eeg.nChans),'dim'); %[pp ip/con chan tr]
    betaTopoLat = sq(diff(nanmean(betaTopoHemis,4),[],2)); % mean over trials, contra - ipsi. [pp chans]
    

    %% plot
    
    chanClusters = {'D19','D28','D18' 'B22','B21','B18'}'; % 

    [~, chanClusterInds] = ismember(chanClusters, eeg.chanNames);
    
    mapLims = [-1.3 1.3];
    topoplot(sq(nanmean(betaTopoLat,1)),...
        eeg.chanlocs, 'colormap',colours.cmap, 'mapLimits', mapLims, ...
         'electrodes','off', 'emarker2', {chanClusterInds, '.','k',20,1});
    c = colorbar; % there are no negatives




end

function PlotBetaTopoRegs(topoWins, colours, subplotInds, iFs)
factors = {'acc','certainty','CoM','confInR1','conf3'};

    % set these
    useCSD = 1;
    excludeBadPps = 1; % <640 trials
    excludeTooFew = 1; % remove if <20 trials per conf3
    excludeByRT = 1; % remove outside [100 1500]
    excludeCoMFromCert = 0; % remove CoM trials from behData.certainty
    doBaseline = 0; % baseline to lsat window before evonset, log
    
    
    % load
    outFolder = './Saves';
    load(fullfile(outFolder, 'ExtractEpochs.mat'));
    load(fullfile(outFolder, 'FlagArtefacts3.mat'),'isFlagged','ppsToExclude');
    load(fullfile(outFolder, 'BehDataLoad.mat'),'behData','labels','behTab');
    
    if useCSD
        f = load('./Saves/DoFFTCSD.mat');
        load(fullfile(outFolder, 'GetMuBetaTopoCSD.mat'),'wins','betaTopoWins');
    else
        f = load('./Saves/DoFFT.mat');
        load(fullfile(outFolder, 'GetMuBetaTopo.mat'),'wins','betaTopoWins');
    end
    
    labels.hemispheres = {'ipsi','contra','contra-ipsi'};
    
    
    %% exlude bad?
    
    toRemove = sq(any(isnan(betaTopoWins(:,:,3,10)),2)); % flagged trials or cpp outliers
    
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
    end
    
    if excludeTooFew
    
        % only remove those trials, not the entire person
        behData.certainty(toRemove) = NaN;
        toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.cond);
    
        toRemove(toRemove2) = 1;
    
    end
    
    % actually exclude now
    betaTopoWins(repmat(toRemove,1,1,eeg.nChans,size(betaTopoWins,4))) = NaN;

    %% calculate lat index

    iT = ismember(wins(2:end), topoWins); % all windows to average over
    betaTopoMean = nanmean(betaTopoWins(:,:,:,iT),4); % av over time %[pp tr ch]

%     % swap hemisphere electrodes over, or calc lat index per trial here
%     % get coords
%     coords = [[eeg.chanlocs.X]; [eeg.chanlocs.Y]; [eeg.chanlocs.Z]]'; % [back->front, right->left, down->up]
%     coords = round(coords, 4); % round for easier matching later
% 
%     isLeft = find(coords(:,2) > .01);
%     isRight = find(coords(:,2) < -.01); % small deviations from zero in midline elecs
%     isMid = find(coords(:,2)==0);
% 
%     % find matching 
%     coords2 = coords; % lateral ones only
%     coords2(:,2) = abs(coords2(:,2)); % reflect aroud midline
%     
%     [~, inds(:,1)] = sortrows(coords2(isLeft,:)); % sort, will be in same order then
%     [~, inds(:,2)] = sortrows(coords2(isRight,:));
% 
%     mirrorInds = [isLeft(inds(:,1)), isRight(inds(:,2))]; % [left chans, right chans]
% 
%     b = reshape(betaTopoMean,[],eeg.nChans); %[ [pp x tr] chans]
% 
%     b1 = b; % copy
%     b1(col(behData.respLR2==2),mirrorInds(:,1)) = b(col(behData.respLR2==2),mirrorInds(:,2)); % swap left to right chans
%     b1(col(behData.respLR2==2),mirrorInds(:,2)) = b(col(behData.respLR2==2),mirrorInds(:,1)); % and right to left
% 
%     betaTopoMean = reshape(b1, size(betaTopoMean)); % [pp tr chans]
    

    % just include respLR2 as covar

%     % get contra, ipsi
%     betaTopoHemis = groupMeans(betaTopoMean, 2, repmat(behData.respLR2,1,1,eeg.nChans),'dim'); %[pp ip/con chan tr]
% 
%     betaTopoLat = sq(diff(betaTopoHemis,[],2)); % mean over trials, contra - ipsi. [pp chans tr]
    

    %% set up for lme

    if ~exist('regTab','var')
        regTab = table;
        regTab.pp = nanzscore(col(behData.pp)); % participant for RE
        regTab.cond = col(behData.cond);
        regTab.respLR2 = nanzscore(col(behData.respLR2));
    end
    iF=4;
    regTab.(factors{iF}) = nanzscore(col(behData.(factors{iF}))); % add this factor
    %             x = movmean(times,2,2,'endpoints','discard');
    formula = sprintf('amplWin ~ 1 + %s*respLR2 + (1 | pp)', factors{iF}); % formula to run, in this cond
    
    % run one reg to get number of terms
    regTab.amplWin = nanzscore(col(betaTopoMean(:,:,1))); % 
    fit = fitglme(regTab(regTab.cond==1,:), formula);
%     ind = strcmp(fit.CoefficientNames,factors{iF});
    ind=3;% 3 is interaction i.e. factor's effect on ipsi/contra
    % pick respLR2 for pure motor, and other one is just factor but across both hands
    %% plot
    
    chanClusters = {'D19','D28','D18' 'B22','B21','B18'}'; % 

    [~, chanClusterInds] = ismember(chanClusters, eeg.chanNames);
    
    mapLims = [-.05 .05];

    stats = NaN(length(fit.CoefficientNames)-1,eeg.nChans,length(iFs),2);
    for i = 1:length(iFs)
        iF = iFs(i);
        regTab.(factors{iF}) = nanzscore(col(behData.(factors{iF}))); % add this factor
        formula = sprintf('amplWin ~ 1 + %s*respLR2 + (1 | pp)', factors{iF}); % update per factor

        for iC = 1:2
            


            for iCh = 1:eeg.nChans
                regTab.amplWin = nanzscore(col(betaTopoMean(:,:,iCh))); % each channel, within post-dec window
                if sum(~isnan(regTab.amplWin)) > 100 % if at least 100 data
                    fit = fitglme(regTab(regTab.cond==iC,:), formula);
                    stats(:,iCh,i,iC) = fit.Coefficients.Estimate(2:end); % beta
    %                 stats(:,iT,2) = fit.Coefficients.pValue(2:end); % p-value
                end
            end

            subplot(subplotInds(i, iC))
            topoplot(stats(ind,:,i,iC),...
                eeg.chanlocs, 'colormap',colours.cmap, 'mapLimits',mapLims, ...
                 'electrodes','off', 'emarker2', {chanClusterInds, '.','k',20,1});
            c = colorbar; % there are no negatives
        end
    end



end

function PlotPhaseTopoBeta(subplotInds, topoWins, colours, mapLim, iFs)
% do regressions at each channel, plot beta

    factors = {'acc','certainty','CoM','confInR1','conf3'};
    

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
        topoName = fullfile(outFolder, 'GetPhaseTaggingTopo_CSD.mat');
    else
        topoName = fullfile(outFolder, 'GetPhaseTaggingTopo.mat');
    end
    
    load(topoName, 'projWaveRespMeanWindows', 'respWins');
    projWaveRespMeanWindows = real(projWaveRespMeanWindows); % just the real component
    
           
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
    
    projWaveRespMeanWindows(repmat(permute(toRemove,[1,3,4,2]), 1, size(projWaveRespMeanWindows,2), eeg.nChans)) = NaN; % remove
    %% plot
    
    mapLims = [-1 1] .* mapLim;
    iT = ismember(respWins(2:end), topoWins); % all windows to average over
    meanTopo = nanmean(projWaveRespMeanWindows(:,iT,:,:),2); % av over windows

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


    visChanNames = {'A22','A23','A24','A25'}; % midline only

    [~, visChanInds] = ismember(visChanNames, eeg.chanNames);
    nF = length(iFs);
    for i = 1:nF
        iF = iFs(i);
        regTab.(factors{iF}) = nanzscore(col(behData.(factors{iF}))); % add this factor
        formula = sprintf('amplWin ~ 1 + %s + (1 | pp)', factors{iF}); % update with factor name

        stats = NaN(length(fit.CoefficientNames)-1,eeg.nChans,2,2);
    
        for iCh = 1:eeg.nChans
    
            % take mean within window
            regTab.amplWin = nanzscore(col(meanTopo(:,1,iCh,:))); % each channel, within post-dec window
    
    
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