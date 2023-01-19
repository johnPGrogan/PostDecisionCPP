% BetaLatStuff
% try to see if beta changes earlier for high certainty
clear; clc; close all
load('./Saves/MuBetaAnalysis_CSD__.mat','betas','behData','f','factors','regTab','labels');

%% find when CoM lateralisation starts to fall?

% invert CoM trials so that it reflects final response
betas1 = betas;
betas1 = reshape(betas,[],3);
c = col(repmat(permute(behData.CoM==1,[1,3,2]),1,length(f.respWindows),1)); % in same shape
betas1(c,[1 2]) = betas1(c,[2 1]); % swap hemis
betas1 = reshape(betas1, size(betas));
betas1(:,:,:,3) = diff(betas1(:,:,:,1:2),[],4); % contra-ipsi
betas = betas1;
clear betas1;


rInds = isBetween(f.respWindows,[-500 1000]);
nT = sum(rInds);
   

iF=4; % confinr1
facByCond = groupMeans(behData.(factors{iF}), 2, behData.cond,'dim');
facByCond = repmat(permute(facByCond,[1,4,2,5,3]),1,nT,1,3,1);
   
resplocked = groupMeans(betas(:,rInds,:,:),3,repmat(permute(behData.cond,[1,3,2]),1,nT,1,3),'dim'); %[pp t cond 3 tr]
resplocked1 = groupMeans(resplocked,5,facByCond,'dim'); %[pp t cond hemi fac tr]

y =  sq(resplocked1(4,:,2,3,1,1:10))'; % lat, continued, certIncorr
x = f.respWindows(rInds);

%% or can just find when it crosses zero?

% use av traces per pp for now
resplockedMean = nanmean(resplocked1,6); 

t = f.respWindows(rInds);
t2 = t(t>=0);

window = [0 1000]; % only in this

iH = 3;
% crossLat = NaN(30,6,2);
% for iPP = 1:30
%     for iC = 1:2
%        for iL = 1:size(resplockedMean,5)
%            tInd = find(resplockedMean(iPP,t>=0,iC,iH,iL) <= 0, 1);
%            if tInd
%             crossLat(iPP,iL,iC) = t2(tInd); % index
%            end
%        end
%     end
% end

%% can i fit a polynomial to each?
% 
% % instead find first difference from value at init resp
% initMeans = nanmean(resplockedMean(:,isBetween(t, [-50 50]),:,:,:),2);
% 
% diffLat = NaN(30,6,2);
% for iPP = 1:30
%     for iC = 1:2
%        for iL = 1:size(resplockedMean,5)
% %            if ~isnan(finalMeans(iPP,1,iC,iH,iL))
% %            p = polyfit(t2,resplockedMean(iPP,t>=0,iC,iH,iL),2); % gives [x^2 x c]
%            % now find turning/max point?
% 
%            % isn't that similar to finding max in trace?
% 
%            % I want to find when the descent begins
% 
%            % or when it first becomes indistinguishable from it's final
%            % value
% 
%            % so take mean in last 100ms
%            % do t-tests to see earliest consec point of difference
%            y = sq(resplocked1(iPP,t>=0,iC,iH,iL,:)); %[t tr]
%            [~,p] = ttest(y', initMeans(iPP,1,iC,iH,iL)); % [1 t]
% 
%            % find beginning of 200ms sig period
%            pReg = findregions(p<.05); % find is sig
%            d = find(diff(pReg,[],2)>=10,1); % length
%            tInd = pReg(d,1); % start of that period (10 samples x 20ms = 200ms)       
% 
% %            tInd = find(p<=.05,1,'first');
%            if tInd
%             diffLat(iPP,iL,iC) = t2(tInd); % find last difference
%            end
%        end
%     end
% end

%% find when becomes n.s. from -1, per pp

t2 = t(t>=0);

endLat = NaN(30,6,2);
for iPP = 1:30
    for iC = 1:2
       for iL = 1:size(resplockedMean,5)
%            if ~isnan(finalMeans(iPP,1,iC,iH,iL))
%            p = polyfit(t2,resplockedMean(iPP,t>=0,iC,iH,iL),2); % gives [x^2 x c]
           % now find turning/max point?

           % isn't that similar to finding max in trace?

           % I want to find when the descent begins

           % or when it first becomes indistinguishable from it's final
           % value

           % so take mean in last 100ms
           % do t-tests to see earliest consec point of difference
           y = sq(resplocked1(iPP,t>=0,iC,iH,iL,:)); %[t tr]
           [~,p] = ttest(y', -1, 'tail', 'both'); % [1 t]
           
           tInd = find(p>.05,1,'first');
           if tInd
               endLat(iPP,iL,iC) = t2(tInd); % find last difference
           else
               if all(p<=.05)
                   endLat(iPP,iL,iC) = NaN; %t2(end);
               elseif all(p>.05)
                   endLat(iPP,iL,iC) = NaN; %t2(1);
               end
           end
       end
    end
end

% %% find when becomes n.s. from zero for 200ms
% 
% t2 = t(t>=0);
% 
% endLat = NaN(30,6,2);
% for iPP = 1:30
%     for iC = 1:2
%        for iL = 1:size(resplockedMean,5)
% %            if ~isnan(finalMeans(iPP,1,iC,iH,iL))
% %            p = polyfit(t2,resplockedMean(iPP,t>=0,iC,iH,iL),2); % gives [x^2 x c]
%            % now find turning/max point?
% 
%            % isn't that similar to finding max in trace?
% 
%            % I want to find when the descent begins
% 
%            % or when it first becomes indistinguishable from it's final
%            % value
% 
%            % so take mean in last 100ms
%            % do t-tests to see earliest consec point of difference
%            y = sq(resplocked1(iPP,t>=0,iC,iH,iL,:)); %[t tr]
%            [~,p] = ttest(y',0,'tail','left'); % [1 t]
% 
%            % find beginning of 200ms sig period
%            pReg = findregions(p>.05); % find is sig
%            d = find(diff(pReg,[],2)>=10,1,'last'); % length
%            tInd = pReg(d,1); % start of that period (10 samples x 20ms = 200ms)       
% %            tInd = find(p>.05,1,'first');
%            if tInd
%                endLat(iPP,iL,iC) = t2(tInd); % find last difference
%            else
%                if all(p<=.05)
%                    endLat(iPP,iL,iC) = NaN; %t2(end);
%                elseif all(p>.05)
%                    endLat(iPP,iL,iC) = NaN; %t2(1);
%                end
%            end
%        end
%     end
% end

%% find n.s. times across pps (vs -1)

t2 = t(t>=0);

endLat = NaN(3,2);
for iC = 1:2
    for iL = 1:size(endLat,1)

       % so take mean in last 100ms
       % do t-tests to see earliest consec point of difference
       y = sq(resplockedMean(:,t>=0,iC,iH,iL)); %[t tr]
%        [~,p] = ttest(y, -1, 'tail','both'); % [1 t]

       v(iL,iC) = -0.9881; % mean from figure
       % or vs mean final value?
%        v(iL,iC) = nanmean(nanmean(y(:,t(t>=0)>=900),2),1); % mean final value for this level
%        v(iL,iC) = nanmean(nanmean(nanmean(resplockedMean(:,isBetween(t, [-150 -50]),iC,iH,4:6),5),2),1); % mean initResp level
       [~,p] = ttest(y, v(iL,iC), 'tail','both'); % [1 t] 

       pReg = findregions(p>.05); % find is sig
       d = find(diff(pReg,[],2)>=1,1,'first'); % length
       tInd = pReg(d,1); % start of that period (10 samples x 20ms = 200ms)       
%            tInd = find(p>.05,1,'first');
       if tInd
           endLat(iL,iC) = t2(tInd); % find last difference
       else
           if all(p<=.05)
               endLat(iL,iC) = NaN; %t2(end);
           elseif all(p>.05)
               endLat(iL,iC) = NaN; %t2(1);
           end
       end
   end
end

figure();
colours.confInR1 =  [flipud(crameri('roma',6)); .2 .2 .2];
set(gca,'ColorOrder',colours.confInR1,'nextplot','replacechildren');
h = errorBarPlot(sq(resplockedMean(:,:,2,3,1:3)),'area',1,'xaxisvalues',t);
xline(0);yline(0);
xlabel('time from init RT');
ylabel('beta lateralisation');
if length(unique(endLat(1:3,2))) < 3
    endLat(1:3,2) = endLat(1:3,2) + [0 -5 5]';
end
h1 = xlines(endLat(1:3,2),'--','LineWidth',2);
for i = 1:3
    h1(i).Color = colours.confInR1(i,:);
end
if length(unique(v(~isnan(v))))>1
    h2 = ylines(v(1:3,2),':','LineWidth',1.5);
    for i = 1:3
        h2(i).Color = colours.confInR1(i,:);
    end
else
    yline(unique(v(~isnan(v))),':k');
end
legend(([h{:,1}]),labels.confInR1,'Location','Best');
%% just find <= -1 in each trace, for fitlgme

endLatTr = sq(apply(2,@(X) max([NaN find(X<=-1,1,'first')]), betas(:,f.respWindows>0,:,3))); %[pp tr]

endLatCondConf = sq(apply(2,@(X) max([NaN find(X<=-1,1,'first')]), resplockedMean(:,t>0,:,3,:))); %[pp tr]

behDataByCond.confInR1= groupMeans(behData.confInR1,2,behData.cond,'dim');
latByCond = groupMeans(endLatTr,2,behData.cond,'dim');
latByCondConf = groupMeans(latByCond,3,behDataByCond.confInR1,'dim');
regTab.lat = nanzscore(col(endLatTr));
fitglme(regTab(regTab.cond>=0 & regTab.CoM>=0,:),'lat ~ 1 + certainty + (1 | pp)') % isn't sig tho

%% cond plot init RT vs later confidence

behDataByCond = structfun(@(x) groupMeans(x,2,behData.cond,'dim'), behData,'Uni',0);

facs = {'RT', 'certainty'; 'confRT', 'certainty'};

figure();
for i = 1:2
    subplot(1,2,i);
    [~,~,~,~,h] = conditionalPlot(permute(behDataByCond.(facs{i,1}),[3 1 2]),permute(behDataByCond.(facs{i,2}),[3 1 2]));
    xlabel(facs{i,1});
    ylabel(facs{i,2});
    if i==1; legend(h(2:2:end), labels.cond,'Location','Best'); end
end

% also split by CoM

figure();
colours.CoM = [0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2];
for i = 1:2
    
    x = groupMeans(behDataByCond.(facs{i,1}),3,behDataByCond.CoM,'dim');
    y = groupMeans(behDataByCond.(facs{i,2}),3,behDataByCond.CoM,'dim'); %[pp cond com tr]
    
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
%         set(gca,'ColorOrder',colours.CoM,'nextplot','replacechildren');
        [~,~,~,~,h] = conditionalPlot(permute(x(:,:,j,:),[4 1 2 3]),permute(y(:,:,j,:),[4 1 2 3]));
        xlabel(facs{i,1});
        ylabel(facs{i,2});
        if i==1; legend(h(2:2:end), labels.cond,'Location','Best'); end
        title(labels.CoM(j));
        axis([0 1500 1.5 2.7]);
    end
end
%% cond plot beta vs cpp
% 
% figure;
% betaByCond = groupMeans(betaVars.betaAmplPost(:,:,3),2,behData.cond,'dim');
% cppByCond = groupMeans(cppVars.amplPost,2,behData.cond,'dim');
% conditionalPlot(permute(betaByCond,[3 1 2]),permute(cppByCond,[3 1 2]));


%% can i just use slopes?

% if use whole window?

betaVars.betaSlopePost  = FitCPPSlope(betas, [0 1000], f.respWindows); % get slopes
regTab.slopePost = nanzscore(col(betaVars.betaSlopePost(:,:,3)));
fitglme(regTab(regTab.cond>=0 & regTab.CoM>=0,:),'slopePost ~ 1 + certainty + (1 | pp)') % isn't sig tho


%% does rate of change slow down?

dResplocked1 = diff(resplocked1,[],2); % [pp t cond hemi conf tr]
x = f.respWindows(rInds);

colours.certainty = [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];
% plot
clf;
set(gca,'ColorOrder',flipud(crameri('roma',6)),'nextplot','replacechildren');
h = errorBarPlot(sq(nanmean(abs(dResplocked1(:,:,2,3,:,:)),6)),'area',1,'xaxisvalues',x(2:end));


%% print means into csv for jasp

% this is using inverted betas
amplWindows = [-300 0; 700 1000]; % ms
betaVars.betaAmplPost = sq(nanmean(betas(:,isBetween(f.respWindows,amplWindows(2,:)),:,:),2));
amplByCond = groupMeans(betaVars.betaAmplPost(:,:,3),2,behData.cond,'dim');

factor = 'confInR1';
data = groupMeans(amplByCond,3,behDataByCond.(factor)); %[pp cond factor]

data = reshape(permute(data,[1 3 2]), 30, []); %[pp [cond then cert]
t = array2table(data,'VariableNames',col(strcat(repmat(labels.cond',1,length(labels.(factor))),{'_'} ,repmat(labels.(factor),2,1))')');
t.pp = (1:height(t))';

writetable(t, './anovaJaspBetaConfInR1.csv');
