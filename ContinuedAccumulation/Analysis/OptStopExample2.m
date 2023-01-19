function OptStopExample2()
% do a toy simulation to illustrate that early stopping + decay can give
% lower final CPP amplitudes for high-confidence
%
% idea is conf-bound hitting = stop accumulating + decay away
% while signal is on average rising until this
% and so early conf-crossing = greater decay
%

clc; %clear; %close all;
gf = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',14);
%% do each one separately?

tLims = [0 1.5]; % post-dec time
dt = .002;
t1 = tLims(1):dt:tLims(2);   % Full time range to simulate. 
nT1 = length(t1);
t2Lims = [0 1];
t2 = t2Lims(1):dt:t2Lims(2);   % Full time range to simulate. 
nT2 = length(t2);

c = 0.1; % confidence-threshold
c2 = -c; % certain CoM thresh

z = .0; % starting point of post-choice accum

gamma = 4; % decay per sec

rng(1);

% have noise?
s = 0.0; % noise SD
n1 = randn(3,nT2).*sqrt(s).*dt;
n2 = randn(3,nT2).*sqrt(s).*dt; % indep
n3 = randn(3,nT2).*sqrt(s).*dt;
n4 = randn(3,nT2).*sqrt(s).*dt; % indep

slopesMaybe = -z + [.06 .07 .085]'; % per sec
slopesCertain = -z + [.12 .7 .8]'; % per sec

% make negative, and a bit weaker than no-CoM
slopesMaybeCoM = z -[.05 .068 .075]'; % per sec
slopesCertainCoM = z - [.115 .3 .6]'; % still has to hit -c


% labels.bounds = {'confidence boundary','initial decision boundary', 'confidence boundary'};
labels.bounds = {'-c','0', 'confidence threshold'};
labels.certainty = {'maybe','certain'};
labels.confidence =  {'Certain CoM','maybe CoM', 'maybe no-CoM', 'certain no-CoM'};
cols = [  0    0.4470    0.8; 0.8500    0.3250    0.0980; ]; % maybe; certain

%% maybe trials

dvMaybe = slopesMaybe.*t2 + z + cumsum(n1,2);

%% certain trials

dvCertain = slopesCertain.*t2 + z + cumsum(n2,2);
% apply decay - linear for now
for i = 1:3
    crossInd = find(dvCertain(i,:) >= c, 1, 'first');
    if isempty(crossInd); continue; end
    % linear
%     y2(i,crossInd+1:end) = gamma.*t2(1:(nT2-crossInd)) + y2(i,crossInd) + n2(i,crossInd+1:end)-n2(i,crossInd);

    % exponential
    for j = 1:(nT2-crossInd)
        dvCertain(i,crossInd+j) = dvCertain(i,crossInd+j-1) - dvCertain(i,crossInd+j-1) .*gamma.*dt + n2(i,crossInd+j); 
    end
end



%% now simulate CoM + no-CoM too

dvMaybeCoM = slopesMaybeCoM.*t2 + z - cumsum(n3,2);

%% certain CoM

dvCertainCoM = slopesCertainCoM.*t2 + z + cumsum(n4,2);
% apply decay - linear for now
for i = 1:3
    crossInd = find(dvCertainCoM(i,:) <= c2, 1, 'first');
    if isempty(crossInd); continue; end
    % linear
%     y4(i,crossInd+1:end) = gamma.*t2(1:(nT2-crossInd)) + y4(i,crossInd) + n2(i,crossInd+1:end)-n2(i,crossInd);

    % exponential
    for j = 1:(nT2-crossInd)
        dvCertainCoM(i,crossInd+j) = dvCertainCoM(i,crossInd+j-1) - dvCertainCoM(i,crossInd+j-1) .*gamma.*dt + n2(i,crossInd+j); 
    end
end


%% plot all 4 on one

allY = cat(3, dvCertainCoM, dvMaybeCoM, dvMaybe, dvCertain); % [certain CoM, maybe CoM, maybe no-CoM, certain no-CoM

nCols = 4;
cols2 = flipud(crameri('roma',nCols));

figure();
% subplot(2,1,1);
set(gca,'ColorOrder', cols2,'nextplot','replacechildren');
h = plot(t2, sq(nanmean(abs(allY),1)), '-','LineWidth', 4);
yline(c, ':k');
xlabel('time from initial decision (ms)')
xticks([0 .5 1]); xticklabels(xticks*1000);
ylabel('CPP = absolute decision variable')
% ylim([0 c+.02]);
legend(h, labels.confidence,'Location','Best','AutoUpdate','off');
yticks([-c z c]); yticklabels(labels.bounds); ytickangle(0);

% super the others?
hold on;
for i = 1:4
    h1 = plot(t2, abs(allY(:,:,i)), 'Color', [cols2(i,:) 1]);
    for j = 1:size(h1,1)
%         h1(j).Color(4) = .2 / ceil(j/2) + .1;
%         h1(j).LineWidth =  1-  .5/ceil(j/2);
    end
end

%%

set(0,'DefaultAxesFontSize',gf); % reset

end