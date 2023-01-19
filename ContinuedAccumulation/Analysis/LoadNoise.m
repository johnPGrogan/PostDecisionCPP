% LoadNoise

clear;
load('./Saves/BehDataLoad.mat');
load('./Saves/FlagArtefacts3.mat', 'isFlagged','ppsToExclude');

noiseData.doNorm = 0;
noiseData.invertIncorrect = 0;
noiseData.ignoreNames = {'thresh','ignoreNames','deltaC','deltaCSigma','doNorm','invertIncorrect','times','nPre'}; % ignore these later

%% for each pp/block, load up raw data

noiseData.pre = NaN(30,1280, 92);
noiseData.post = NaN(30,1280, 60);

for iPP = 1:30
    disp(iPP);
    thisPP = find(cellRegexpi(fileInfo.fNames, fileInfo.ppID(iPP)) == 1);
    thisBlocks = fileInfo.blockNums(thisPP);
    
    for iB = 1:16
        
        fileNames{iB} = fullfile(fileInfo.files(thisPP(thisBlocks(iB))).folder, fileInfo.files(thisPP(thisBlocks(iB))).name);
        
        matName = [fileNames{iB}(1:end-4) '.mat'];
        mat = load(matName);
        
        if length(mat.Confidence.scoretemp)<80
            mat.Confidence.scoretemp(length(mat.Confidence.scoretemp):80) = NaN;
            mat.Confidence.Keypressed{length(mat.Confidence.scoretemp):80} = [];
        end
        
        
        % match by RT first, then confRT
        rt = sq(dataMat(iPP,(iB*80-79):iB*80,7)) / 1000;
        rtConf = sq(dataMat(iPP,(iB*80-79):iB*80,9)) / 1000;
        
        
        % get noise stuff
        ev = NaN(80,92);
        evPost = NaN(80,60); % always 60
        for i = 1:80
            if isfield(mat.par.stim_ev{i}, 'deltaC_presented')
                ev(i,1:size(mat.par.stim_ev{i}.deltaC_presented,2))= [mat.par.stim_ev{i}.deltaC_presented] - mat.par.deltaC(1);  % subtract mean
            end
            
            if isfield(mat.par.stim_ev{i}, 'deltaC_presentedafter')
                evPost(i,:) = [mat.par.stim_ev{i}.deltaC_presentedafter] - mat.par.deltaC(1);
            end
        end
        deltaC(iPP,iB) = mat.par.deltaC(1);
        deltaCSigma(iPP,iB) = mat.par.deltaCsigma;
%         ev = sq(nancat(ev))'; %[tr samples]
%         evPost = sq(nancat(evPost))';
        
        
        matchedTrials = (rt - mat.respT') == 0 | isnan(rt);
        
        if all(matchedTrials)
            noiseData.pre(iPP,(iB*80-79):iB*80,1:size(ev,2)) = ev;
            noiseData.post(iPP,(iB*80-79):iB*80,:) = evPost;
        else
            keyboard;
        end

        %% also get which sitmulus was flickered first

        % get which stim was shown first in ssvep flicker
        firstTilt(iPP,(iB*80-79):iB*80) = mat.par.firstTilt;

    end
end

%% need to discard every 2nd sample as they weren't shown

noiseData.pre(:,:,2:2:end) = [];
noiseData.post(:,:,2:2:end) = [];


%% sum

noiseData.preSum = nansum(noiseData.pre,3);
noiseData.preSum(all(isnan(noiseData.pre),3)) = NaN;

noiseData.postSum = nansum(noiseData.post,3);
noiseData.postSum(all(isnan(noiseData.post),3)) = NaN;

noiseData.nPre = sum(~isnan(noiseData.pre),3);

%% get it into resp-locked epoch, along with times,

preResp = reshape(alignRight(reshape(noiseData.pre, [], size(noiseData.pre,3))), size(noiseData.pre));

noiseData.respLocked = cat(3, preResp, noiseData.post);
noiseData.times = (-90:2:60) * (1000/60);

%% reshape deltaC & deltaCSigma#

% get into [pp tr]
deltaC = repmat(deltaC, 1,1,80);
noiseData.deltaC = reshape(permute(deltaC, [1,3,2]),fileInfo.nPP, fileInfo.maxTr);

deltaCSigma = repmat(deltaCSigma, 1,1,80);
noiseData.deltaCSigma = reshape(permute(deltaCSigma, [1,3,2]),fileInfo.nPP, fileInfo.maxTr);


%% invert to make it evidence towards the 1st chosen option
% i.e. flip incorrect trials

if noiseData.invertIncorrect
    noiseNames = fieldnames(noiseData);
    noiseNames(ismember(noiseNames, noiseData.ignoreNames)) = []; % not the number of frames or sample times
    nNoise = length(noiseNames);
    for i = 1:nNoise
        noiseData.(noiseNames{i}) = noiseData.(noiseNames{i}) .* (2*behData.acc-1);
    end
end
%% shall I normalise them all by deltaC?
% could make it so 0 is equal, +1 is deltaC, -1 is incorrect option
% but then cumsum is always positive
% instead make 0 =deltaC, -1 = -deltaC

if noiseData.doNorm
    noiseNames = fieldnames(noiseData);
    noiseNames(ismember(noiseNames, noiseData.ignoreNames)) = []; % not the number of frames or sample times
    for i = 1:length(noiseNames)
        noiseData.(noiseNames{i}) = (noiseData.(noiseNames{i})) ./ noiseData.deltaC;
    end
    noiseData.thresh = 0;
else
    noiseData.thresh = 0;
end


%% split by pre/sum

noiseData.preSumStrong = double(noiseData.preSum > noiseData.thresh);
noiseData.postSumStrong = double(noiseData.postSum > noiseData.thresh);
noiseData.totalSum = nansum(cat(3, noiseData.preSum,  noiseData.postSum),3); % interrupts = just pre
noiseData.totalSumStrong = double(noiseData.totalSum > noiseData.thresh);


% accumulated noise over frames
noiseData.respLockedCumul = cumsum(noiseData.respLocked,3,'omitnan');
% set nan to nan
noiseData.respLockedCumul(isnan(noiseData.respLocked)) = NaN;


%% get slopes of cumulative noise in each period

noiseData.preSlope = FitCPPSlope(permute(noiseData.respLockedCumul,[1,3,2]), [min(noiseData.times), 0], noiseData.times);
noiseData.postSlope = FitCPPSlope(permute(noiseData.respLockedCumul,[1,3,2]), [1 max(noiseData.times)], noiseData.times);
noiseData.totalSlope = FitCPPSlope(permute(noiseData.respLockedCumul,[1,3,2]), minMax(noiseData.times,2), noiseData.times);

noiseData.preSlopeStrong = double(noiseData.preSlope > noiseData.thresh);
noiseData.postSlopeStrong = double(noiseData.postSlope > noiseData.thresh);
noiseData.totalSlopeStrong = double(noiseData.totalSlope > noiseData.thresh);

%% early post-dec noise only - do same metrics

tLims = [1 500]; % ms
tIndsPost = isBetween(noiseData.times(noiseData.times>0), tLims); % just in post-dec times
tInds = isBetween(noiseData.times, tLims); % in total times

noiseData.postEarly = noiseData.post(:,:,tIndsPost); % just keep early samples
noiseData.postSumEarly = nansum(noiseData.postEarly,3);
noiseData.postSumEarly(all(isnan(noiseData.postEarly),3)) = NaN;
noiseData.postSumEarlyStrong = double(noiseData.postSumEarly > noiseData.thresh);

noiseData.postSlopeEarly = FitCPPSlope(permute(noiseData.respLockedCumul,[1,3,2]), tLims, noiseData.times); % just in early window
noiseData.postSlopeEarlyStrong = double(noiseData.postSlopeEarly > noiseData.thresh);



%% make Strongs have NaN where their inputs did (i.e. missings, or interrupted for postSumStrong)

noiseNames = fieldnames(noiseData);
noiseNames = noiseNames(cellRegexpi(noiseNames, 'Strong')>0); % only do this for the strongs

for i = 1:length(noiseNames)
    noiseData.(noiseNames{i})(isnan(noiseData.(noiseNames{i}(1:end-6)))) = NaN;
end


%% remove flagged

% remove from noise 
noiseNames = fieldnames(noiseData);
noiseNames(ismember(noiseNames, noiseData.ignoreNames)) = []; % not the number of frames or sample times
nNoise = length(noiseNames);

for i = 1:nNoise
    
    noiseData.(noiseNames{i})(repmat(isFlagged',1,1,size(noiseData.(noiseNames{i}),3))) = NaN;
    noiseData.(noiseNames{i})(ppsToExclude,:,:) = NaN;
    
    % also remove bad RTs?
    rtLims = [100 1540];
    noiseData.(noiseNames{i})(repmat(~isBetween(behData.RT, rtLims), 1,1,size(noiseData.(noiseNames{i}),3))) = NaN;
    
end
%% save

save('./Saves/NoiseData.mat','noiseData','firstTilt')