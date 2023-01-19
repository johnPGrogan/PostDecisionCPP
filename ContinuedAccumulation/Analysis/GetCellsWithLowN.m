function toRemove = GetCellsWithLowN(thresh, varargin)
% function toRemove = GetCellsWithLowN(thresh, varargin)
% Remove cells with low numbers, while keeping toRemove to [30 1280] shape.
% Will take in any number of [30 1280] matrices, and count the number of
% unique combinations of all their values (cells), to find which trials to
% remove.
% 
% Inputs:
%   thresh: scalar, lowest number of trials within a set of conditions to
%       keep, i.e. n < thresh are flagged toRemove
%   varargin: any number of [30 1280] matrices, each corresponding to a
%       factor you wish to affect the cell counting. e.g. if passing in
%       behData.cond and behData.confInR1, it will remove any trials within
%       each combination of condition & confInR1 where there are <thresh
%       number of those trials.
% 
% Outputs:
%   toRemove = [30 1280] logical matrix, where 1 represents a trial that is
%       in a low-count cell (unique combination of different factors passed
%       in). You can then set those trials to NaN to exclude them
% 
% Requires: CountUnique.m, sq.m, 
% John Grogan, 2022
% 

%% get how many factors used, and unique values of them

nV = length(varargin);
sz = cellfun(@size, varargin, 'UniformOutput',0);
if nV > 1
    if any(diff(cat(1, sz{:})),'all')
        disp(sz);
        error('varargin have different sizes');
    end
end

% get unique values in each
u = cell(1,nV);
for i = 1:nV
    [~,u{i}] = CountUnique(varargin{i});
end
nU = cellfun(@length, u); % count then


scaling = 10.^(0:(nV-1)); % will divide each factor by this, so they remain unique (assuming they are integers)
% check they are integers?
if ~all(cellfun(@(x) all(mod(x,1)==0 | isnan(x),'all'), varargin))
    error('some factors are not integers, so the decimilisation used here may mess things up, please check this');
end

% combine them all, each getting a smaller decimilisation
x = zeros(sz{1}); % preset
for i = 1:nV
    x = x + varargin{i} ./ scaling(i);
end

%% count them
[counts, uniqVals] = CountUnique(x,2);

% check we got the expected number out
if numel(uniqVals) ~= prod(nU) || size(counts,2) ~= prod(nU)
    error('unexpected number of unique combination values returned, please check inputs. n = %d', size(counts,2));
end

counts = reshape(counts, [sz{1}(1),nU(end:-1:1)]); % reshape into [pp factors...], in reverse order as the smallest decimilisation is first
uniqVals = reshape(uniqVals, [1,nU(end:-1:1)]); % also reshape this

lowInds = find(counts<thresh); % find low

% get it into same shape as counts
lowIndsMat = []; 
[lowIndsMat(:,1), lowIndsMat(:,2), lowIndsMat(:,3), lowIndsMat(:,4)] = ind2sub(size(counts), lowInds); % get indices of each


% now for each row in lowIndsMat, set those values in toRemove to 1
toRemove = isnan(x); % preset
for i = 1:size(lowIndsMat,1)
    isLow = x(lowIndsMat(i,1),:) == uniqVals(1, lowIndsMat(i,2), lowIndsMat(i,3), lowIndsMat(i,4));
    toRemove(lowIndsMat(i,1), isLow) = 1;
end


end