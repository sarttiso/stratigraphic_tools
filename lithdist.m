% This function computes the distance of query points in a stratigraphy to
% the nearest bed of a certain type of lithology. Beds of a minimum
% thickness can be specified as well.
%
% IN:
% beds: vector of thicknesses of all beds in stratigraphy
% lith: Vector of strings or integers that encode the lithologies along the
%   stratigraphy; must be same length as beds
% lithcmp: Lithology against which to compare, must be a member of lith
% samples: vector of heights for which to compute distance to bed of given
%   lithology
% 'min_thickness': Minimum thickness of a bed of the given lithology for
%   which to compute the distance to
% 'sided': (default 'both') Compute minimum distance to top, bottom, or
%   either of a bed of a given lithology; can be 'top','bottom','both'
%
% OUT:
% dist: vector of distances, one distance for each member of samples. a
%   distance of zero indicates that a queried height is actually in a bed
%   of given lithology (if user specified 'both' for 'sided'). NaNs may
%   occur if a user requests a sample height above a bed of a given
%   lihology when they specified 'bottom', or if they request a sample
%   height below a bed of a given lithology when they specified 'top'
%
% TO DO:
% - allow for computing distance from arbitrarily many lithologies
%
%
% Adrian Tasistro-Hart, 30.08.2018

function dist = lithdist(bedthick,lith,lithcmp,samples,varargin)

% parser
parser = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(parser,'beds',@isnumeric)
addRequired(parser,'lith',@(x) isnumeric(x) || iscell(x))
addRequired(parser,'lithcmp',@(x) isnumeric(x) || ischar(x))
addRequired(parser,'samples',@isnumeric)
addParameter(parser,'min_thickness',0,validScalarPosNum)
addParameter(parser,'sided','both',@ischar)

parse(parser,bedthick,lith,lithcmp,samples,varargin{:})

bedthick = parser.Results.beds;
lith = parser.Results.lith;
lithcmp = parser.Results.lithcmp;
samples = parser.Results.samples;
minthick = parser.Results.min_thickness;
side = parser.Results.sided;

% validate beds and lith
assert(length(bedthick) == length(lith),...
    'must specify lithology for each bed')

% make column
bedthick = bedthick(:);

% validate lithcmp
assert(ismember(lithcmp,lith),...
    'requested lithology must be present in vector of lithologies (lith)')

% validate side
side = validatestring(side,{'both','top','bottom'});

% if lithologies given as strings, encode as integers
if iscell(lith)
    [lithunique,~,lithidx] = unique(lith);
    nlith = length(lithunique);
    lithnum = (1:nlith)';  
    % integer encoding of lithologies here
    lith =  lithnum(lithidx);
    % get corresponding integers of desired lithologies for distance
    % computation (i.e. integers for lithcmp)
    lithcmp = find(strcmp(lithcmp,lithunique));
end

% make column
lith = lith(:);

% compute basic quantities
nbeds = length(bedthick);   % number of beds
nsamples = length(samples); % number of samples
height = cumsum(bedthick);  % total stratigraphic height at top of each bed
lithidx = lith == lithcmp;  % indices into beds with lithology of interest

% filter out beds of lithology of interest that are too thin, and then
% encode them as zero
lith(bedthick < minthick & lithidx) = 0;
% update indices into beds with lithology of interest
lithidx = lith == lithcmp;
nbedslith = sum(lithidx); % number of beds with lithology of interest

% compute heights of lower and upper bounds of all beds matching the 
% lithology of interest
lb = height(lithidx)-bedthick(lithidx); % lower bounds
ub = height(lithidx);  % upper bounds


% the minimal distance will always be the smallest distance from the sample
% height to one of the lower or upper bounds of a bed with the lithology of
% interest, unless the sample height is actually in such a bed

% first compute distances of all upper bounds from each sample height
distub = repmat(samples',nbedslith,1) - repmat(ub,1,nsamples);
% if the upper boundary is the closest to a given height, then the distance
% between them will be positive, so we can ignore negative distances
distub(distub < 0) = NaN;

% now compute distances of all lower bounds from each sample height
distlb = repmat(lb,1,nsamples)-repmat(samples',nbedslith,1);
% if the lower boundary is the closest to a given height, then the distance
% between them will be positive, so we can ignore negative distances
distlb(distlb < 0) = NaN;

% return a distance based on what user wants (closest to top, bottom, or 
% either of a bed of certain lithology)
switch side
    case 'both'
        % create distance vector by taking minimum of all distances for each sample
        % height from both distlb and distub
        dist = min([min(distub);min(distlb)])';
        % if a sample is in a bed of the lithology of interest, then it will have a
        % corresponding nan in both distlb and distub, so we need to correct for
        % these points
        inlithidx = logical(sum(isnan(distlb) & isnan(distub)))';
        dist(inlithidx) = 0;
        
    case 'top'
        dist = min(distub)';
        
    case 'bottom'
        dist = min(distlb)';
        
end

end