% This function computes the location and uncertainty in time corresponding 
% to a known stratigraphic height. Two ages are required, and interpolation
% of the time is linear. If the temporal uncertainty of the heights is
% requested, then uncertainties on the input ages and ages locations are
% required.
%
% IN:
% height: vector of heights within sequence of interest to convert to ages
% ages: (2x1 vector) begin and end ages for the sequence of interest
% locs: (2x1 vector) begin and end stratigraphic heights of sequence
% 'ntrial': (default 10,000) number of times to monte carlo sedimentation 
%   rates
% 'stdages': (2x1 vector) standard deviations of above ages
% 'stdlocs': (2x1 vector) standard deviations of heights
% 'stdheight': (default 0) standard deviation of the stratigraphic height.
%   Can be constant for all heights, or vector equal in length to height
%
% OUT:
% age: vectors of ages corresponding to requested heights
%
% TO DO:
% 
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 21.09.2018

function [age,stdage] = height2age(height,ages,locs,varargin)

% parse inputs
parser = inputParser;
addRequired(parser,'height',@isnumeric)
addRequired(parser,'ages',@isnumeric)
addRequired(parser,'locs',@isnumeric)
addParameter(parser,'stdages',[0;0],@isnumeric)
addParameter(parser,'stdlocs',[0;0],@isnumeric)
addParameter(parser,'stdheight',0,@isnumeric)
addParameter(parser,'ntrial',10000,@isscalar)

parse(parser,height,ages,locs,varargin{:});

height = parser.Results.height;
ages    = parser.Results.ages;
stdages = parser.Results.stdages;
locs    = parser.Results.locs;
stdlocs = parser.Results.stdlocs;
stdheight = parser.Results.stdheight;
nt      = parser.Results.ntrial;

% make sure user specified uncertainties
if nargout > 1
    assert(any(stdages ~= 0) || any(stdlocs ~= 0),...
        'user must specify uncertainties in ages and/or age locations')
end

% validate ages and locations
assert(length(ages)==2,'must provide 2 ages')
assert(length(stdages)==2,'must provide uncertainties for both ages')
assert(length(locs)==2,'must provide 2 locations of ages')
assert(length(stdlocs)==2,...
    'must provide uncertainties for both age locations')

% check number of requested heights
nheight = length(height);

% validate stdheight
if length(stdheight) > 1
    assert(length(stdheight) == nheight, ...
        'stdheight must be scalar or equal in length to height')
else
    stdheight = stdheight*ones(nheight,1);
end

% generate many ages with locations
age2 = normrnd(ages(2),stdages(2),nt,1);
age1 = normrnd(ages(1),stdages(1),nt,1);
loc2 = normrnd(locs(2),stdlocs(2),nt,1);
loc1 = normrnd(locs(1),stdlocs(1),nt,1);
m = (age2-age1)./(loc2-loc1);
b = age1 - loc1.*m;

m = repmat(m,1,nheight);
b = repmat(b,1,nheight);

curheight = normrnd(repmat(height,1,nt),repmat(stdheight,1,nt))';
curage = m.*curheight + b;
age = mean(curage);
stdage = std(curage);

end