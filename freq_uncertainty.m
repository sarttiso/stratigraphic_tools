% This function computes the uncertainty in a temporal frequency given two
% ages that constrain a portion of a stratigraphy and the uncertainties in
% both the ages and their locations in the stratigraphy. A constant or
% roughly constant sedimentation rate is assumed.
%
% IN:
% f: temporal frequencies of interest
% ages: (2x1 vector) begin and end ages for the sequence of interest
% stdages: (2x1 vector) standard deviations of above ages
% locs: (2x1 vector) begin and end stratigraphic heights of sequence
% stdloc: (2x1 vector) standard deviations of heights
% 'ntrial': (default 10,000) number of times to monte carlo sedimentation 
%   rates
% 'stdfreq': uncertainties on frequencies; either a constant which is
%   applied to all frequencies, or a vector of equal length to f
%
% OUT:
% fstd: standard deviation of temporal frequency due to uncertainty in
%   constant sedimentation rate
%
% TO DO:
% - add option of including uncertainty in frequency itself
% 
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 16.08.2018

function [fstd] = freq_uncertainty(f,ages,stdages,locs,stdlocs,varargin)

% parse inputs
parser = inputParser;
addRequired(parser,'f',@isnumeric)
addRequired(parser,'ages',@isnumeric)
addRequired(parser,'stdages',@isnumeric)
addRequired(parser,'locs',@isnumeric)
addRequired(parser,'stdlocs',@isnumeric)
addParameter(parser,'ntrial',10000,@isscalar)
addParameter(parser,'stdfreq',0,@isnumeric)

parse(parser,f,ages,stdages,locs,stdlocs,varargin{:});

f       = parser.Results.f;
ages    = parser.Results.ages;
stdages = parser.Results.stdages;
locs    = parser.Results.locs;
stdlocs = parser.Results.stdlocs;
nt      = parser.Results.ntrial;
stdfreq = parser.Results.stdfreq;

% make column
f = f(:);
nf = length(f);

% validate ages and locations
assert(length(ages)==2,'must provide 2 ages')
assert(length(stdages)==2,'must provide uncertainties for both ages')
assert(length(locs)==2,'must provide 2 locations of ages')
assert(length(stdlocs)==2,...
    'must provide uncertainties for both age locations')

% validate stdfreq
assert(length(stdfreq) == 1 || length(stdfreq) == nf,...
    'stdfreq must be constant or of equal length to f')

% if constant, just make into vector
if length(stdfreq) == 1
    stdfreq = stdfreq*ones(nf,1);
end

% now generate random sed rates. we need two sets: 
% 1) one for converting from time to space. This set reflects the fact 
% that, in our model, an unknown but constant sedimentation rate roughly 
% constrained by our age data takes a temporal frequency of interest and
% maps it to a spatial frequency.
% 2) one from converting from space to time. This set reflects the fact
% that the spatial frequency from above is then converted back to a
% temporal frequency using another estimate of the sedimentation rate,
% which, unless we are very lucky, will not be the same as the true
% sedimentation rate that original encoded the temporal frequency as a
% spatial one.

% these sed rates are for converting the frequency of interest to a spatial
% frequency
sr2space = sedrate_gen(ages,stdages,locs,stdlocs,'nt',nt);
     
% these sed rates are for converting the spatial frequency back into a
% temporal frequency
sr2time = (normrnd(locs(2),stdlocs(2),nt,1) - ...
            normrnd(locs(1),stdlocs(1),nt,1)) ./ ...
           (normrnd(ages(2),stdages(2),nt,1) - ...
            normrnd(ages(1),stdages(1),nt,1));
     
% now we convert from time to space back to time
fobs = zeros(nt,nf);
for ii = 1:nf
    % spatial frequency
    fx = normrnd(f(ii),stdfreq(ii),nt,1)./sr2space;   
    % recovered temporal frequency
    fobs(:,ii) = fx.*sr2time; 
end
 
% the resulting distribution of frequencies is Gaussian in this model, so
% the standard deviation is a useful statistic
fstd = std(fobs);
fstd = fstd(:);

end