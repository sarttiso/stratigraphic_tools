% This function computes the uncertainty in a constant sedimentation rate
% given two ages and their locations in a stratigraphy, along with the
% corresponding uncertainties. Assumes Gaussian uncertainties.
%
% IN:
% ages: (2x1 vector) begin and end ages for the sequence of interest
% stdages: (2x1 vector) standard deviations of above ages
% locs: (2x1 vector) begin and end stratigraphic heights of sequence
% stdloc: (2x1 vector) standard deviations of heights
% 'ntrial': (default 10,000) number of times to monte carlo sedimentation 
%   rates
%
% OUT:
% srstd: standard deviation of temporal frequency due to uncertainty in
%   constant sedimentation rate
%
% TO DO:
% 
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 21.09.2018

function [srstd] = sedrate_uncertainty(ages,stdages,locs,stdlocs,varargin)

% parse inputs
parser = inputParser;
addRequired(parser,'ages',@isnumeric)
addRequired(parser,'stdages',@isnumeric)
addRequired(parser,'locs',@isnumeric)
addRequired(parser,'stdlocs',@isnumeric)
addParameter(parser,'ntrial',10000,@isscalar)

parse(parser,ages,stdages,locs,stdlocs,varargin{:});

ages    = parser.Results.ages;
stdages = parser.Results.stdages;
locs    = parser.Results.locs;
stdlocs = parser.Results.stdlocs;
nt      = parser.Results.ntrial;

% validate ages and locations
assert(length(ages)==2,'must provide 2 ages')
assert(length(stdages)==2,'must provide uncertainties for both ages')
assert(length(locs)==2,'must provide 2 locations of ages')
assert(length(stdlocs)==2,...
    'must provide uncertainties for both age locations')

% now generate random sed rates. we need two sets: 
sr = sedrate_gen(ages,stdages,locs,stdlocs,'nt',nt);
     
srstd = std(sr);
        
end