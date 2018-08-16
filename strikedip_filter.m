% x, y, z: 3D coordinates of points to fit plane to
% e      : Eccentricity threshold of the two major axes of ellipsoid fit 
%   to the give points, ensuring that enough spread in space is achieved to 
%   get a reliable strike and dip. Ratio of major to minor axes; a larger
%   number indicates a higher tolerance for ellipticity.
% o      : Oblateness threshold of the fit points. Oblateness is the aspect
%   ratio between the second and third axes of the ellipsoid: o = 1-z/y,
%   where y is the second axis and z is third (shortest) axis and y is the
%   second axis, roughly equal to x (ideally). If the ellipsoid is not
%   sufficiently oblate, then the MSE of the points will suffer and a fit
%   plane will not reliably capture the attitude of the bed. Values closer 
%   to zero indicate less oblateness (worse fit) while values closer to one
%   indicate more oblateness (better fit).
% s      : Suppress output?
% 
% This function assumes that points representative of a single strike and
% dip have already been selected, i.e. the spatial coverage (but not 
% geometry) of the points is deemed appropriate for the measurement to be 
% made. The point of this function is to filter out points whose geometries
% are not conducive to structural measurements.
% Requires: strikedip.m

function [strike, dip] = strikedip_filter(x,y,z,e,o,s)
    [coeff,score,latent] = pca([x,y,z]);
    e_points = latent(1)/latent(2);
    if e_points > e
        if ~s
            disp(['Given points are not spread enough to meet ' ...
            'the given eccentricity threshold'])
        end
        strike = NaN; dip = NaN;
        return
    end
    o_points = 1 - latent(3)/latent(2);
    if o_points < o
        if ~s
            disp(['Given points are not oblate enough to meet ' ...
            'the given oblateness threshold'])
        end
        strike = NaN; dip = NaN;
        return
    end
    
    [strike, dip] = strikedip(x,y,z);
    
end