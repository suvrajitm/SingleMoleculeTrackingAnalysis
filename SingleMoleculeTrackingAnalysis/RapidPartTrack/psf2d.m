% psf2d.m
% function to calculate the PSF -- returns a 2D matrix.
% psf2d returns a square array: central point +/- np pixels
%
% scale is *microns / pixel* (changed from earlier version)
% NA is the numerical aperture
% nw is the index of refraction of the medium
% lambda is the wavelength, microns
%
% PSF is *normalized* to sum==1
%
% Raghuveer Parthasarathy 9 Feb 2004
% Revised Aug. 28, 2011 (more efficient)


function [psf] = psf2d(np, scale, NA, nw, lambda)

if (scale > 1.0)
    % unlikely, so warn user that he/she may be using pixels/micron rather
    % than microns/pixel
    disp('** ');
    disp('WARNING:  scale seems high -- Are you sure the units are correct?')
    pause(1)
end

c = repmat(-np:np,2*np+1,1);  % Array of column numbers, relative to center
   % e.g.    -1 0 1
   %         -1 0 1
   %         -1 0 1    for np=1
r = repmat((-np:np)',1,2*np+1);  % Array of row numbers, relative to center
   % e.g.    -1 -1 -1
   %          0  0  1
   %          1  1  1  for np=1

d = sqrt(r.*r + c.*c);  % matrix of distances from the center, in pixels

rpos = d*scale;  % matrix of distances, microns


% suppress the divide-by-zero error
warning off MATLAB:divideByZero;

v = (2*pi/lambda)*(NA/nw)*rpos;  % array of 'v' values
besv = besselj(1,v)./v;
psf = 4.0*(besv.*conj(besv));
psf(np+1,np+1) = 1.0;  % correct central value
psf = psf / sum(psf(:)); % normalize

% Restore the warning state.
warning on MATLAB:divideByZero;
