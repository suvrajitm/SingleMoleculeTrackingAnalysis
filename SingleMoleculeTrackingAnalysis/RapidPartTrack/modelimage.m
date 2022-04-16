% modelimage.m
% 
% Creates a model of a CCD image of a fluorescent particle
% Create a "high resolution" image, (a point source convolved
% with the point-spread function), pixelate, and incorporate noise.
% Can create a stack of (independent) images
%
% *Does not* create an image with two particles (for testting
% robustness of tracking of the central particle) -- see modelimage_2p.m,
% which uses the slower convolution-based method of generating a simulated
% image.
%
% Calls psf2d.m;
% Uses the poisson-distributed random number generator poissrnd from 
% the MATLAB Statistics toolbox
%
% Inputs
%  SNr  : Signal/noise ratio.  Assumes Poisson noise with mean =
%        sqrt(intensity), so that the simulated CCD image is scaled to have
%        the appropriate peak amplitude . 
%        Default 10.
%  N : size of CCD image (NxN pixels).  Make an odd integer to ensure
%      centering in the middle of a pixel (not required) (default 11)
%  xc, yc : object center, x- and y-coordinate (nm) (default 0.0)
%      If xc and yc are (1D) arrays, create an independent image for each
%      desired center
%      [As in make2Dgaussian.m, y increases with increasing row number
%      (i.e. "downward")]
%  bkg : background (dark) intensity.  (default 10.0)
%  scale : CCD pixel scale (nm/pixel) (default 100 nm/px)
%  lambda : wavelength of light (nm) (default 530 nm)
%  NA : numerical aperture (default 1.3)
%  dhr : grid size of the "high resolution" image (default 2 nm) -- note
%      that this may determine the simulated tracking precision for high
%      Signal/Noise
%  maxx0nm : the maximum possible offset in object center, *nm* -- used for
%      calculating the extent of the psf.  Note that maxx0nm must be larger
%      than any of the xc, yc !
%  bigpsf : (optional) the point spread function, calculated on a grid 
%      with resolution dhr and an extent of N*scale + 2*maxx0nm nm.
%      If empty, calculate by calling psf2d.
%      If possible, calculate beforehand for speed!
% Outputs
%   im : simulated 2D image, or stack of images.  Size NxNxlength(xc).
%  
% Raghuveer Parthasarathy
% begun August 28, 2011
% Significant modification (proper Poisson statistics) February 8, 2012
% Significant modification (simple shifted PSF, delete 2-particle) February 13, 2012
% last modified February 14, 2012


function im  = modelimage(SNr, N, xc, yc, bkg, scale, lambda, NA, dhr, maxx0nm, bigpsf)

% Initialize random number generator (for noise)
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

%% defaults for input parameters
if ~exist('SNr', 'var') || isempty(SNr)
    SNr = 10.0; 
end
if ~exist('N', 'var') || isempty(N)
    N = 11;
end
if ~exist('xc', 'var') || isempty(xc)
    xc = 0.0; %nm
end
if ~exist('yc', 'var') || isempty(yc)
    yc = 0.0; %nm
end
if ~exist('bkg', 'var') || isempty(bkg)
    bkg = 10.0; 
end
if ~exist('lambda', 'var') || isempty(lambda)
    lambda = 530; % nm
end
if ~exist('NA', 'var') || isempty(NA)
    NA = 1.3; 
end
if ~exist('scale', 'var') || isempty(scale)
    scale = 100.0; %nm
end
if ~exist('dhr', 'var') || isempty(dhr)
    % resolution of the "high resolution" image, nm/px
    dhr = 1;
end
if ~exist('maxx0nm', 'var') 
    maxx0nm = 0.5*scale;  % max displacement, *nm*
end

%% Check that max. is large enough
if maxx0nm < max([xc(:); yc(:)])
    size(xc)
    size(yc)
    size(maxx0nm)
    fs = sprintf('Max xc: %.3f, Max yc: %.3f, maxx0nm: %.3f',...
        max(xc(:)), max(yc(:)), maxx0nm); disp(fs);
    errordlg('Error! maxx0nm in modelimage.m must be larger than the x, y center positions!');
end

%%
Nimages = length(xc);  % Should be the same as length(yc); not bothering to check

im = zeros([N N Nimages]);

% To be used for the "high resolution" image
xmax = ((N-1)/2)*scale;  % max x, nm

onecolxhr = -xmax:dhr:xmax;
xhr = repmat(onecolxhr,length(onecolxhr),1);
if ~exist('bigpsf', 'var') || isempty(bigpsf)
    % Calculate the psf over a "large" area, if it wasn't input
    maxsize = (N-1)*scale + 2.0*maxx0nm;  % maximum needed size for PSF, nm
    npsf = round(maxsize/dhr/2.0);
    bigpsf = psf2d(npsf, dhr/1000, NA, 1.3, lambda/1000);  % the point-spread function
end

% For the "CCD image"
x = ((1:N)-((N+1)/2))*scale; % array of x-positions of the CCD pixels
y = ((1:N)-((N+1)/2))*scale; % ... y ...

%% tabulate indices once, to avoid repetition. 
% Could input this, to further improve speed.
%
% Here's a very slow way:
% imidx = repmat(struct('idx', []), [N N]);
% for j=1:N
%     for k=1:N
%         imidx(k,j).idx = find(abs(yhr-y(k))<scale/2.0 & abs(xhr-x(j))<scale/2.0);
%         % positions in the high resolution array that correspond to this CCD pixel
%     end
% end
%
% A faster way (could probably improve further with replicating arrays)
imidx1 = repmat(struct('xidx', [], 'yidx', []), [N 1]);
for k=1:N
    % rows
    imidx1(k).yidx = find(abs(onecolxhr-y(k))<scale/2.0);
    % columns
    imidx1(k).xidx = find(abs(onecolxhr-x(k))<scale/2.0);
end
[nRowshr, nColshr] = size(xhr);
imidx = repmat(struct('idx', []), [N N]);
for j=1:N
    for k=1:N
        allyind = repmat((imidx1(k).yidx)',1,length(imidx1(j).xidx));
        allxind = repmat(imidx1(j).xidx,length(imidx1(k).yidx),1);
        imidx(k,j).idx = sub2ind([nRowshr, nColshr], allyind(:), allxind(:));
        % positions in the high resolution array that correspond to this CCD pixel
    end
end

% For the image scaling
% peak signal intensity, based on desired SNr
S = SNr*SNr;

% Create each model image
for imj = 1:Nimages
    
    %% Create the "high resolution" image
    
    % The convolution of the PSF with a point source is simply the psf,
    % spatially offset.
    % (The previous version of this function, which exists as
    % modelimage_2p.m, convolved a 10nm disk with the psf, as in Cheezum 2001,
    % which was very slow and which gave ~identical output to this simpler 
    % and faster approach.
    % Take the negative of xc, yc so that positive direction is right, up
    imhr = bigpsf(round((-yc(imj)+maxx0nm)/dhr) + 1 : round((-yc(imj)+maxx0nm+(N-1)*scale)/dhr) + 1, ...
                  round((-xc(imj)+maxx0nm)/dhr) + 1 : round((-xc(imj)+maxx0nm+(N-1)*scale)/dhr) + 1);
              
    %% Create the "CCD image"
    imccd = zeros(N);
    for j=1:N
        for k=1:N
%            subim = imhr(abs(yhr-y(k))<scale/2.0 & abs(xhr-x(j))<scale/2.0);
            tempimhrsub = imhr(imidx(k,j).idx);
            imccd(k,j) = mean(tempimhrsub(:));  % mean rather than sum in case the pixelation 
                        % has unequal numbers of high-res bins per CCD bin,
                        % due to rounding issues
        end
    end
    imccd = imccd*S/max(imccd(:)) + bkg;  % rescale, and add const. dark level (background)
    
    %fs = sprintf('Mean total non-background photon count: %.1f', sum(imccd(:))-bkg*numel(imccd));
    %disp(fs);
    
    %% Noise
    % Poisson-distributed noise, so draw intensities from a distribution
    % with mean = "imccd" value
    im(:,:,imj) = poissrnd(imccd);  
end

