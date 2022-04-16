% tracking_tests_RP_Apr2012.m
%
% Characterization of particle fitting algorithms
% Comparing:
% -- (1) radial symmetry fitting: distance to gradient line method
%    (radialcenter.m)
% -- (2) Gaussian fitting via nonlinear least squares (gaussfit2Dnonlin.m)
% -- (3) Gaussian fitting via maximum likelihood estimation (gaussfit2DMLE.m)
% -- (4) centroid calculation
% -- (5) weighted linearized Gaussian fitting gauss2dcirc.m (from Stephen M. Anthony / Granick group)
% -- (6) FluoroBancroft fitting (fluorobancroft.m)
%
% uses
%  fitline.m (RP), psf2d.m, modelimage.m [revised, fast version], in addition to
%  radialcenter.m and other tracking functions (radialcenter.m , 
%  gaussfit2Dnonlin.m, gaussfit2DMLE.m, gauss2dcirc.m (from Granick group),
%  fluorobancroft.m)
%
% Inputs
%    Ntrials : number of tests per signal/noise (default 1000)
%    maxx0 : max. shift of true particle image center in each dimention, px (default 0.5)
%    SNrparams : Signal-to-noise ratios to test; either 
%        -- a single number (default; value = 20)
%        -- a 3-element array containing min(SNr), max(SNr), and the number of
%           SNr to examine.
%        -- an array that is not 3-elements in size that contains the 
%           SNr values to examine
%    outfilename : name of MAT file to output (default
%        'track_test_output.mat')
% Outputs
%    sigma: Precision (px) (error minus linear fit)
%    time : execution time per fit (s)
%    bias : bias (slope of error vs. true position)
%    sigbias  : uncertainty in bias
%    toterror : total error of fit (px)
%
% Raghuveer Parthasarathy
% based on tracking_tests_28Aug2011.m (original version)
%    and ... Feb2012
% Apr. 6, 2012
% last modified Apr. 6, 2012

function [sigma time bias sigbias toterror] = tracking_tests_RP_Apr2012(Ntrials, maxx0, SNrparams, outfilename)

% Initialize random number generator
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

% Parameters for the simulated images
scale = 100; % CCD pixel scale (nm/pixel)
lambda = 530; % wavelength of light (nm)
NA = 1.3;  % numerical aperture
N = 9; % size of the simulated image (NxN pixels).
bk = 10.0;  % background (dark) intensity.
dhr = 1.0;  % high resolution grid size for simulated images, nm
fs = sprintf('Using dhr = %.1f', dhr); disp(fs)

% For fluorobancroft, must input width and background
sigmaFB = 0.61*lambda/sqrt(2)/NA/scale;  % width, px

% Parameters for the tests
if ~exist('Ntrials', 'var') || isempty(Ntrials)
    Ntrials = 1000;  % number of tests
    fs = sprintf('Ntrials = %d',Ntrials); disp(fs)
end
if ~exist('maxx0', 'var') || isempty(maxx0)
    maxx0 = 0.5;  % max. shift of true particle image center, px
    fs = sprintf('maxx0 = %.1f',maxx0); disp(fs)
end
if ~exist('SNrparams', 'var') || isempty(SNrparams)
    SNrparams = 20;  % SNr to examine
    fs = sprintf('Single SNr: %.2f',SNrparams); disp(fs)
end
if ~exist('outfilename', 'var') || isempty(outfilename)
    outfilename = 'track_test_output.mat';  % output file name
end

% Pre-calculate the point spread function (over a "large" area)
maxsize = ((N-1) + 2.0*maxx0)*scale;  % maximum needed size for PSF, nm
npsf = round(maxsize/dhr/2.0);
bigpsf = psf2d(npsf, dhr/1000, NA, 1.3, lambda/1000);  % the point-spread function

% least-squares fitting options for the Gaussian fitting
% get here and send to functions, for speed
lsqoptions = optimset('lsqnonlin');

if length(SNrparams)==1
    SNrarray = SNrparams;  % Just one value
elseif length(SNrparams)==3
    % the range is specified
    SNrmin = SNrparams(1);  
    SNrmax = SNrparams(2);  
    Ntests = SNrparams(3);  % 20 % Number of signal-noise ratios to examine
    SNrarray = logspace(log10(SNrmin), log10(SNrmax), Ntests);  % All signal-to-noise ratios to examine
else
    % the array is specified
    SNrarray = SNrparams;
end

Nmethods = 4;  % number of methods to test

% Allocate memory
sigma  = zeros(Nmethods, length(SNrarray));
time = zeros(Nmethods, length(SNrarray));
bias = zeros(Nmethods, length(SNrarray));
sigbias = zeros(Nmethods, length(SNrarray));
toterror = zeros(Nmethods, length(SNrarray));

for j=1:length(SNrarray)
    fs = sprintf('j= %d of %d; Signal/noise %.3f', j, length(SNrarray), SNrarray(j)); disp(fs)
    
    disp('Creating model images...');
    % ... each with random centers, Poisson noise
    x0 = 2*maxx0*rand(Ntrials,1)-maxx0;
    y0 = 2*maxx0*rand(Ntrials,1)-maxx0;
    %x0 = maxx0*rand(Ntrials,1);  % one quadrant only
    %y0 = maxx0*rand(Ntrials,1);  % one quadrant only
    %y0 = x0;  % for testing "displacements" at 45 degree angle to x,y axes
    z = modelimage(SNrarray(j), N, x0*scale, y0*scale, bk, scale, lambda, NA, dhr, maxx0*scale, bigpsf);
    % save tempz z
    keyboard;
    imshow(z)
    disp('1 Starting radial symmetry fitting')
    x0_r = zeros(1,Ntrials); y0_r = zeros(1,Ntrials);
    tic
    for k=1:Ntrials
        [x0_r(k) y0_r(k)] = radialcenter(z(:,:,k));
    end
    time(1,j) = toc;

    disp('2 Starting Gaussian fitting via least-squares')
    x0_g = zeros(1,Ntrials);     y0_g = zeros(1,Ntrials);
    tic
    for k=1:Ntrials
        [Afit, x0_g(k), y0_g(k), sigmafit, offset] = gaussfit2Dnonlin(z(:,:,k), [], [], [], [], lsqoptions);
    end
    time(2,j) = toc;
    
    disp('3 Starting Gaussian fitting via maximum likelihood')
    tic
    x0_m = zeros(1,Ntrials);
    y0_m = zeros(1,Ntrials);
    for k=1:Ntrials
        [Afit, x0_m(k), y0_m(k), sigmafit, offset] = gaussfit2DMLE(z(:,:,k));
    end
    time(3,j) = toc;
    
    disp('4 Starting centroid fitting')
    tic
    x0_c = zeros(1,Ntrials); y0_c = zeros(1,Ntrials);
    for k=1:Ntrials
        lx = 1:N;
        ly = 1:N;
        sumz = sum(sum(z(:,:,k)));
        x0_c(k) = sum(sum(z(:,:,k)) .* lx) / sumz;
        y0_c(k) = sum(sum(z(:,:,k),2) .* ly') / sumz;
    end
    time(4,j) = toc;
    
    disp('5 Starting weighted linearized Gaussian fitting')
    % background is bk
    tic
    x0_w = zeros(1,Ntrials);
    y0_w = zeros(1,Ntrials);
    for k=1:Ntrials
        lx = 1:N;
        ly = 1:N;
        [x0_w(k),y0_w(k)] = gauss2dcirc(z(:,:,k),repmat(lx,N,1),...
            repmat(ly',1,N),bk);
    end
    time(5,j) = toc;

    disp('6 Starting FluoroBancroft')
    % background is bk
    tic
    x0_fb = zeros(1,Ntrials);
    y0_fb = zeros(1,Ntrials);
    if bk>0
        disp('*WARNING* Extremely inaccurate if realistic background is used.  Using bk=0')
    end
    for k=1:Ntrials
        [x0_fb(k),y0_fb(k)] = fluorobancroft(z(:,:,k), sigmaFB, 0.0);
    end
    time(6,j) = toc;
    
    % Characterize the various methods (needlessly redundant code...)
    [bias(1,j) sigbias(1,j) sigma(1,j) toterror(1,j)] = calcbias(x0, y0, x0_r, y0_r, N);
    [bias(2,j) sigbias(2,j) sigma(2,j) toterror(2,j)] = calcbias(x0, y0, x0_g, y0_g, N);
    [bias(3,j) sigbias(3,j) sigma(3,j) toterror(3,j)] = calcbias(x0, y0, x0_m, y0_m, N);
    [bias(4,j) sigbias(4,j) sigma(4,j) toterror(4,j)] = calcbias(x0, y0, x0_c, y0_c, N);
%     [bias(5,j) sigbias(5,j) sigma(5,j) toterror(5,j)] = calcbias(x0, y0, x0_w, y0_w, N);
    [bias(6,j) sigbias(6,j) sigma(6,j) toterror(6,j)] = calcbias(x0, y0, x0_fb, y0_fb, N);

    % Time per fit
    time(:,j) = time(:,j)/Ntrials;
    
    fs = sprintf('Total error (bias & precision)'); disp(fs)
    fs = sprintf('   Radial:  %.3f px\t Gaussian NLLS %.3f px\t Gaussian MLE %.3f px', ...
        toterror(1,j), toterror(2,j), toterror(3,j)); disp(fs)
    fs = sprintf('    Centroid %.3f px\t Weighted Lin. Gauss %.3f px\t FluoroBancroft %.3f px', ...
        toterror(4,j), toterror(5,j), toterror(6,j)); disp(fs)
end

% Correlations
if length(SNrarray)==1; 
    % Correlation of residuals for radial symmetry and Gaussian methods
    % (To be plotted later)
    % Calculate...
    dx = x0_r-(N+1)/2-x0'; dy = y0_r-(N+1)/2-y0'; err_r = sqrt(dx.*dx + dy.*dy);
    dx = x0_g-(N+1)/2-x0'; dy = y0_r-(N+1)/2-y0'; err_g = sqrt(dx.*dx + dy.*dy);
    dx = x0_m-(N+1)/2-x0'; dy = y0_m-(N+1)/2-y0'; err_m = sqrt(dx.*dx + dy.*dy);
    disp('To make plots, run make_singleSNR_plots.m')
else
    disp('To make plots, run make_multiSNR_plots.m')
end

% Save output 
% Clear all the test images (z), and other things
clear j k s z lx ly

save(outfilename);


%% --------------

    function [bias sigbias sigma toterror] = calcbias(x0, y0, x0_fit, y0_fit, N)
        
        % The error of the fits
        
        % total error in x, y
        % Note that the fitting functions return the center relative to the
        % upper left pixel, so that pixel (N+1)/2 is the image center.
        dx = x0_fit-x0'-(N+1)/2;
        dy = y0_fit-y0'-(N+1)/2;
        
        % Bias
        [Ax, sigA, biasx, sigbiasx] = fitline(x0, dx, false);
        [Ay, sigA, biasy, sigbiasy] = fitline(y0, dy, false);
        bias = -0.5*(biasx+biasy);  % average in x, y
        sigbias = 0.5*(sigbiasx+sigbiasy) + 0.5*abs(sigbiasx-sigbiasy);
        % estimate uncertainty as sum of fitting and uncertainty and
        % different in directions; not perfect
        
        % Precision
        std_x = std(dx - Ax - biasx.*x0');
        std_y = std(dy - Ay - biasy.*y0');
        sigma = sqrt(std_x*std_x + std_y*std_y); % precision of the tracking error
        
        % total error (bias and precision)
        toterror = sqrt(std(dx)*std(dx) +std(dy)*std(dy));
        
    end

end

