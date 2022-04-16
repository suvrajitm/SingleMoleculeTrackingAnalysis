% fluorobancroft.m
% 
% Find the center of a 2D intensity distribution using the FluoroBancroft
% algorithm -- Andersson, S. Localization of a fluorescent source without 
%    numerical fitting. Opt. Express 16, 18714-18724 (2008)
%
% Inputs
%   I  : 2D intensity distribution (i.e. a grayscale image)
%        Must be square
% 
%   sigma : width of the intensity distribution (px)
%   bkg   : background intensity
%           *WARNING* Extremely inaccurate if realistic background is used;
%           use bkg = 0 !
%
% Outputs
%   xc, yc : the center of radial symmetry,
%            px, from px #1 = left/topmost pixel
%            So a shape centered in the middle of a 2*N+1 x 2*N+1
%            square (e.g. from make2Dgaussian.m with x0=y0=0) will return
%            a center value at x0=y0=N+1.
%            (i.e. same convention as gaussfit2Dnonlin.m and my other 
%            particle finding functions)
%            Note that y increases with increasing row number (i.e. "downward")
%
% Raghuveer Parthasarathy
% Feb. 13, 2012
% last modified Feb. 17, 2012


function [xc yc] = fluorobancroft(I, sigma, bkg)

% Number of grid points
[Ny Nx] = size(I);
if Nx~=Ny
    errordlg('Image must be square! -- cancel')
end
N = Nx;
if bkg>0
    disp('*WARNING* Extremely inaccurate if realistic background is used')
end

Q = [1 0 0 ; 0 1 0]; % for a symmetric distribution

% array of x and y values (faster than repmat):
x_onerow = -(N-1)/2.0:(N-1)/2.0;
x = x_onerow(ones(N, 1), :);
y_onecol = (-(N-1)/2.0:(N-1)/2.0)';
y = y_onecol(:, ones(N, 1));

B = [x(:) y(:) ones(N*N,1)];
Isub = I - bkg*ones(size(I));
Isub(Isub<=0) = 0.001*bkg;  % avoid negative numbers, for log below
P2 = 2*sigma*sigma*log(Isub);
alpha = 0.5*(x(:).*x(:) + y(:).*y(:) + P2(:));

Bdagger = pinv(B);  % Moore-Penrose inverse
temp =  Q*Bdagger*alpha;

xc = temp(1);
yc = temp(2);

% Return output relative to upper left coordinate
xc = xc + (Nx+1)/2.0;
yc = yc + (Ny+1)/2.0;
