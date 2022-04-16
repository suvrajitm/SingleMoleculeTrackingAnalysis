% fitline.m
% routine to fit a line y = A + Bx to data sets [x], [y] 
% also returns uncertainties in A, B (based on fit)
% also returns correlation coefficient (R2)
% 
% Inputs:
%     x : x array
%     y : y array
%     sigy : array of uncertainties in y; can omit (either omit input, or
%            enter empty array)
%     plotopt : IF plotopt==1, plot the points, and the line x array
%               Can omit plotopt argument; then will not plot
%
% considers uncertainties in y in the determination of A, B -- Aug. 6, 2009
% uncertainties (sigy) an optional input.
% For uncertainties in A, B, use deviation from line as uncertainty in y if
% these are not input.
% See discussion in Bevington (1969); uncertainties in A, B from Eq. 6-25
%
% Raghuveer Parthasarathy
% modified Nov. 5, 2008 to calc. correlation coefficient (R2)
% last modified November 23, 2010 (faster; calc. R2 only if desired)

function [A, sigA, B, sigB, R2] = fitline(x, y, sigy, plotopt)

if (length(x(:)) ~= length(y(:)))
    disp('Error!  x, y are not the same length!')
    input('Recommend Control-C to abort. [Enter]');
end
sigyinput = true;  % is the function called with uncertainty values? Default: true

if (nargin < 4)
    plotopt = false;
    if (nargin < 3)
       sigy = ones(size(x));  % irrelevant, uniform uncertainty
       sigyinput = false;
    else
       % There are three input arguments
       if (length(sigy(:))==1)
          % three arguments, and the third is just one number -- so this argument
          % is probably "plotopt", and the function is being called by something
          % written before the consideration of sigy
          plotopt = sigy;
          sigy = ones(size(x));
          sigyinput = false;
       end
    end
end
if isempty(sigy)
    % in case it's empty
    sigy = ones(size(x));
    sigyinput = false;
end

% To ensure arrays are same shape
x = x(:);
y = y(:);
sigy = sigy(:);
N = length(x);

if (N<2)
    % Zero or 1 elements in array
    A = NaN; sigA = NaN; B = NaN; sigB = NaN;
else
    % Least squares linear fit: y = A + Bx
    sx = sum(x./sigy./sigy);
    sxx = sum(x.*x./sigy./sigy);
    sy = sum(y./sigy./sigy);
    sxy = sum(x.*y./sigy./sigy);
    ssig = sum(1./sigy./sigy);
    D = ssig*sxx - (sx*sx);
    A = (sxx*sy - sx*sxy)/D;
    B = (ssig*sxy - sx*sy)/D;
    if nargout==5
        % user wants R2 to be calculated; the slowest step in the function,
        % so do only if requested
        meanx = sum(x)/N;  % faster than mean(x);
        meany = sum(y)/N;  % faster than mean(y);
        R2num1 = sum(x.*y) - N*meanx*meany;
        R2 = R2num1*R2num1 / ...
            ((sum(x.*x) - N*meanx*meanx)*(sum(y.*y) - N*meany*meany));
    else
        R2 = NaN;
    end
    % Uncertainties
    if (N > 2)
        if (sigyinput == false)
            % use deviation from the line as the uncertainty
            devfromline = y - A - B*x;
            sigyfitsq = sum(devfromline.*devfromline)/(N-2);
            sigA = sqrt(sigyfitsq*sxx/D);
            sigB = sqrt(sigyfitsq*N/D);
        else
            sigA = sqrt(sxx/D);
            sigB = sqrt(ssig/D);
        end
    else
        sigA = NaN;
        sigB = NaN;
    end
    
    if plotopt
        figure; plot(x, A + B*x, '-', 'Color', [0.5 0.5 0.5]);
        hold on;
        plot(x, y, 'ko', 'markerfacecolor', [0.2 0.6 1.0]);
    end

end