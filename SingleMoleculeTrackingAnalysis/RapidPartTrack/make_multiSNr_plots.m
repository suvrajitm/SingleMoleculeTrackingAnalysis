% make_multiSNr_plots.m
% Tracking test plots for a range of SNr values
% Run after running tracking_tests_RP_Feb2012.m, which calculates tracking accuracies
%
% Input:
%    matfilename : file name of the .mat file containing tracking test
%      results
%    methodstoplot : array containing the "id numbers" of the methods whose
%      results will be plotted (default [1:6]):
%      -- (1) radial symmetry fitting; (radialcenter.m)
%      -- (2) Gaussian fitting via nonlinear least squares (gaussfit2Dnonlin.m)
%      -- (3) Gaussian fitting via maximum likelihood estimation (gaussfit2DMLE.m)
%      -- (4) centroid calculation
%      -- (5) weighted linearized Gaussian fitting gauss2dcirc.m (from Stephen M. Anthony / Granick group)
%      -- (6) FluoroBancroft fitting (fluorobancroft.m)
%    CRBmatfilename : file name of the .mat file with Cramer-Rao lower
%      bound values corresponding to the SNr values examined (variable: CRB)
%      Obtained from the Ward / Ober Lab's FandPLimitTool software
%      *Leave empty* ("[]") to avoid plotting these points
%    printfilenamebase : if not empty, print figures to file with this
%       "base" file name
%
% Raghuveer Parthasarathy
% Feb. 8, 2012
% last modified April 6, 2012

function make_multiSNr_plots(matfilename, CRBmatfilename, methodstoplot, printfilenamebase)

if nargin<4
    printfilenamebase = [];
end
if ~exist('methodstoplot', 'var') || isempty(methodstoplot)
    methodstoplot = 1:6;
end
if methodstoplot(1) ~= 1
    disp('** Warning!  Plot symbols etc. will be "wrong" if radial symmetry tracking isn''t method #1! **')
end
load(matfilename);
if ~isempty(CRBmatfilename)
    load(CRBmatfilename, 'CRB', 'SNrarray_CR');
end

%% 
% Colors for plots
c = zeros(Nmethods,3);
c(1,:) = [1.0 0.4 0.0];  % For radial symmetry fit: red (slightly orange)
c(2,:) = [0.0 0.0 0.0];  % Gaussian fit via NLLS: black
c(3,:) = [0.1 0.6 0.1];  % Gaussian fit via MLE: green
c(4,:) = [0.0 0.3 0.8];  % Centroid fit: blue
c(5,:) = [0.8 0.3 0.8];  % Weighted linearized Gaussian: violet
c(6,:) = [0.1 0.8 0.7];  % FluoroBancroft: cyan

% Markers for plots
mp = repmat('o', [Nmethods 1]);
mp(1) =  'o';  % For radial symmetry fit: circle
mp(2) =  's';  % Gaussian fit via NLLS: square
mp(3) =  'x';  % Gaussian fit via MLE: x
mp(4) =  'd';  % Centroid fit: diamond
mp(5) =  '+';  % Weighted Linearized Gaussian: +
mp(6) =  '*';  % FluoroBancroft: *


%% Make figures

% For axis limits
aSNrmin = 0.8*SNrmin;
aSNrmax = 1.2*SNrmax;

% For legend
Mlegendall = {'Radial Sym.', 'Gauss. NLLS', 'Gauss. MLE', 'Centroid', ...
        'W. Lin. Gauss.', 'FluoroBancroft'};
Mlegend = Mlegendall(methodstoplot);

% Note: I want the radial fit to be first in the legend, but also to be "on
% top," and so I re-plot it in most of the plots below

% Total error vs SNr
htot = figure('name', 'Total error vs SNr', 'Position',[100 50 650 600]);
loglog(SNrarray, toterror(1,:), mp(1), 'color', c(1,:), 'markersize', 11)
hold on
for j=methodstoplot(2:end)
    loglog(SNrarray, toterror(j,:), mp(j), 'color', c(j,:), 'markersize', 11)
end
if ~isempty(CRBmatfilename)
    loglog(SNrarray_CR(1:17), CRB(1:17), '-', 'color', 0.7*[1 1 1], 'linewidth', 2.0);
end
loglog(SNrarray, toterror(1,:), mp(1), 'color', c(1,:), 'markersize', 11) % plot again, to be on top
set(gca,'fontsize', 24)
xlabel('SNr')
ylabel('Total error, pixels')
a = axis;
axis([aSNrmin aSNrmax 3e-3 a(4)])
% The following will change the tick label sizes, not the previously-made axis labels
set(gca,'fontsize', 19)
set(gca,'FontWeight', 'normal')
if ~isempty(CRBmatfilename)
    Mlegend_toterror = [Mlegend, {'CRB'}];
else
    Mlegend_toterror = Mlegend;
end
legend(Mlegend_toterror, 'FontSize', 18)
legend boxoff

% Precision vs SNr
hp = figure('name', 'Precision vs SNr', 'Position',[100 50 650 600]);
loglog(SNrarray, sigma(1,:), 'o', 'color', c(1,:), 'markersize', 11)
hold on
for j=[methodstoplot(2:end) 1]
    loglog(SNrarray, sigma(j,:), mp(j), 'color', c(j,:), 'markersize', 11)
end
if ~isempty(CRBmatfilename)
    loglog(SNrarray_CR(1:17), CRB(1:17), '-', 'color', 0.7*[1 1 1], 'linewidth', 2.0);
end
set(gca,'fontsize', 24)
xlabel('SNr')
ylabel('Precision \sigma , pixels')
a = axis;
axis([aSNrmin aSNrmax 3e-3 a(4)])
% The following will change the tick label sizes, not the previously-made axis labels
set(gca,'fontsize', 19)
set(gca,'FontWeight', 'normal')

% Bias vs SNr
% need to plot and replot one set of symbols to get them to be on top of
% the error bar lines.
hb = figure('name', 'Bias vs SNr', 'Position',[50 50 650 600]);
hold on
for j=methodstoplot
    errorbar(SNrarray, bias(j,:), sigbias(j,:), mp(j), 'color', 0.5*(1+c(j,:)), 'markeredgecolor', c(j,:), 'markersize', 11)
end
errorbarlogx;
set(gca,'fontsize', 24)
xlabel('SNr')
ylabel('Bias')
axis([aSNrmin aSNrmax -0.1 1.1])
% The following will change the tick label sizes, not the previously-made axis labels
set(gca,'fontsize', 19)
set(gca,'FontWeight', 'normal')
set(gca,'box','on')

% Fitting time vs SNr
ht = figure('name', 'Time vs SNr', 'Position',[150 50 650 600]);
loglog(SNrarray, time(1,:)*1000, 'o', 'color', c(1,:), 'markersize', 11)
hold on
for j=methodstoplot(2:end)
   loglog(SNrarray, time(j,:)*1000, mp(j), 'color', c(j,:), 'markersize', 11)
end
set(gca,'fontsize', 24)
xlabel('SNr')
ylabel('Time per fit, ms')
a = axis;
axis([aSNrmin aSNrmax a(3) a(4)])
% The following will change the tick label sizes, not the previously-made axis labels
set(gca,'fontsize', 19)
set(gca,'FontWeight', 'normal')


%% Print figures

if ~isempty(printfilenamebase)
    set(htot, 'PaperPosition', [0.25 2.5 6 6])  % to get a decent aspect ratio 
    print(htot, '-dpng', strcat(printfilenamebase, '_totalerror.png'), '-r300')
    set(hp, 'PaperPosition', [0.25 2.5 6 6])  % to get a decent aspect ratio 
    print(hp, '-dpng', strcat(printfilenamebase, '_precision.png'), '-r300')
    set(hb, 'PaperPosition', [0.25 2.5 6 6])  % to get a decent aspect ratio 
    print(hb, '-dpng', strcat(printfilenamebase, '_bias.png'), '-r300')
    set(ht, 'PaperPosition', [0.25 2.5 6 6])  % to get a decent aspect ratio 
    print(ht, '-dpng', strcat(printfilenamebase, '_time.png'), '-r300')
else
    disp('Not printing graphs.')
end
