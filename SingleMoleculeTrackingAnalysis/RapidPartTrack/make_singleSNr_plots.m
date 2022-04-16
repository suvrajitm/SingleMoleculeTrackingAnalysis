% make_singleSNr_plots.m
%
% Tracking test plots for a single SNr value
% Run after running tracking_tests_RP_Feb2012.m, which calculates tracking
% accuracies.
%
% Revised version of make_singleSNr_plots.m (the old version of which is 
% renamed make_singleSNr_plots_notinterleaved.m) in which points from the
% various data sets are interleaved to avoid obscuring some sets.
% Calls interleaveplot.m to do interleaving.
% Calls usual_labels_for_tracking_tests.m to make labels, set font size,
% etc.
%
% Rather redundant code; could be made much shorter and more efficient,
% e.g. by moving interleaving to a separate function.
% 
% Input:
%    matfilename : file name of the .mat file containing tracking test
%    results
%    printfilenamebase : if not empty, print figures to file with this
%       "base" file name
%
% Raghuveer Parthasarathy
% Feb. 2012; 
% last modified April 6, 2012

function make_singleSNr_plots(matfilename, printfilenamebase)

if nargin<2
    printfilenamebase = [];
end

load(matfilename);

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

%% 

%% Make figures

h1x = figure;
hold on
disp('Wait: making interleaved plot')
h1x = interleaveplot(h1x, x0, [x0_r' x0_g' x0_m' x0_c'] -(N+1)/2-repmat(x0, [1 4]), mp, c);
usual_labels_for_tracking_tests(h1x, 'True x_0 (px)', 'Error: x_{fit} - x_0 (px)',...
    [-maxx0 maxx0 -maxx0 maxx0], {'Radial Sym.', 'Gaussian NLLS',  'Gaussian MLE','Centroid'})
legend off
axis(maxx0*[-1 1 -0.6 0.8])

h1y = figure;
hold on
disp('Wait: making interleaved plot')
h1y = interleaveplot(h1y, y0, [y0_r' y0_g' y0_m' y0_c'] -(N+1)/2-repmat(y0, [1 4]), mp, c);
usual_labels_for_tracking_tests(h1y, 'True y_0 (px)', 'Error: y_{fit} - y_0 (px)',...
    [-maxx0 maxx0 -maxx0 maxx0], {'Radial Sym.', 'Gaussian NLLS',  'Gaussian MLE','Centroid'})
legend off
axis(maxx0*[-1 1 -0.6 0.8])

% without centroid
h1xb = figure;
hold on
disp('Wait: making interleaved plot')
h1xb = interleaveplot(h1xb, x0, [x0_r' x0_g' x0_m'] -(N+1)/2-repmat(x0, [1 3]), mp, c);
usual_labels_for_tracking_tests(h1xb, 'True x_0 (px)', 'Error: x_{fit} - x_0 (px)',...
    [-maxx0 maxx0 -maxx0 maxx0], {'Radial Sym.', 'Gaussian NLLS',  'Gaussian MLE'})
legend off
axis(maxx0*[-1 1 -0.2 0.2])

h1yb = figure;
hold on
disp('Wait: making interleaved plot')
h1yb = interleaveplot(h1yb, y0, [y0_r' y0_g' y0_m'] -(N+1)/2-repmat(y0, [1 3]), mp, c);
usual_labels_for_tracking_tests(h1yb, 'True y_0 (px)', 'Error: y_{fit} - y_0 (px)',...
    [-maxx0 maxx0 -maxx0 maxx0], {'Radial Sym.', 'Gaussian NLLS',  'Gaussian MLE'})
legend off
axis(maxx0*[-1 1 -0.2 0.2])

% Correlation of residuals for radial symmetry and Gaussian methods
hcorr_g = figure('Position', [50 50 600 600]);
plot(err_g, err_r, 'ko', 'markerfacecolor', 0.8*[1 1 1])
grid on
usual_labels_for_tracking_tests(hcorr_g, 'Gaussian NLLS: total error (px)', ...
    'Radial Sym.: total error (px)', [0 0.1 0 0.1], [])

hcorr_m = figure('Position', [50 50 600 600]);
plot(err_m, err_r, 'ko', 'markerfacecolor', 0.8*[1 1 1])
grid on
usual_labels_for_tracking_tests(hcorr_m, 'Gaussian MLE: total error (px)', ...
    'Radial Sym.: total error (px)', [0 0.1 0 0.1], [])


%% Print figures

if ~isempty(printfilenamebase)
    disp('Wait: writing figure files...')
    set(h1x, 'PaperPosition', [0.25 2.5 6 3.75])  % to get a decent aspect ratio
    print(h1x, '-dpng', strcat(printfilenamebase, '_x.png'), '-r300')
    set(h1y, 'PaperPosition', [0.25 2.5 6 3.75])  % to get a decent aspect ratio
    print(h1y, '-dpng', strcat(printfilenamebase, '_y.png'), '-r300')
    set(h1xb, 'PaperPosition', [0.25 2.5 6 3])  % to get a decent aspect ratio
    print(h1xb, '-dpng', strcat(printfilenamebase, '_xb.png'), '-r300')
    set(h1yb, 'PaperPosition', [0.25 2.5 6 3])  % to get a decent aspect ratio
    print(h1yb, '-dpng', strcat(printfilenamebase, '_yb.png'), '-r300')
    set(hcorr_g, 'PaperPosition', [0.25 2.5 6 6])  % to get a decent aspect ratio
    print(hcorr_g, '-dpng', strcat(printfilenamebase, '_corr_g.png'), '-r300')
    set(hcorr_m, 'PaperPosition', [0.25 2.5 6 6])  % to get a decent aspect ratio
    print(hcorr_m, '-dpng', strcat(printfilenamebase, '_corr_m.png'), '-r300')
else
    disp('not printing')
end