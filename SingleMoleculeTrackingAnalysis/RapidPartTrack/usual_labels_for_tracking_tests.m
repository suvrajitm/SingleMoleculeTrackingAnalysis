% usual_labels_for_tracking_tests.m
% 
% function to avoid redundant code for setting axis labels, etc., for
% tracking tests
% 
% Inputs: 
%   h = figure handle
%   x and y label strings
%   ax : axis limits (leave empty for defaults)
%   array of strings for legend (leave empty for none)
%
% Raghuveer Parthasarathy
% Feb. 18, 2012

function usual_labels_for_tracking_tests(h, xlabelstring, ylabelstring, ax, Mlegend)

figure(h)
set(gca,'fontsize', 22)
xlabel(xlabelstring); ylabel(ylabelstring);
if ~isempty(ax)
    axis(ax)
end
% The following will change the tick label sizes, not the previously-made axis labels
set(gca,'fontsize', 18)
set(gca,'FontWeight', 'normal')
set(gca,'box','on')
if ~isempty(Mlegend)
   legend(Mlegend)
end
legend boxoff
