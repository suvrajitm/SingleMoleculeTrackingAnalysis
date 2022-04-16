% interleaveplot.m
% 
% simple function to make plots of y1(), y2(), y3(), ... yN() vs the same
% x(), in which the display points of each set are interleaved, so that 
% y(k+1) does not obscure y(k).
% The first points are in the order 1, 2, 3, ..., so 'legend' should work
% properly
%
% inputs
%    h = figure handle
%    x = the "x" array (of length M)
%    y = all the "y" arrays, arranged as [y1(1) y2(1) y3(1) ... yN(1);
%                                         y1(2) y2(2) y3(2) ... yN(2);
%                                         ...
%                                         y1(M) y2(M) y3(M) ... yN(M)]
%    mp = markers; an array of N strings; default all are 'o'
%    c  = colors; an array of N x 3 rgb values; default is from the 'jet'
%    colormap
%
% Raghuveer Parthasarathy
% Feb. 18, 2012

function h = interleaveplot(h, x, y, mp, c)

% Parameter defaults
if ~exist('mp', 'var') || isempty(mp)
    mp = repmat({'o'}, [1 size(y,2)]);  % markers
end
if ~exist('c', 'var') || isempty(c)
    c = jet(size(y,2));  % colors
end

%%

x = x(:);
if length(x) ~= size(y,1)
    errordlg('Error: array sizes incorrect in interleaveplot.m');
end

Nplots = size(y,2);
M = size(y,1);  % number of elements
figure(h)
hold on
% Plot in ever-permuting order, to interleave
for j=1:M
    ytoplot = y(j,:);
    which_to_plot = circshift((1:Nplots),[0 j-1]);
    for k=which_to_plot
        plot(x(j), ytoplot(k), mp(k), 'color', c(k,:));
    end
end
