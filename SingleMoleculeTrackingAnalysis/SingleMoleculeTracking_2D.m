% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : March 1, 2012
% Date Updated : April 25, 2013

% code for 2d tracking of fluorescence spots 
% base tracking code
% http://glinda.lrsm.upenn.edu/~weeks/idl


Mag = 60;
DetectorPix = 6.5; % 16 micron
PixelSize = DetectorPix/Mag;
dT = 0.8; % s,

positionlist = [xout yout TimePoints];
%positionlist = [xout yout zout TimePoints];
maxdisp = 4; %3, max allowed distance between successive frames

param.mem = 3; %2
param.good = 4;%3
param.dim = 3;
param.quiet = 0;

trajectories = track(positionlist, maxdisp ,param);

% no of the unique trajectories
numUniqueObjects  = max(trajectories(:,end));

% view the trajectory movies and also do the MSD analyis after that
plot_dynamic_trajectory;

% 
% sigma_counts = std2(countObjTrajFrames);
% bin_sz_counts = 3.49*sigma_counts/(length(sigma_counts))^(1/3);
% nbins = ceil((max(countObjTrajFrames)-min(countObjTrajFrames))/bin_sz_counts);

togglefig hist_trajlength;
hist(countObjTrajFrames,20);
xlabel('Frames','FontSize',fontSize);
ylabel('Count','FontSize',fontSize);
set(gca,'FontSize',fontSize,'Linewidth',lwidth);
hl = findobj(gca,'Type','patch');
set(hl,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0 0 0],'Linewidth',lwidth);
print('-dtiff','-zbuffer',printResolution,['hist_trajlen_' num2str(maxdisp) '.tif']);


togglefig hist_trajlength_zoom;
hist(countObjTrajFrames,20);
xlabel('Frames','FontSize',fontSize);
ylabel('Count','FontSize',fontSize);
xlim([0 100]);
ylim([0 250]);
set(gca,'FontSize',fontSize,'Linewidth',lwidth);
hl = findobj(gca,'Type','patch');
set(hl,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0 0 0],'Linewidth',lwidth);
print('-dtiff','-zbuffer',printResolution,['hist_trajlen_zoom' num2str(maxdisp) '.tif']);

