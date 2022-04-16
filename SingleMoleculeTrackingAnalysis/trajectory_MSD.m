% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : March 1, 2012
% Date Updated : May 27, 2014

% MSD code : msdanalyzer
% Jean-Yves Tinevez - Institut Pasteur, 2013 <tinevez at pasteur dot fr>

SPACE_UNITS = 'µm';
%SPACE_UNITS = 'nm';
TIME_UNITS = 's';

N_PARTICLES = size(xytrajectory,2);
N_TIME_STEPS = length(img_frames);
N_DIM = 2; % 2D


% copy into the format for msdanalyzer
% try to remove very heterogeneous spots (most likely spots with very long
% trajectories relative to other spots)

traj_len = sort_traj_len;
mean_traj_len = mean(traj_len);
traj_len_std =  1*std(traj_len);

% tracks = cell(N_PARTICLES, 1);
k =1;
tracks = [];
for j = 1:size(xytrajectory,2)
    if size(xytrajectory{j},1) < mean_traj_len + traj_len_std
        tracks{k} = [ xytrajectory{j}(:,1)*dT xytrajectory{j}(:,2)*PixelSize xytrajectory{j}(:,3)*PixelSize];
        k = k + 1;
    end
end

% Typical values taken from studies of proteins diffusing in membranes:
% Diffusion coefficient
D_typical  = 1e-3; % µm^2/s

% Time step between acquisition; fast acquisition!


% Area size, just used to disperse particles in 2D. Has no impact on
% analysis.
SIZE = 2; % µm
lwidth = 0.8;
% time
time = (0 : N_TIME_STEPS-1)' * dT;

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);

togglefig all_trajectory_plots
clf;
ma.plotTracks;
ma.labelPlotTracks;
hl = findobj(gca,'Type','line');
set(hl,'Linewidth',lwidth*0.75,'color','black');
set(gca,'FontSize',fontSize,'Linewidth',lwidth);
print('-dtiff','-zbuffer',printResolution,['all_trajectory_' num2str(maxdisp) '.tif']);


matraj_sel=[];
k=1;
sc=1;
for j = 1:size(ma.tracks,1)
    trj = ma.tracks(j);
    if size(trj{:},1) >9
        matraj_sel{k,1} =  trj{:}/sc;
        k=k+1;
    end
end

for j=1:size(matraj_sel,1)
    togglefig sel_trajectory_plots
    clf;
    masel = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
    masel = masel.addAll(matraj_sel(j));
    masel.plotTracks;
    masel.labelPlotTracks;
    hl = findobj(gca,'Type','line');
    set(hl,'Linewidth',lwidth*0.75,'color','black');
    set(gca,'FontSize',fontSize,'Linewidth',lwidth);
    xv = get(gca,'xlim');
    yv = get(gca,'ylim');
    axis([xv(1)-0.1 xv(1)+0.5 yv(1)-0.1 yv(1)+0.5]);
    axis equal;
    print('-dtiff','-zbuffer',printResolution,['sel_trajectory_' num2str(j) '.tif']);
   
    
end


% MSD analysis
ma = ma.computeMSD;
ma.msd;

togglefig MSD_all
clf;
ma.plotMSD;
set(gca,'FontSize',fontSize,'Linewidth',lwidth);
hl = findobj(gca,'Type','line');
set(hl,'Linewidth',lwidth);
print('-dtiff','-zbuffer',printResolution,['MSD_All_' num2str(maxdisp) '.tif']);

togglefig Ensemble_MSD
clf;
ma.plotMeanMSD(gca, true);
mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
errorbar(t, x, dx, 'k');
[fo, gof] = ma.fitMeanMSD;
plot(fo)
ma.labelPlotMSD;
legend off
set(gca,'FontSize',fontSize,'Linewidth',lwidth);
warning off;
Ensemble_MSD = [t x dx];

%individual D
ma = ma.fitMSD;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

fprintf('Estimation of the diffusion coefficient from linear fit of the MSD curves:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));


togglefig D_dist
clf;
D = ma.lfit.a(good_enough_fit)/2 / ma.n_dim;
hist(D,20);
hl = findobj(gca,'Type','patch');
set(hl,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0 0 0],'Linewidth',lwidth);
xlabel('Diffusion coefficient','FontSize',fontSize);
ylabel('Count','FontSize',fontSize);
set(gca,'FontSize',fontSize,'Linewidth',lwidth);
print('-dtiff','-zbuffer',printResolution,['D_dist_' num2str(maxdisp) '.tif']);
 