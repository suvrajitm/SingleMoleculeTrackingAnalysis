% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : Oct 14, 2013
% Date Updated : June 05, 2014


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% save all the tracking parameters to mat file
param_list.acquisition.Mag = Mag;
param_list.acquisition.DetectorPix = DetectorPix;
param_list.acquisition.PixelSize = PixelSize;
param_list.acquisition.dT = dT;

param_list.detection.thresh_sig = thresh_sig;
param_list.detection.pix_size_lim = pix_size_lim;

param_list.diffusion.maxdisp = maxdisp;
param_list.diffusion.param = param;

ds = datestr(now,30);
[pathstr, fname, fext] = fileparts(fileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% save the workspace to a .mat file
%trackmatfile = ['Track_' fname '-' ds];
%save(trackmatfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% move the results to separate folders
mkdir('TrackingResults\figures\');
mkdir('TrackingResults\matfiles\');
mkdir('TrackingResults\logfiles\');

[st1,mess1,messid1]=movefile('*.tif', 'TrackingResults\figures\');
[st2,mess2,messid2]=movefile('*.mat', 'TrackingResults\matfiles\');

% rename the log file to include the last time stamp
[st3 mess3 messid3] = movefile(logfile, ['TrackingResults\logfiles\' fname ds '.log'], 'f');

