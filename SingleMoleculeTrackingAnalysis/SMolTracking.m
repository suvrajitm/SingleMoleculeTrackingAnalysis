%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for analysing single molecule images. It performs tracking and
% diffusion MSD analysis and also does intensity profile analysis of
% the single molecules with ON/OFF time statistics etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : Aug 5, 2011
% Date Updated : June 12, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
%
% 1. Localization: The single molecule centroid detection algorithm
% (radial symmetry center) that I have used for localization and tracking,
% (see the paper on Rapid Particle Tracking by R. Parthasarathy 2012)
% is pretty robust and the accuracy is near the theoretical limits of
% Cramer Rao Lower Bound (CRLB), comparable to that of
% Gaussian Maximum Likelihood Estimation (MLE) algorithm.
%
% 2. Error checking: Currently there are not many error checking
% for the different parts of the code.I plan to implement the various
% error checkings as I find them and from user feedbacks. If you
% encounter any problems executing the code, please send me an email
% to my gmail address listed above with the exact error or
% warning message, I will try to respond and if necessary correct
% the problem as soon as possible.
%
% 3. Extra Features: If you would like to do some specific analysis
% or want some particular kind of figure, just let me know,
% if it is feasible, I'll implement it.
%
% 4. We can also perform 3D diffusion analysis, but we have to use 3D
% localization code to get the x,y,z coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
% Put the full file path name under 'the file path' section
% e.g fileName = 'I:\username\imagefile.tif'
% choose proper parameter values for detection, tracking by going into the
% corresponding section and scripts e.g. SingleMoleculeTracking_2D.m
% and trajectory_MSD.m
%
% Output:
% The figures are saved in 'TrackingResults\figures\' and the .mat files
% are saved in 'TrackingResults\matfiles\'.
%
% The .mat files are saved with different time stamps each time the
% code is run. The size of the mat files could be large. So, if you are
% just testing the data or running several times with different
% parameters, be sure to check if you need to keep the files otherwise go
% and delete the .mat files you will not need in order to save space.
%
% If you don't want the whole workspace to be saved to .mat file for some
% reason, just go to SaveTracking.m file and comment out the
% save(trackmatfile). you can always run this command at the end if you
% need to save the workspace.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear workspace
clear global;
clear;

%add the external code directories to the path
addpath 'tracking\';
addpath 'RapidPartTrack\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% File Name : input the file path in this section

fileName = 'I:\Code\SingleMoleculeTracking\7.tif';
% fileName = 'I:\Code\QDtracking\New QD cell data\QD655_MG_405-triple(2)4.tif';
% fileName = 'I:\.tif';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set some global parameters if needed
global filedir img_frames;
global ObjIdx ObjIdxWbox;
global GaussShapeFilter; % use gaussian shape filter for single molecules
filedir = fileName;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Camera Parameters :
% set the camera gain
camera_gain = 300;
fps = 10; % 10 frame per sec ?

% Detection Parameters :  
% Set some initial Parameter values for object detection and thresholding
MaxIntensity = 30000;%%% lower for normal FAP objects,
step = 1; % step for frame reading, 1 = read all frames


% initialize the x,y,z coordinate arrays
xout = [];
yout = [];
zout = [];
TimePoints =[];
numZSlice = 1;
t = 0;

% Background intensity / can be also determined automatically.there is a
% script Bkgdetect.m for each frame
% Bkgdetect uses very basic thresholding from histogram of intensity values
% for the frame
bkgval = 2; % a fixed background intensity value

thresh_sig = 3.0; % 3 sigma used before 
pix_size_lim = [4 35]; % [4 25] used before

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setting up figure printing
fontSize = 22;
printResolution = '-r300';
markersize = 10;
lwidth = 1.5;
paperposition = [0.25 2.5 8.0 6.0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start time for the analysis
tt=tic;

% generate a log file for the data analysis
fprintf('Data File : %s\n',fileName);
[pathstr, fname, fext] = fileparts(fileName);
if exist(fname, 'var')
    delete(fname);
end
logfile = [fname '.log'];

% save the diary for the session to the log file
diary(logfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obtain the image stack information
%for large stacks it causes out of memory problem, use a reduced set
%otherwise don't use a stack variable, read frame by frame.
imginfo = imfinfo(fileName); 
img_frames = 1:numel(imginfo);
img_width = imginfo.Width;
img_height = imginfo.Height;
num_frames = numel(imginfo);

% initialize a stack structure variable 
Istack = struct('I',0);
imgframes = 1:num_frames;%%%
Ibgrnd = struct('m',0,'s',0);

% background value, same for each frame or use Bkgdetect.m to calculate
% frame by frame
Ibgrnd.m = bkgval;

tstart = tic;
Ilastframe = imread(filedir,'tiff',num_frames);

% Creating a stack from input images
fprintf('\nStack Reading ');
sequence = zeros(img_width,img_height,num_frames);
for frnum = 1:length(imgframes)
    frn = imgframes(frnum);
    Iframe = imread(filedir,'tiff',frn);
    
    Ifr = double(Iframe);
    Istack(frn).I = Ifr;
    sequence(:,:,frn) = Ifr;
    percComplete(frnum,length(imgframes));
end

fprintf('\nImage File Information Obtained ...\n');


%Ifr_std = std(Intstack,[],3);  % for blinking particles
%Ifr_std = sequence(:,:,1);
%togglefig stdimage
%imshow(Ifr_std,[min(min(Ifr_std)) max(max(Ifr_std))]);

% start the single molecule particle detection/ thresholding in each frame
fprintf('\nObject Detection ... ');
paramfitall = cell(num_frames);

GaussShapeFilter = 1; % use gaussian shape filter for single molecules

if ~exist('NumShapeFiltObj','var')
    for frnum = 1:step:length(imgframes)
        frn = imgframes(frnum);
        if mod(frn,step)==0
            m = frn/step;
        else
            m = floor(frn/step)+1;
        end
       
        % image reading
        Image = Istack(frn).I;
        I = double(Image);
        
        [Ibgrnd(frn).m Ibgrnd(frn).s] = Bkgdetect(Iframe,2);%%% background
        % check the file nobjR_thresh for proper background level detection
        I(I > MaxIntensity) = Ibgrnd(frn).m ; %%% removing cosmic spots
        
        Imagebs =  I - Ibgrnd(frn).m;
        Iadj = Imagebs;
        
        % object detection
        [Labelframe  NumShapeFiltObj ObjIdxFrame  ObjIdxWboxFrame ObjectCentroids ObjectArea ObjEcc paramfit] = ...
            SingleMoleculeObjectDetection(I,frn,thresh_sig,pix_size_lim,0,3.0,[64 64],0,0,Ibgrnd); % 3 ,[4 30]
        % control data try: 2.5, [15 120]
        
        for obj_id = 1:NumShapeFiltObj
            ObjIdx{frn,obj_id} = ObjIdxFrame{obj_id};
            ObjIdxWbox{frn,obj_id} = ObjIdxWboxFrame{obj_id};
            
            if GaussShapeFilter ==1
                paramfitall{frn,obj_id} = paramfit{obj_id};
                objparam = paramfit{obj_id};
                Ibg{frn,obj_id} = objparam(1);
                I0{frn,obj_id} = objparam(2);
                x0{frn,obj_id} = objparam(3);
                y0{frn,obj_id} = objparam(4);
                sigmax{frn,obj_id} = objparam(5);
                sigmay{frn,obj_id} = objparam(6);
            end
        end

        %for objects on the edge put the x,y coordinates of the centroids
        %to zero/image width/height if the original x/y <0 or x/y > image
        %width/height
        deledgeobj = 1;
        
        % if not needed remove all the edge objects
     
        if deledgeobj > 0
            xct = ObjectCentroids(:,1);
            yct = ObjectCentroids(:,2);
            xnegid = xct<0;
            ynegid = yct<0;
            xwid = xct>img_width;
            yhid = yct>img_height;
              
            xct(xnegid |ynegid | xwid | yhid )=[]; 
            yct(xnegid |ynegid | xwid | yhid )=[];
     
        else
            xct(xct<0)=0;
            yct(yct<0)=0;
            xct(xct>img_width)=img_width;
            yct(yct>img_height)=img_height;
        end
        
        ObjectCentroids = [xct yct];
        stackObjectCentroids{m} = ObjectCentroids;
        
        % storing values for tracking in 2D / 3D
        xout = [xout ; ObjectCentroids(:,1)];
        yout = [yout ; ObjectCentroids(:,2)];
        
        
        %zout_s = [];
        %zout_s = mod(j+1,numZSlice*2)/2 *ones(size(ObjectCentroids,1),1);
        %zout_s(zout_s==0) = numZSlice;
        %zout = [zout; zout_s];
        
        % frame number
        t = t + 1;
        TimePoints = [TimePoints; t*ones(size(ObjectCentroids,1),1)];
        
        % dispaly thersholded object
        displayThresholdedObjects;
        m = m+1;% m = m + 1 is not doing anything here .
        
        percComplete(frnum,ceil(length(imgframes)/step));
    end
    
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform tracking in 2D
fprintf('\nObject Tracking in 2D ... \n');
SingleMoleculeTracking_2D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Intensity profile analysis of the tracked objects through the
% image stack and also perform ON state statistics
TrackObjectStackAnalysis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysing the intensity profile of the objects for extracting meaningful
% information
ObjectSignalAnalysis;

% turn off the diary after finishing the analysis
diary off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the parameter values in .mat file
% fprintf('\nSaving the parameter values to file');
fprintf('\nSaving workspace and parameter values to .mat files ... \n');
SaveTrackingParams;

toc(tt)