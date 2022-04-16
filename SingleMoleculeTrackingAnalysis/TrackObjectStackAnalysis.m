% Script for analysis of object Intensity profiles in all the
% ROIs through stack
%
% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : April 29, 2011
% Date Updated : June 06, 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%global  uniqueObjIdx uniqueObjStatus  uniqueObjID stuniqueObjCentroid stuniqueObjIdx stuniqueObjIntensityMean stuniqueObjIntensityPeak ...
%    stuniqueObjIntensitySum stuniqueObjArea PeakIntensityDissipationRate uEventArea;


%if(~exist('ObjectROI','var'))
ObjIntensityMean = zeros(numUniqueObjects);
ObjIntensityPeak = zeros(numUniqueObjects);
ObjIntensitySum = zeros(numUniqueObjects);
ObjArea = zeros(numUniqueObjects);

%for tracking all objects appearing in all the frames through the image stack
fprintf('\n\n\nTotal object tracked : %d ...\n Intensity Profile Analysis Next ...',numUniqueObjects);
for k =1:numUniqueObjects
    % find the trajectory number
    trackobjid = find(trajectories(:,end)==k);
    objframes_k = trajectories(trackobjid,3);
    
    % circular ROI setup
    theta = 0:pi/20:2*pi;
    R0 = 3.0; % radius of the circular ROI
    
    % get the information of each object from the each of the frame
    % in its trajectory
    
    for f=1:size(trackobjid,1)
        tfrn = objframes_k(f);
        
        % centroid
        x0 = trajectories(trackobjid(f),1);
        y0 = trajectories(trackobjid(f),2);
        
        xi = R0*cos(theta)+x0;
        yi = R0*sin(theta)+y0;
        
        %Iframe = imread(fileName,'tiff',frn);
        Iframe = Istack(tfrn).I;
        ObjectROImask = poly2mask(xi,yi, size(Iframe,1),size(Iframe,2));
        ObjectROI = find(ObjectROImask);
        bkgsub = 1;
        
        
        ObjIntensityMean = mean(Iframe(ObjectROI));
        ObjIntensityPeak = max(Iframe(ObjectROI));
        ObjIntensitySum = sum(Iframe(ObjectROI));
        ObjArea = size(ObjectROI,1);
             
        objImean = ObjIntensityMean';
        objImax = ObjIntensityPeak';
        objIsum = ObjIntensitySum';
        objArea = ObjArea';
        
        % with the extended ROI
        stuniqueObjIntensityMean(tfrn,k) = objImean';
        stuniqueObjIntensityPeak(tfrn,k) = objImax';
        stuniqueObjIntensitySum(tfrn,k) = objIsum';
        stuniqueObjIdx{k} = ObjectROI;
        stuniqueObjArea(tfrn,k) = objArea';
        
        % frame where the object first appeared, not necessarily frame 1
        if f==1
            stuniqueObjCentroids_1f(k,:) = [x0 y0];
        end
    end
    
    percComplete(k,numUniqueObjects);
end

% end


% This part is for analysing photobleaching effect of the objects in the
% first frame
%
% for tracking objects appearing in the first frame through the stack
framenum=1;
trackids_f1 = trajectories(:,3)==framenum;
objidx_f1 = trajectories(trackids_f1,4);

% Objects in the trajectories variable are recorded if the trajectory
% lengths are > 2 ( 3 or 4 as set for 'param.good' in the SinglemoleculeTracking_2D.m )
% get the objects in the first frame from the saved list for all frames

ObjectCentroidframe1 = stackObjectCentroids{1};
fprintf('\n First frame Objects Tracked (length >= %d) : %d ...',param.good,size(objidx_f1,1));
fprintf('\n Analysing First frame Objects (Also includes those \n\twith trajectory length < %d): %d ...',param.good,size(ObjectCentroidframe1,1));

rf1=0;
for j = 1:size(ObjectCentroidframe1,1)
    
    x01 = ObjectCentroidframe1(j,1);
    y01 = ObjectCentroidframe1(j,2);
    
    xj = R0*cos(theta)+x01;
    yj = R0*sin(theta)+y01;
    
    ObjectROImask_f1 = poly2mask(xj,yj, size(Iframe,1),size(Iframe,2));
    ObjectROI_f1 = find(ObjectROImask_f1);
    bkgsub = 1;
    % after getting the circular ROI mask for the object ,
    % go through the stack
    [ObjIntensityMean_f1 ObjIntensityPeak_f1 ObjIntensitySum_f1 ObjArea_f1] = ...
        arrayfun(@(x,y)(findObjectAreaThruStack(x.I,ObjectROI_f1,y,bkgsub)),...
        Istack,Ibgrnd, 'UniformOutput', false);
    
    objImean_f1 = cell2mat(ObjIntensityMean_f1');
    objImax_f1 = cell2mat(ObjIntensityPeak_f1');
    objIsum_f1 = cell2mat(ObjIntensitySum_f1');
    objArea_f1 = cell2mat(ObjArea_f1');
    
    % with the extended ROI
    stuniqueObjIntensityMean_f1(:,j) = objImean_f1';
    stuniqueObjIntensityPeak_f1(:,j) = objImax_f1';
    stuniqueObjIntensitySum_f1(:,j) = objIsum_f1';
    stuniqueObjIdx_f1{j} = ObjectROI_f1';
    stuniqueObjArea_f1(:,j) = objArea_f1';
    
    % now replace the profile of the objects which have already been
    % tracked with the corresponding frame intensity values, leave the
    % frames which are not tracked as the roi from first frame
    % 
    ids_f1 = find((stuniqueObjCentroids_1f(:,1) == x01) & (stuniqueObjCentroids_1f(:,2) == y01));
    if ~isempty(ids_f1)
        intf1 = stuniqueObjIntensityMean(:,ids_f1);
        stuniqueObjIntensityMean_f1(intf1 > 0,j) = intf1(intf1 >0) ; 
        rf1 = rf1 + 1;
    end
    percComplete(j,size(ObjectCentroidframe1,1));
end

%end

% threshold the signals for ON/OFF state analysis
% these are from tracked frames , so anything positive is ON frame
SigTrack = double(stuniqueObjIntensityMean > 0);

