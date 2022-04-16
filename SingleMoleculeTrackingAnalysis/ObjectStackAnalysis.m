% Script for analysis of puffs Intensity profiles in all the 
% ROIs through stack
%  
% Can be executed along with 'PuffAnalysis' or independently

% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : April 29, 2011
% Date Updated : Oct  21, 2011

tic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% parsing all the possible puff events

global  uniqueObjIdx uniqueObjStatus  uniqueObjID stuniqueObjIdx  stuniqueObjCentroid stuniqueObjIntensityMean stuniqueObjIntensityPeak ... 
    stuniqueObjIntensitySum stuniqueObjArea PeakIntensityDissipationRate uEventArea;


if(~exist('D','var'))
    stuniqueObjIntensityMean = zeros(size(uniqueObjIntensity));
    stuniqueObjIntensityPeak = zeros(size(uniqueObjIntensity));
    stuniqueObjIntensitySum = zeros(size(uniqueObjIntensity));
    stuniqueObjArea = zeros(size(uniqueObjIntensity));
    
    NumLikelyEvents = 0;
    
    for k =1:numUniqueObj
        fprintf('\nObject %d\n',k);
        
        % circular ROI setup
        t = 0:pi/20:2*pi;
        
        R0 = 3.0; % radius of the circular ROI
        % centroid
        x0 = uniqueObjCentroids{k}(1) ; 
        y0 = uniqueObjCentroids{k}(2); 
        
        xi = R0*cos(t)+x0;
        yi = R0*sin(t)+y0;
        
        ObjectROImask = poly2mask(xi,yi, size(Iframe,1),size(Iframe,2));
        ObjectROI = find(ObjectROImask);
        bkgsub = 1;
        % after getting the circular ROI mask for the object , 
        % go through the stack
        tic
        objno = k;
        [ObjIntensityMean ObjIntensityPeak ObjIntensitySum ObjArea] = ...
            arrayfun(@(x,y)(findObjectAreaThruStack(x.I,ObjectROI,y,bkgsub)),...
            Istack,Ibgrnd, 'UniformOutput', false);
        toc
        
        objImean = cell2mat(ObjIntensityMean');
        objImax = cell2mat(ObjIntensityPeak');
        objIsum = cell2mat(ObjIntensitySum');
        objArea = cell2mat(ObjArea');
        
        % with the extended ROI 
        stuniqueObjIntensityMean(:,k) = objImean';
        stuniqueObjIntensityPeak(:,k) = objImax';
        stuniqueObjIntensitySum(:,k) = objIsum';
        stuniqueObjIdx{k} = ObjectROI;
        stuniqueObjArea(:,k) = objArea';
        
    end
    
    stuniqueObjCentroid = cell2mat(uniqueObjCentroids');
     
end




