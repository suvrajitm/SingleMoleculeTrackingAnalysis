% Script for tracking fixed single molecules for calculating the 'ON' and
% 'OFF' time statistics
% This script is called from the script 'SMObjectKinetics.m'

% [Labelframe LabelUniqueObj uniqueObjIdx uniqueObjCentroids uniqueObjID uniqueObjStatus ...
%             uniqueObjIntensity numUniqueObj newObj] = ...
%             SMuniqueStaticObjects(Ifr,Ibgrnd,frn,NumObjects,numUniqueObj,ObjIdx,ObjectCentroids,...
%             uniqueObjIdx,uniqueObjCentroids,uniqueObjID,uniqueObjStatus,uniqueObjIntensity,...
%             Labelframe,LabelUniqueObj,OverlapPercentCriteria,newObj);
%
% Author       : Suvrajit Maji, CPCB PhD Student 
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : July 22, 2008
% Date Updated : April 29, 2011

function [varargout] = SMuniqueStaticObjects(varargin)
 
Iframe = varargin{1};
Ibgrnd = varargin{2};
frn = varargin{3};
NumObjects = varargin{4};
numuobj = varargin{5};
ObjIdx = varargin{6};
ObjectCentroids = varargin{7};
uniqueObjIdx = varargin{8};
uniqueObjCentroids = varargin{9};
uniqueObjID = varargin{10};
uniqueObjStatus = varargin{11};
uniqueObjIntensity = varargin{12};
Labelframe  = varargin{13};
LabelUniqueObj = varargin{14};
OverlapPercentCriteria = varargin{15};
newObj = varargin{16};

%%% try to figure out a way to do it more efficiently , if possible without
%%% using object by object checking , instead directly through manipulation
%%% of matrices, i.e through vectorization instead of for-loops !
%%% A proper tracking code might be more appropriate instead like KALMAN FILTER BASED TRACKING. Current code
%%% might have many superfluous/redundant checks,conditions,constraints !
%%% need to clean up later if it works fine

%%% current method is naive and not robust, it will work fine
%%% for very clean data, well separated uniform sized regular objects,with
%%% no displacements,drifts etc., The pixel overlap percentage will affect
%%% the 'on' and 'off' times for real data
%%% 

% check this script to make sure it is producing expected results as there
% could be discrepencies in detecting the object pixels and determining 
% through pixel overlap if we are seeing the same object in successive
% frames %%%

L_UniqueObj = LabelUniqueObj > 0;

Uid = zeros(NumObjects(frn),3);

%%% for each object in the frame = frn
for j = 1: NumObjects(frn)
   
    if ~isempty(ObjIdx{frn,j}) %%% to ignore blank cells

        %%% finding the intersecting pixels of object j of the current
        %%% frame in the uniqueulative frame (till previous frame)
        objid_cur = L_UniqueObj(ObjIdx{frn,j});
        objid_cur = objid_cur(:);
        
        
        %         curObj = repmat(ObjIdx(frn,j),1,length(uniqueObjIdx));
        %
        %         overlap = cellfun(@ismember,uniqueObjIdx,curObj,'UniformOutput', false);
        %         overlapPercent = cell2mat(cellfun(@(x) sum(x==1)/length(x),overlap,'UniformOutput', false)); % convert to array
        %
        %         uniqueFrameObjId = find(overlapPercent > OverlapPercentCriteria);
        
        %%% if there are intersecting pixels that means object j was
        %%% there in the previous frames , otherwise it is a new object
        %%% which first appears in the current frame
        
        %%% new object are stored in the matrix 'newObj' for
        %%% each frame
        
        ovlp1 = OverlapPercentCriteria*length(objid_cur); % overlap percent of its size

        if sum(objid_cur) >= ovlp1 %%% there are pixels in the area of the 
                               %%% current object in the previous frames
            newObj(frn,j) = 0; %%% so not a new object
                               %%% this object is on in current frame but 
                               %%% we need to know which object is this in
                               %%% the uniqueulative list
                                                                         
                               %%% as the object size may vary with frames
                               %%% the overlap may depend on the object
                               %%% size, so even if an object was 'on' in a 
                               %%% frame , it may not be considered 'on'
                               %%% because the object size may have
                               %%% decreased, but the object size in the 
                               %%% uniqueulative list is still the same and 
                               %%% hence the overlap percent could be small 
                               %%% ,so to take care of this, one more check 
                               %%% should be done for the status of
                               %%% uniqueulative object in the current frame
            
                               
                               %%% find the object in the uniqueulative list
                               %%% corresponding to this object in the
                               %%% current frame and put it's status 'on'
                               
                               
            %%% find which object is this in the previous (or first ) frame 
            %%% if not a new object later ...(to do)
      
            prevobjid = LabelUniqueObj(ObjIdx{frn,j});
            
            prevobjid = prevobjid(prevobjid>0);
            Uid(j,1) = frn;
            Uid(j,2) = j;
            Uid(j,3) = prevobjid(1);
            
        else
            newObj(frn,j) = 1; %%% no overlaping pixels means new object

            %%% object j is a new object so added to the uniqueulative
            %%% list of objects
            numuobj = numuobj + 1;
            uniqueObjIdx{numuobj} = ObjIdx{frn,j};
            uniqueObjCentroids{numuobj} = ObjectCentroids(j,:);
            %num_newobj(frn) = num_newobj(frn) + 1;

            %%% in the uniqueulative object list the new object is the
            %%% latest object hence it's id is the last element of the
            %%% list
       
            newid_in_uniqueobjidx = numuobj;
       
            
            Uid(j,1) = frn;
            Uid(j,2) = j;
            Uid(j,3) = numuobj;
            
            %%% new objects in the current frame are 'on' and they
            %%% are 'off' in the previous frame . jth object is the
            %%% new object in the current frame 'frn' and so jth
            %%% object in all the previous frame should be 'off';
            %%% this jth object corresponds to last object in the uniqueObjIdx

            uniqueObjStatus(frn,newid_in_uniqueobjidx)=1; %%% this will be again 
                                                    %%% checked along with
                                                    %%% all other objects 
                                                    %%% in the list
           % uniqueObjIntensity(frn,newid_in_uniqueobjidx)=mean2(Iframe(ObjIdx{frn,j}));
           LabelUniqueObj(ObjIdx{frn,j}) = numuobj;
                                                    
        end
    end
   
end


L_cur = Labelframe >0 ;%%% Lbwshed or Lbwobj; be careful to use L_cur
                  %%% as if there is object labels like 1,2,3,4,5,etc 
                  %%% instead of 1's in the matrix, then while intersection
                  %%% it would give the label 1,2,3,4 and not 1's , so
                  %%% summing the number of intersecting pixels would not
                  %%% give the actual number of common pixels rather the
                  %%% sum of the object labels , and results will be
                  %%% erroneous ! so use the logical instead

                  
%%% current frame added to unique frame through 'or' operation
% L_UniqueObj = or(L_UniqueObj,L_cur);

num_totObj = numuobj;           %%% unique object list contains all the 
                                %%% objects till the current frame, but the
                                %%% uniqueulative frame contains objects till
                                %%% the previous frame. after the last
                                %%% frame is proccessed, uniqueulative frame
                                %%% will contain all the objects 
                                
%%% This is a double check   
%%% checking the status of all the objects till current frame, in the
%%% current frame . The new objects will automatically have their 
%%% status 'ON' in the current frame


for k = 1:num_totObj
    if ~isempty(uniqueObjIdx{k})
        objidstat = L_cur(uniqueObjIdx{k});%%% be careful with the overlap of
        %%% object indices in the unique list
        %%% with the current frame as the list
        %%% indices might have lower overlap in
        %%% in the cuurent frame,
        %%% whereas the list was built with
        %%% common indices of the objects
        objidstat = objidstat(:);
        
        ovlp2 = OverlapPercentCriteria*length(objidstat); %%% 60 percent overlap in current
        %%% frame is sufficient as it is not crowded
        %%% so, even there is a low overlap ,
        %%% that means it is present in  current
        %%% frame, as size of objects can
        %%% vary with frame
        
        %this check is for position
        if sum(objidstat) >= ovlp2 
            uniqueObjStatus(frn,k) = 1;
            %uniqueObjIntensity(frn,k)=mean2(Iframe(uniqueObjIdx{k}));
        end
    end
    
    % want to have the intensity values for all frames not just the 'ON'
    % ones
    uniqueObjIntensity(frn,k)= mean2(Iframe(uniqueObjIdx{k}));
    
    % this check is for intensity , the object intensity could be low but
    % still not 'off' so make a double check to ensure if it is 'ON'
    if uniqueObjIntensity(frn,k) > Ibgrnd(frn).m + 1.5*Ibgrnd(frn).s
        uniqueObjStatus(frn,k) = 1;
    else
        uniqueObjStatus(frn,k) = 0;
    end
    
end


uniqueObjID = [uniqueObjID ; Uid];

varargout{1} = Labelframe; 
varargout{2} = LabelUniqueObj;
varargout{3} = uniqueObjIdx;
varargout{4} = uniqueObjCentroids;
varargout{5} = uniqueObjID;
varargout{6} = uniqueObjStatus; 
varargout{7} = uniqueObjIntensity; 
varargout{8} = numuobj; 
varargout{9} = newObj;

end
