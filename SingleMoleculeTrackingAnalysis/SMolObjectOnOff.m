% Script  for calculating the 'ON' and 'OFF' time statistics

% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : July 22, 2008
% Date Updated : June 17, 2014


AllObjStatus = SigDigitizedAllObjStatus;
AllObjIdx = stuniqueObjIdx;

[nf nob] = size(AllObjStatus);


frames = 1:nf;

ontimeL = zeros(nf,nob);
offtimeL = zeros(nf,nob);

N_ontime = zeros(1,nob);
N_offtime = zeros(1,nob);

obj_ontime = zeros(nf/2,nob);%%% half the size of total events is enough
obj_offtime = zeros(nf/2,nob);



%%%%% calculating the 'on' and 'off' states for each object



cObjStatus = AllObjStatus;
%cObjStatus = sObjStatus;


cObjIdx = AllObjIdx;
%cObjIdx = sObj;


fprintf('ON time calculations from the object intensity profile status');
for j =1:size(cObjIdx,2)
    obj_onstat = cObjStatus(:,j);
    [ontimeL(:,j) N_ontime(j)] = bwlabel(obj_onstat,4);
    
    for k1 = 1:N_ontime(j)
        obj_ontime(k1,j) = sum(ontimeL(:,j) == k1);
    end
    
    obj_offstat = ~AllObjStatus(:,j);
    [offtimeL(:,j) N_offtime(j)] = bwlabel(obj_offstat,4);
    
    for k2 = 1:N_offtime(j)
        obj_offtime(k2,j) = sum(offtimeL(:,j) == k2);
    end
    
    percComplete(j,size(cObjIdx,2));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% find objects which have atleast one ontime 'ton' >= **

global objon2fIdx;
global objon2fstatus;

on2fIDs = [];

minONtime = 2;
maxONtime = length(frames);


% Anything which appears after firstNFrames is'nt considered for
% statistics
startFrame = 1;
firstNframes = num_frames; %600
minONtimeObjCount  = 0;
FirstNoff = 0;

% object with events more than 'minTon' frames
TonCutoff = 3;

fprintf('Removing Objects which appear after first %d frames and \nless than %.1f ON time frames for getting statistics',firstNframes,TonCutoff)
if firstNframes == num_frames
    fprintf('All frame chosen, so keeping valid Objects from  all frames ...');
end
for j = 1:size(AllObjIdx,2)
    sizefirstNframes = zeros(1,1);
    if (max(obj_ontime(:,j)) >= minONtime) % minimum frames for being considered to be a valid object
        if(max(obj_ontime(:,j)) <= maxONtime)
            minONtimeObjCount  =  minONtimeObjCount + 1;
            if sum(AllObjStatus(1:firstNframes,j)) >= TonCutoff
                
                firstNL = bwlabel(AllObjStatus(startFrame:startFrame+firstNframes-1,j),4);
                
                for t = 1:max(firstNL)
                    sizefirstNframes(t) = length(find(firstNL==t));            
                end
                
                if  min(sizefirstNframes) >= TonCutoff
                    %sizefirstNframes
                    FirstNoff  =   FirstNoff + 1;
                    on2fIDs  = [on2fIDs j] ;
                end
                
            end
        end
    end
    percComplete(j,size(AllObjIdx,2));
end


%%%%% object with events more than 2 frames
objon2fIdx = AllObjIdx(on2fIDs);
objon2fstatus = AllObjStatus(:,on2fIDs);

num_ons = N_ontime(on2fIDs);

ontime2f = ontimeL(:,on2fIDs);
offtime2f = offtimeL(:,on2fIDs);

%%% Ton & Toff of objects with more than 2 frame events
Ton = sum(obj_ontime(:,on2fIDs),1); % summing on times for each objects
Toff = sum(obj_offtime(:,on2fIDs),1);


%%%%%% calculating the 'ons' cycles;
zont = find(N_ontime < 1);
N_ontimep = N_ontime; %% these are the 'ons' ;
N_ontimep(zont) = N_ontime(zont) + 1 ;


