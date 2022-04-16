
% Script to produce single molecule Object ON, OFF state blinking graphs


% Author       : Suvrajit Maji, CPCB PhD Student 
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : Sep 12, 2008
% Date Updated : Feb 23, 2009


frames = 1:size(AllObjStatus,1);

%%% objects with ontime of atleast 2 frames 
global objon2fIdx;
global objon2fstatus;

%%%%% sorting on the basis of most longest ON time
% [mxon sortidx] = sort(max(Ton),'descend');
% sObjstatus = objon2fstatus(:,sortidx);
% sObj = objon2fIdx(sortidx);

sObjstatus = objon2fstatus;
sObj = objon2fIdx;
num_totObj = length(sObj);
ObjInt = zeros(length(frames),num_totObj);

ObjectNums = 1:length(sObj);

for k = ObjectNums
    
    % can use already calculated mean intensities of the objects instead
    % of reading the image stack again 
    for frn = frames
        Iob = imread(filedir,'tiff',frn);
        ObjInt(frn,k) = mean2(Iob(sObj{k}));% can be a box so use 2d mean
    end

    figure(1);
    subplot(2,1,1);
    plot(ObjInt(:,k));
    title(['Obj # ',num2str(k)]);

    subplot(2,1,2);
    stairs(sObjstatus(:,k),'r');
    ylim([0 1.5]);
    title(['Obj # ',num2str(k)]);
 
    drawnow;
end





