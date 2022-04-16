
% Script to produce single molecule Object ON, OFF state blinking graphs


% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : Sep 12, 2008
% Date Updated : Feb 23, 2010


frames = 1:size(cumObjStatus,1);

%%% objects with ontime of atleast 2 frames
global objon2fIdx;
global objon2fstatus;
% 
% 
sObjstatus = objon2fstatus;
sObj = objon2fIdx;
num_totObj = length(sObj);
ObjInt = zeros(length(frames),num_totObj);


startObj = 1160;

for filenum =1:numfiles
    
    filedir = filedirs{filenum};
    
    numObjects = 1%numCumObjects(filenum);
    
    endObj = startObj + numObjects - 1;
    
    % endObj = startObj;
    
    for k = startObj : endObj
        
        % can use already calculated mean intensities of the objects instead
        % of reading the image stack again
        for frn = frames
            Iob = imread(filedir,'tiff',frn);
            ObjInt(frn,k) = mean2(Iob(sObj{k}));% can be a box so use 2d mean
        end
        
        figure(1);
        subplot(2,1,1);
        plot(ObjInt(:,k));
        hold on
        obn = on2fIDs(k);
        plot(denoisedSignal(:,obn),'r');
        hold off;
        xlabel('Frame -->');
        ylabel('Signal Count -->')
        legend('original','denoised');
        title(['Original Object Signal , ',' Obj # ',num2str(k)]);
        
        subplot(2,1,2);
        stairs(sObjstatus(:,k),'r');
        ylim([0 1.5]);
        title(['Digitized Signal ,',' Obj # ',num2str(k)]);
        legend('Digitized');
        drawnow;
    end
    
    startObj = startObj + numObjects;
    
end




