% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : March 1, 2012
% Date Updated : April 29, 2014

countObjTrajFrames = zeros(numUniqueObjects,1);
xytrajectory = [];
xoffset = 0.0;
yoffset = 0.0;
frameoffset = 0;

traj_jump = zeros(1);
traj_jump_micron=[];
traj_lag = 1;

for objID = 1:numUniqueObjects
    objectIdx = trajectories(:,end)==objID;
    countObjTrajFrames(objID) = sum(objectIdx);
    x = trajectories(objectIdx,1) + xoffset;
    y = trajectories(objectIdx,2) + yoffset;
    frame =  trajectories(objectIdx,3);
    xytrajectory{objID} = [frame x y ];
    
    
    % calculate the jumps 
    % 
    k=1;
    for j=1:length(frame)-1
        df=frame(j+traj_lag)-frame(j);
        if df==traj_lag
            traj_jump(k,1)=sqrt((x(j+traj_lag)-x(j))^2+(y(j+traj_lag)-y(j))^2);
            k=k+1;
        end
    end
        
    traj_jump_micron = [traj_jump_micron ; traj_jump*PixelSize];
end

%fontSize = 10;

% object # to look at the trajectory
[sort_traj_len sort_objID ]= sort(countObjTrajFrames,'descend');
objectNo = sort_objID(1); % 68 ?
roi_size = 40; % typical 20 pixels
objectFrame = xytrajectory{objectNo}(:,1) + frameoffset;
centrex =  round(mean(xytrajectory{objectNo}(:,2)));
centrey =  round(mean(xytrajectory{objectNo}(:,3)));
xrange_pix =  [max(1,centrex - roi_size)  min(centrex + roi_size,img_width)];
yrange_pix =  [max(1,centrey - roi_size)  min(centrey + roi_size,img_height)];

for k = 1:length(objectFrame)
  togglefig ObjectTrajectory
  clf % to clear the frame for each plot otherwise it becomes large size to visualize
  %ImageFrame = (sequence(yrange_pix,xrange_pix,objectFrame(k))); % roi, to
  %fix
  ImageFrame = (sequence(:,:,objectFrame(k)));
  imshow(ImageFrame,[min(ImageFrame(:)) max(ImageFrame(:))]);
  hold on;
  plot(xytrajectory{objectNo}(1:k,2),xytrajectory{objectNo}(1:k,3),'y-','LineWidth',2)
  title(sprintf('Frame(%d:%d, Total : %d) : %d ',objectFrame(1),objectFrame(end),length(objectFrame),objectFrame(k)),'FontSize',fontSize)
  set(gca,'XLim',xrange_pix,'YLim',yrange_pix);
end
%impixelinfo;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n');
% MSD Analysis
trajectory_MSD;