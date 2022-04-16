numUniqueObjects  = max(trajectories(:,end));
countObjTrajFrames = zeros(numUniqueObjects,1);

xoffset = 0.0;
yoffset = 0.0;
frameoffset = 0;

for objID = 1:numUniqueObjects
    objectIdx = trajectories(:,end)==objID;
    countObjTrajFrames(objID) = sum(objectIdx);
    
    xytrajectory(objID).x = trajectories(objectIdx,1) + xoffset;
    xytrajectory(objID).y = trajectories(objectIdx,2) + yoffset;
    xytrajectory(objID).frame =  trajectories(objectIdx,3);
end

% object # to look at the trajectory
[sort_len sort_objID ]= sort(countObjTrajFrames,'descend');
objectNo = sort_objID(1);
roi_size = 20;
objectFrame = xytrajectory(objectNo).frame + frameoffset;
centrex =  round(mean(xytrajectory(objectNo).x));
centrey =  round(mean(xytrajectory(objectNo).y));
xrange_pix =  max(1,centrex - roi_size) : min(centrex + roi_size,img_width);
yrange_pix =  max(1,centrey - roi_size) : min(centrey + roi_size,img_height);

for k = 1:length(objectFrame)

  togglefig ObjectTrajectory
  clf
  %ImageFrame = (sequence(yrange_pix,xrange_pix,objectFrame(k)));
  ImageFrame = (sequence(:,:,objectFrame(k)));
  hold on;
  imshow(ImageFrame,[min(ImageFrame(:)) max(ImageFrame(:))]);
  plot(xytrajectory(objectNo).x(1:k),xytrajectory(objectNo).y(1:k),'y-','LineWidth',2);
  zoom xon
end
impixelinfo;
hold off;
