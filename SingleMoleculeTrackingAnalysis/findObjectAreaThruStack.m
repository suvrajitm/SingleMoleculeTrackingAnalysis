function [IntensityMean IntensityPeak IntensitySum ObjArea] = findObjectAreaThruStack(I,ObjROI,Ibg,backgroundsubtract)

global objno;

if nargin<3
    mI = 0;
    sI = 0;
else
    mI = Ibg.m;
    sI = Ibg.s;
end
Ih = I;

%%%%%% thresholding
% mI = mean2(I); % careful with Ih & I
% sI = std2(I);


%%% thresholding factor times std.dev above the mean
%%% try to make it adaptive, for various example cases
sigfactor = 3;
th = mI-sigfactor*sI;

if backgroundsubtract == 0
    mI = 0;
end
% ObjectPixels = find(Ibw >0);
ObjectPixels = ObjROI;

% check if background subtraction is performed here 
if ~isempty(ObjectPixels) && length(ObjectPixels) >= 4
    ObjArea = length(ObjectPixels);
    Iobject = Ih(ObjectPixels)-mI; % careful with Ih & I 
else
    % center pixel of the box,it should be a background pixel   
    ObjArea = 0;%length(ObjectPixels);
    Iobject = Ih([ObjectPixels;round(length(Ih)/2)])-mI; 
end

IntensityMean = mean(Iobject(:));
IntensityPeak = max(Iobject(:));
IntensitySum = sum(Iobject(:));

end