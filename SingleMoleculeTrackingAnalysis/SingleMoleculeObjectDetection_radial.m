% Script for tracking the single molecule FAPs through
% object detection and calculating the density of objects
% called from : SMObjectKinetics.m/ SMObjectKineticsMultiFile.m
% Scripts used : imgwatershed
% requires Image frame 'Iframe' and the current frame no 'frn'
% to run the script independently
%
% Input :
%       Iframe      : Image frame
%       frn         : frame number
%       sigfact     : number of std. dev for thresholding
%       objpixrange : size of objects in pixels [min max]
%       useLocalThreshold : 1 or 0
%       localSigmaFactor : 2.0
%       blockSize   : [16 16]
%       addpixels   : no of pixels to add ,to extend the box around the
%       object
%       show        : 1 = showing the intermediate figures
%
% Output :
%
%

% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : Dec 8, 2008
% Date Updated : July 28, 2011

function varargout = SingleMoleculeObjectDetection(varargin)


global ObjIdx;
global ObjIdxWbox;
%global useWaterShed;
%global useObjectBox;

if nargin < 2
    error('Not enough Input arguments');
elseif nargin>=2
    % Image frame to be thresholded
    Iframe = varargin{1};
    %the current frame number
    frn = varargin{2};
    sigfactor = 2;
    show = 0;
 
end

if nargin >= 3
    sigfactor = varargin{3};
end

if nargin >= 4
    objpixrange = varargin{4};
    minpix = objpixrange(1);
    maxpix = objpixrange(2);
else
    minpix = 10;
    maxpix = 100;
end

if nargin>=8
    useLocalThreshold = varargin{5};
    localSigmaFactor = varargin{6};
    blockSize = varargin{7};
    addpixels = varargin{8};
    show = varargin{9};
    Ibgrnd = varargin{10};
else
  useLocalThreshold = 0;
  localSigmaFactor = 2.0;
  blockSize = [16 16];  
  addpixels = 3; % add pixels to extend the box 
end



%%% some processing options
useFilter = 0; % not required for clean image
useWaterShed = 0;%% by default don't use watershed
morphoperations =0;% sometimes this is necessary when finding edges
% not required for clean single molecule image

GaussShapeFilter = 1; % use gaussian shape filter for single molecules
NyquistCriteria =  0; % Nyquist Sampling
useObjectBox = 0;


%%% use filtering for noisy image, modify accordingly as needed
if useFilter ==1
    Ih = Iframe;
    H = fspecial('average',8); %%%% for cell images it should be 15/20
    
    Iftmean = imfilter(Iframe,H,'replicate'); % mean filter
    %Iftmed = medfilt2(Iframe,[20 20]); % median filter
    
 
    %se = strel('disk',12);
    %Ih = imtophat(Ih,se);
    
    Ift = Iftmean;
    %Ift =  Iftmed;
    Ih = imsubtract(Iframe,Ift);
   
    % %%% using a bandpass filter instead of the above filtering
    % Ih = bpass(Id,1,10,30);
    
    if show >0
        figure(2);
        imshow(imadjust(Ih));
        title('filtered');
    end
else
    
    Ih = Iframe; %%% no change , no filter used
end

%%%%%% thresholding
mI = mean2(Ih);
sI = std2(Ih);

%%% thresholding factor times std.dev above the mean

%%% try to make it adaptive and more robust, for various example cases
%%%

if useLocalThreshold==1
    Ith = thresholdLocally(double(Ih),blockSize,'SigmaFactor',localSigmaFactor);
else
    th = mI + sigfactor*sI;
    Ith = (Ih >= th);
end



if show > 0
    figure(3);
    imshow(Ith);
    title('thresholded');
end

%Ith = thresholdLocally(double(Ih));

% removing objects with pixel area < 7
Ibw = bwareaopen(Ith,minpix);

% object finding by bwlabel
Lbw  = bwlabel(Ibw,4);
s  = regionprops(Lbw, 'Area'); % logical Lbw to reduce memory

%%%% pixel area will depend on the single molecule dataset being used

%%%%% removing smaller groups of pixels
%%% area > 15-20 pixels for newer objective+lense
idxl = find([s.Area] >= minpix);
idxh = find([s.Area] <= maxpix);
idx_filt = intersect(idxl,idxh);
Ibw = ismember(Lbw,idx_filt);

%%% object label before watershed algorithm
[Lbwobj Nbwobj] = bwlabel(Ibw,4);
Ibwobj = Lbwobj > 0;

if show > 0
    figure(4);
    imshow(Ibwobj);
    title('before watershed');
end

if useWaterShed ==1
    Ibw = Ibwobj;
    
    %%% use some morphorlogical operations if necessary
    %%% modify accordingly as needed
    
    if morphoperations ==1
        % Ibw = img_dilate(Ibw,1); % for cell objects,dilation by 2 reqd
        Ibw = bwmorph(Ibw,'fill');
        Ibw = bwmorph(Ibw,'dilate',1);
        Ibw_nws = Ibw;
        
        if show > 0
            figure;
            imshow(Ibw);
            title('before watershed modified');
        end
        
    end
    
    %%%% applying watershed algorithm to distinguish overlapping objects
    [Lbws0 Nbws0] = imgwatershed(Ibw);
    
    Iws = ~Lbws0;
    if show > 0
        figure;
        imshow(Iws);
        title('after watershed');
    end
    
    %%% used for applying correct labels to the object label matrix after
    %%% applying watershed to correct for any missed object pixel
    Lbwshed = nbrfillLabel(Lbws0);
    
    %%% pixels after watershed
    if show > 0
        figure;
        imshow(Lbwshed);
        title('after watershed corrected label');
    end
    Nbws = Nbws0-1;
    maxlabel = Nbws; % no of  object
    
    Num_obj = Nbws;
    Label_obj = Lbwshed;
    
else
    
    Num_obj = Nbwobj;
    Label_obj = Lbwobj;
    maxlabel = Num_obj;
end

mIobj = zeros(1,maxlabel);
NumShapeFiltObj =  0;
ObjGaussCentroids = zeros(Num_obj,2);

% add # pixels on each of the 4 sides of the original box
addpixels = 2;

ObjIdxFrame = cell(1);
ObjIdxWboxFrame = ObjIdxFrame;
 
for j =1:Num_obj
    
    idxL = find(Label_obj == j);
    
    %%%% find a suitable square box covering the object
    [yy xx] = ind2sub(size(Iframe),idxL);
    
    if addpixels >0
        %%%% defines a m x m matrix containing the object
        [xxb yyb] = findsingleObjectROI(xx,yy,size(Iframe),addpixels);
    else
        xxb = xx; yyb = yy;
    end
    %%%% linear indices
    ObjBox = sub2ind(size(Iframe),yyb,xxb);
    Iobjectwbox = double(Iframe(ObjBox));

    if GaussShapeFilter ==1
        %%% fit the object within the box with 2D gaussian function
        %[param2dgfit corrval Ifit] = fit2DGaussianPSF(Iobjectwbox,Ibgrnd);
        [x0_r(j) y0_r(j)] = radialcenter(Iobjectwbox);
        gaussShapeCorr = 1;%corrval(1);
        %gaussShapeCorr = corrval(1);
        %NyquistVal = sqrt(param2dgfit(4)^2 + param2dgfit(6)^2); % nyquist criteria
        %ObjGaussCentroids{j,1} = param2dgfit(3); %X
        %ObjGaussCentroids{j,2} = param2dgfit(5); %Y
        
        NyquistVal = 2.5; 
        ObjGaussCentroids(j,1) = x0_r(j) + min(xxb(:))-1; %X
        ObjGaussCentroids(j,2) = y0_r(j) + min(yyb(:))-1; %Y
        
        
        %togglefig radialcent
        %clf;
        %imshow(Iobjectwbox,[min(Iobjectwbox(:)) max(Iobjectwbox(:))]);
        %hold on;
        %plot(x0_r(j),y0_r(j),'y+');
    else
        gaussShapeCorr = 1; % ignore the shape correlation
    end
    
    ObjPixelIds = idxL;
    ObjPixelIdsWbox = ObjBox; % a rectangular/square ROI for the Object
    
    if gaussShapeCorr > 0.60
        if NyquistCriteria==1 % shape correlation
            if  NyquistVal > 2.3
                NumShapeFiltObj = NumShapeFiltObj  + 1;
                mIobj(NumShapeFiltObj) = mean2(Iframe(ObjPixelIds));
                ObjIdxFrame{j} = ObjPixelIds;
                ObjIdxWboxFrame{j} = ObjPixelIdsWbox;
                %has the shape filtered objects only if GaussShapeFilter == 1
                %Iobjects(frn,NumShapeFiltObj) = mIobj(NumShapeFiltObj);
            end
        else
            NumShapeFiltObj = NumShapeFiltObj  + 1;
            mIobj(NumShapeFiltObj) = mean2(Iframe(ObjPixelIds));
            ObjIdxFrame{j} = ObjPixelIds;
            ObjIdxWboxFrame{j} = ObjPixelIdsWbox;
            %has the shape filtered objects only if GaussShapeFilter == 1
            %Iobjects(frn,NumShapeFiltObj) = mIobj(NumShapeFiltObj);
        end
        
    else
        NumShapeFiltObj = NumShapeFiltObj  + 1;
        mIobj(NumShapeFiltObj) = mean2(Iframe(ObjPixelIds));
        ObjIdxFrame{j} = ObjPixelIds;
        ObjIdxWboxFrame{j} = ObjPixelIdsWbox;
        %has the shape filtered objects only if GaussShapeFilter == 1
        %Iobjects(frn,NumShapeFiltObj) = mIobj(NumShapeFiltObj);
    end
end

sobj  = regionprops(Label_obj,Iframe,'Area','Eccentricity','Centroid'); %%% centroid of the objects, center of mass 
ObjCentroids = cat(1, sobj.Centroid);
ObjectArea = cat(1, sobj.Area);
ObjEccentricty = cat(1, sobj.Eccentricity);
% if frn==16
 %keyboard;
% end
% reject obejcts which are non-elliptical
if GaussShapeFilter ~=1
    shaperejectids = find(ObjEccentricty > 0.9);
    
    ObjIdxFrame(shaperejectids)=[];
    ObjIdxWboxFrame(shaperejectids) =[];
    ObjCentroids(shaperejectids,:) = [];
    ObjectArea(shaperejectids) = [];
    
    NumShapeFiltObj = size(ObjCentroids,1);
end

varargout{1} = Label_obj;
varargout{2} = NumShapeFiltObj;
varargout{3} = ObjIdxFrame;
varargout{4} = ObjIdxWboxFrame;
varargout{5} = ObjGaussCentroids; %ObjCentroids;
varargout{6} = ObjectArea;
varargout{7} = ObjEccentricty;