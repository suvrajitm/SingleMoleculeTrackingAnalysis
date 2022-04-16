% finds a region of interest(square/rectangle) surrounding a object

% Input :
%   xx : x indices 

%   yy : y indices

%   framesize  : image frame size, width and height ..

%   extraPix : no of pixels to added more to the object ROI on all 4 sides
%            (default 2 pixels wider on each side)

% Output :
%   xxb : x indices of the  pixels of a square/rectangular ROI
%         surrounding the specified object

%   yyb : y indices of the pixels of the square/rectangular ROI
%         surrounding the specified object




% Author       : Suvrajit Maji, CPCB PhD Student 
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : July 22, 2009
% Date Updated : Aug 26, 2010

function [xxb yyb] = findsingleObjectROI(xx,yy,framesize,extraPix)

%%%%% this script just serves the purpose of finding a square or a
%%%%% rectangle box surrounding a object, it is not yet robust to all
%%%%% boundary conditions,just need to put few extra conditions,  
%%%%% will implement later !

max_x = max(xx);
min_x = min(xx);
max_y = max(yy);
min_y = min(yy);

if nargin <4
    bdiffxy = 2;
else
    bdiffxy = extraPix;% # pixels to be added more to all 4 sides
end

framesize(1);
framesize(2);

%%%% box length needed to be increased in y-direction

if max_y + bdiffxy <=framesize(1) %%%% within upper limits
    
    yyup = max_y + bdiffxy;
    
else
    %%% else adjust for the boundary condition
%     adjrow = max_y + bdiffxy - framesize;

    yyup = framesize(1);
    
end


if min_y - bdiffxy >= 1 %%%% within lower limits

    yylo = min_y - bdiffxy;
else
    %%% else adjust for the boundary condition
%     adjrow = 1 - min_y + bdiffxy;

    yylo = 1;
end


%%%% box length needed to be increased in x-direction

if max_x + bdiffxy <=framesize(2) %%%% within right limits
    
    xxrt = max_x + bdiffxy;
else
    %%% else adjust for the boundary condition
%     adjcol = max_x + bdiffxy - framesize;

    xxrt = framesize(2);
end



if min_x - bdiffxy >= 1 %%%% within left limits

    xxlt = min_x - bdiffxy;
else
    %%% else adjust for the boundary condition
%       adjcol= 1 - min_x + bdiffxy;

    xxlt = 1;
end


[xxb yyb] = meshgrid(xxlt:xxrt,yylo:yyup);
