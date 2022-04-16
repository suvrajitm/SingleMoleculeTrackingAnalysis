% finds background intensity level from regions where there are no objects
 
% Input :
%   I : Image frame
%   n : how many standard deviation above the mean of the frame
% Output :
%  T : mean intensity as threshold
%  S : standard deviation of background intensities
%


% Author       : Suvrajit Maji, CPCB PhD Student 
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : July 17, 2008
% Date Updated : Aug 25, 2010

function [T S]= Bkgdetect(I,n)

if nargin<2  % assuming first input is always the image
    n=3;
end

sI= std2(I);
mI = mean2(I);

% determining background pixels
mlI = I(I <= mI + n*sI); % you have to adjust the threshold according to 
                         % data manually , you may have to change the 
                         % code here

T = mean2(mlI);
S = std2(mlI);