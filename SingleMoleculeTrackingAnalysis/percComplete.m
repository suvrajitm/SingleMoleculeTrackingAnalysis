function  percComplete(jj,maxjj)
%Reports the percentage of the job done to the screen. 
% 
% PERCOUNT(I,Imax) 
% I is the current iteration between 1 and Imax 
% Imax is the maximum number of iterations 
% 
% Do not print anything to the screen between calls of this function! 
%

% title - s cmspike/perccount vr - 1.2 
% author - bodc/alead date - 2006may16

% Updated 2009May14 
% At suggestion of John D'Errico renamed internal variable "max" to 
% "maxjj" 
% Also following D'Errico's suggestions the following functionality has 
% been added: 
% 1. An invocation check - checks that two input arguments are 
% supplied

% lastCall is the last integer percent displayed 
persistent lastCall; 
% if the number of arguments is 2 
if(nargin == 2) 
% if this is the first run, make lastCall a value that will get the 
% following if loops working from the start. 
if isempty(lastCall) 
lastCall = -1; 
end 
% if the percentage has increased to the next integer, then the 
% displayed percentage will need to be updated 
if(lastCall ~= floor(((jj-1)/maxjj) * 100)) 
% if the progress so far (jj out of maxjj) is not the first 
% representation, then backspace the percentage area to prepare for 
% the next percentage update. 
if(jj ~= 1) 
fprintf(1,'\b\b\b'); 
% if the progress so far IS the first time it is displayed, then 
% display the "context line" of the percentage. 
else 
fprintf(1,'\n\tPercentage complete: '); 
end 
% Calculate the percentage done 
pc_done = num2str(floor(((jj-1)/maxjj) * 100)); 
% if the percentage done so far is in the single digit range 
% (between 0 and 10), then add a zero before it so the percentage 
% displays neatly 
if(length(pc_done) == 1) 
pc_done(2) = pc_done(1); 
pc_done(1) = '0'; 
end 
fprintf(1,'%s%%',pc_done); 
end 

lastCall = floor(((jj-1)/maxjj) * 100); % the last int percent displayed 
% if this is the last part of the total, then display 100% complete and 
% add a few lines of whitespace for neatness 
if(jj == maxjj) 
fprintf(1,'\b\b\b100%%\n\n'); 
end 
else 
error('Error: percComplete needs two input arguments...'); 
end