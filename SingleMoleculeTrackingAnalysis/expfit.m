% Script for exponential/general nonlinear curve fit
% you can change the 'method' to perform different fitting

% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : Feb 19, 2009
% Date Updated : Jan 19, 2010

function [fitted_curve bgof fit_output] =  expfit(X,Y)

% tol =1e-6;
% fopt = optimset('Display','off','TolFun',tol,'LargeScale','off');

x1 = round(0.80*max(X));
y1 = Y(x1);

ymax = mean(Y(x1:end));
ybase = min(Y);
lambda0 = -log(1-(y1-ybase)/ymax)/x1;

method = 'B + A * exp(-lambda * x)';

%Initial-guess here
StartPoint = [ybase ymax lambda0];
Lower = [0 0 0];
Upper = [inf inf inf];

fopt = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',Lower,...
               'Upper',Upper,...
               'Startpoint',StartPoint);
ftype = fittype(method);

[fitted_curve,bgof ,fit_output] = fit(X',Y',ftype,fopt);

end


