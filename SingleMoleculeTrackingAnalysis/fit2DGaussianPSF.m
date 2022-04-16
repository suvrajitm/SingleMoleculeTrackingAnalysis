
% fitting a single object intensity profile to a 2-D gaussian %

% Input :
%   Im : Image frame
%   IntBackground : background intensity
%   PSFsigma :
%   objectNo :

% output :
%
%   params : array of fitting parameters; background level,max
%   amplitude, mean x val, mean y val, std x val,std y val
%   Rsq : coefficient of determination
%   Rsq_adj : adjusted Coefficient of determination
%   corrv : correlation value R of the fit
%   Ifit : fitted object intensity distribution

% Author       : Suvrajit Maji, CPCB PhD Student
% Affiliation  : Carnegie Mellon University
% Date Created : Aug 07, 2009
% Date Updated : Sep 10, 2011

function [params Rsq Rsq_adj standard_errors corrv Ifit] = fit2DGaussianPSF(Im,IntBackground,PSFsigma,objectNo)

[sizey sizex] = size(Im);
[X,Y] = meshgrid(1:sizex,1:sizey);
x = X(:);
y = Y(:);
f = Im(:);


if nargin < 1
    error('atleast one input required');
end

if nargin < 2 % assuming that user enters only the image frame
    Ibg = nobjR_thresh(Im);
end

if nargin >=2
    % Ibg = 0;%%%% for background subtracted image
    Ibg = IntBackground;
end

if nargin < 3
    PSFsigma = 1.5;
end

if nargin < 4
    objectNo = [];
end

f_avg = mean(f); %average intensity
diffsq =(f - f_avg).^2; % squared difference of intensity and the mean intensity
data = [x  y];

options = optimset('MaxFunEvals',6e4,'MaxIter',6e3,'TolFun',1e-9,'TolX',1e-6,'Display','off','LargeScale','on');

%Initial-guess here
param0 = [Ibg   max(f)  mean(x(find(f==max(f))))  mean(y(find(f==max(f))))   PSFsigma  PSFsigma ];

% lower and upper bounds here
base_lb = min(f);
base_ub = mean2(Im);

maxI_lb = 0;
maxI_ub = max(f) - base_lb;

mean_x_lb = 0;
mean_x_ub = 512;

mean_y_lb = 0;
mean_y_ub = 512;

sigma_x_lb = 0.4;
sigma_x_ub = 4.0;

sigma_y_lb = 0.4;
sigma_y_ub = 4.0;

lb = [base_lb  maxI_lb  mean_x_lb   mean_y_lb  sigma_x_lb  sigma_y_lb];
ub = [base_ub  maxI_ub  mean_x_ub   mean_y_ub  sigma_x_ub  sigma_y_ub];


[params,resnorm,residual,exitflag,output,lambda,jacobian]=...
    lsqcurvefit(@gauss2dpsf,param0,data,f,lb,ub,options);

[Xs Ys] = meshgrid(1:0.5:sizex,1:0.5:sizey);
xs(:,1) = Xs(:);
xs(:,2) = Ys(:);
Ifits =  gauss2dpsf(params,xs); %fitted gaussian in vector
Ifits = reshape(Ifits,size(Xs));%gaussian reshaped as matrix

Ifit = gauss2dpsf(params,data);

[R P]= corrcoef(Im(:), Ifit(:));
corrv(1) = R(1,2);
corrv(2) = P(1,2);

%[ypred,delta]=nlpredci(@gauss2dfxn,data,par,residual,'jacobian',jacobian);
ci = nlparci(params,residual,'jacobian',jacobian);
ci_x = ci(3,2);
ci_y = ci(4,2);
ci_sigma_x = ci(5,2);
ci_sigma_y = ci(6,2);

% most probable positions of the gaussian peak
x0 = params(3);
y0 = params(4);
%standard deviations along x and y directions
sigma_x = params(5);
sigma_y = params(6);


% 1.96 is .975 quantile of the normal distribution
quantile_v  = 1.96;

%standard error of the mean positions
sigma_avg_x_fit = abs(ci_x - x0)/quantile_v;
sigma_avg_y_fit = abs(ci_y - y0)/quantile_v;

%standard error of the standard deviation of the mean position
sigma_sigma_x_fit = abs(ci_sigma_x - sigma_x)/quantile_v;
sigma_sigma_y_fit = abs(ci_sigma_y - sigma_y)/quantile_v;

standard_errors = [sigma_avg_x_fit sigma_avg_y_fit sigma_sigma_x_fit sigma_sigma_y_fit];

% sum of squared residuals
SSresid = resnorm;

%total sum of squares
SStotal = sum(diffsq);

% R-square coefficient of determination
Rsq = 1-(SSresid/SStotal);

%adjusted R-square
Rsq_adj = 1-((SSresid/SStotal)*(length(f)-1)/(length(f)-length(params)-1));

global showfit;
showfit =0;

if showfit > 0
    togglefig ObjectIntensityProfile;
    h = surf(Im);
    title(['Object # ' ,num2str(objectNo)]);
    set(h, 'AlphaData', 0.1);
    
    togglefig FittedIntensityProfile;
    fith = surf(Xs,Ys,Ifits);
    set(fith, 'AlphaData', 0.5);
    axis([min(x),max(x),min(y),max(y),0.9*min(f),max(f)*1.1])
    xlabel('x [pixel]');
    ylabel('y [pixel]');
    zlabel('intensity [a.u.]');
    title({['R-square = ',num2str(Rsq),', R-square adjusted = ',num2str(Rsq_adj)],['Object # ' ,num2str(objectNo)]})
end

end