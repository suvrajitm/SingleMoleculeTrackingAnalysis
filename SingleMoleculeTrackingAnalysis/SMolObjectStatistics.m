
% Script to perform single molecuel FAP-Dye ON, OFF state statistics


% Author       : Suvrajit Maji, CPCB PhD Student 
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : Sep 12, 2008
% Date Updated : June 17, 2014


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2);
% Ionp = Ion(Ion>0);
% hist(Ion);
% title('Intensity ''Ion'' histogram');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
mTon = zeros(1,size(Ton,2));
mToff = zeros(1,size(Toff,2));

ObjAllOn = [];

num_images = num_frames;



fprintf('\nObjects with no more than single frame ON times are removed !!')
fprintf('\nGenerating the ON time statistics ');
for j =1:size(Ton,2)
    tonid = find(Ton(:,j) > 0) ;
    ton = Ton(tonid,j)/fps;

    toffid = find(Toff(:,j) > 0) ;
    toff = Toff(toffid,j)/fps;

    if isempty(ton)==0
        mTon(j) = mean(ton);
    end

    if isempty(toff)==0
        mToff(j) = mean(toff);
    else
        fprintf('Object # %d is always ON \n\n',j);
        ObjAllOn = [ObjAllOn j]; %%%% always on , j is the id from Ton/Toff
    end
    percComplete(j,size(Ton,2));
end

%
Tons = Ton(Ton>0)/fps;
Toffs = Toff(Toff>0)/fps;

meanTons=mean(Tons);
medianTons=median(Tons);


%%% for mean of mean
% b1 = min(max(mTon),20);
b1 = 25;
b2 = 25; 

%%% for all events
% bons = min(max(Tons),20);

bons = 25;
boffs = 25;

%xons = [0 (round(max(Tons)/10)+1)*10];
xons = [0 num_images]/fps;
xoffs = [0 num_images]/fps;

xmon =[0 num_images]/fps;
xmoff =[0 num_images]/fps;


fsz = fontSize*0.6;


figure(4);
subplot(1,2,1);
hist(mTon,b1);
title(['<Ton> histogram ','  mean <Ton> : ',num2str(mean(mTon),'%3.1f')],'Fontsize',fsz);
xlabel('Time(s)','FontSize',fsz);ylabel('Counts','FontSize',fsz);
set(gca,'Xlim',xmon,'FontSize',fsz);



subplot(1,2,2);
hist(mToff,b2);
title(['<Toff> histogram ','  mean <Toff> : ',num2str(mean(mToff),'%3.1f')],'FontSize',fsz);
xlabel('Time(Ss)','FontSize',fsz);ylabel('Counts','FontSize',fsz);
set(gca,'Xlim',xmoff,'FontSize',fsz);

[smaxton sonid]= sort(max(Ton),'descend');
[smaxtoff soffid ]= sort(max(Toff),'descend');


togglefig Tons
% subplot(1,2,1);
hist(Tons,bons);
title(['Tons histogram /','  mean Tons: ',num2str(mean(Tons),'%3.1f')],'FontSize',fsz);
xlabel('Time(s)','FontSize',fsz);ylabel('Counts','FontSize',fsz);
set(gca,'Xlim',xons,'FontSize',fsz,'Linewidth',lwidth);
hl = findobj(gca,'Type','patch');
set(hl,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0 0 0],'Linewidth',lwidth);
print('-dtiff','-zbuffer',printResolution,['Tons.tif']);
 

% subplot(1,2,2);
% hist(Toffs,boffs);
% title(['Toffs histogram /','  mean Toffs: ',num2str(mean(Toffs),'%3.1f')],'FontSize',fsz);
% xlabel('Time(S) -->','FontSize',fsz);ylabel('counts -->','FontSize',fsz);
%set(gca,'Xlim',xons,'FontSize',fsz,'Linewidth',lwidth);
%hl = findobj(gca,'Type','patch');
%set(hl,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0 0 0],'Linewidth',lwidth);
%print('-dtiff','-zbuffer',printResolution,['Toffs.tif']);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculating the biochemical reaction rate parameters Kon / Koff

% The ON time is fitted with a single exponential distribution here. If you
% want to fit with a double exponential or different distribution , you can
% write a different script or modify the expfit.m
% 
% Kon = 1/mean(Toffs);
% %Koff = 1/mean(Tons);
% 
% %%% Kd should be given 
% Kd = 0.02; %%% find Kd from other studies
% 
% Koff = Kd/Kon;

for j = 1:size(AllObjStatus,2)
    Obj_Onframes{j} = find(AllObjStatus(:,j)>0); 
    Obj_IntegratedInt(j) = sum(AllObjIntensity(Obj_Onframes{j},j));
    
end


sigma_objInt = std2(Obj_IntegratedInt);
bin_sz_sigma_objInt = 3.49*sigma_objInt /(length(Obj_IntegratedInt)^(1/3));
nbin_sz_sigma_objInt = ceil((max(Obj_IntegratedInt)-min(Obj_IntegratedInt))/bin_sz_sigma_objInt);

%[n_int x_int] = hist(Obj_IntegratedInt,nbin_sz_sigma_objInt);

nbins = 20;
% [ydata x_int] = hist(Obj_IntegratedInt,nbins);
[ydata x_int] = hist(Tons,nbins);
xdata = 1:numel(x_int);
[fitted_curve bgof fit_output] =  expfit(xdata,ydata);

% fontSize = 22;

showfit = 1;

if showfit ==1
    togglefig exp_fit;
    %hist(Obj_IntegratedInt,20);
    plot(x_int,ydata)
    hold on
    yvals = fitted_curve(x_int)*max(x_int)/nbins;
    plot(yvals, 'r');
    
    xpos = 6; ypos = 25;
    text('Position',[xpos,ypos*0.9],'String',['N = ' ,num2str(numel(Obj_IntegratedInt),'%.1f')],...
    'VerticalAlignment','bottom',...
    'HorizontalAlignment','left',...
    'FontSize',fsz*0.9);
    
    
    xlabel('ON Time(s)','FontSize',fontSize)
    ylabel('Counts','FontSize',fontSize)
    title({[' R^2 : ',num2str(bgof.rsquare)]},'FontSize',fontSize);
    set(gca,'FontSize',fontSize);
    hold off
end

showresidual = 1;

single_exp_fit_res = fit_output.residuals;

if showresidual ==1
    togglefig fit_residual
    plot(single_exp_fit_res, 'r');
   
    xlabel('xdata','FontSize',fontSize)
    ylabel('Residual','FontSize',fontSize);
    set(gca,'FontSize',fontSize);
    hold off
end

first_Frame_PeakInt = stuniqueObjIntensityPeak(1,:);

togglefig peakIntensity_firstframe
hist(first_Frame_PeakInt,20);
xlabel('Intensity(a.u)','FontSize',fontSize);
ylabel('Count','FontSize',fontSize);
set(gca,'FontSize',fontSize,'Linewidth',lwidth);
hl = findobj(gca,'Type','patch');
set(hl,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0 0 0],'Linewidth',lwidth);
print('-dtiff','-zbuffer',printResolution,['hist_trajlen_zoom' num2str(maxdisp) '.tif']);


fprintf('\nFPS = %.2f \nmean ON time(s) : %.2f\nmedian ON time(s): %.1f\n',fps,meanTons, medianTons);
