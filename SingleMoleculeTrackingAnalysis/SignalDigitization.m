% Script for digitization of object signals

% This script requires the object intensity matrix as generated by object
% detection and tracking codes to be in the workspace


% Author       : Cheemeng Tan
% Affiliation  : Carnegie Mellon University
% Date Created : Jan 14, 2012
% Date Updated : Jan 28, 2012


% Digitization of the object signal


%%%%%% Cheemeng Tan%%%%%%%
%%%
%%% Modified by Suvrajit Maji Jun 17, 2014
%%% Added*: Added to include individual ON frame status for digitization of
%%% signals and calculate ON time statistics and Photon counts

data = AllObjIntensity;

%plotToggle=1: plot all figures for debugging
plotToggle=1;
if plotToggle==1
    togglefig wvf_figure;
end

bleachInt=[];
denoisedata=[];
ObjbleachON = zeros(1,size(data,2));
SigDigitizedAllObjStatus_wvf=zeros(size(SigDigitizedAllObjStatus));
for i=1:size(data,2)
    %denoise data
    [ad] = wden(data(:,i),'sqtwolog','s','sln',8,'haar');
    denoisedata(:,i) = ad;
    %calculate distribution of events
    low=min(ad);
    high=max(ad);
    int=(high-low)/10;
    x=low-int:int:high+int;
    y=histc(ad,x);
    
    %detect peaks and valley in the histogram
    %
    %peakdet(y,alpha,x);
    %Change alpha to define minimum length of events that are
    %considered to be positive events
    %
    [maxP,minP]=peakdet(y,1,x);
    
    %identify length of ON events
    %    bleachInt contains ON frames
    temp=[];
    
    if size(maxP,1)==2
        index=find(x>minP(1,1));
        temp=sum(y(index));
        bleachInt=[bleachInt;temp];
        
        % Added* to include individual ON frame status
        onwvf_ids = find(ad > minP(1,1));
        % Added* to include individual ON frame status
        for l = 1:size(temp,1)
            ObjbleachON(l,i) = temp(l);
        end
    elseif size(maxP,1)==3
        index1=find(x>minP(1,1)& x<minP(2,1));
        index2=find(x>minP(2,1));
        
        % Added* to include individual ON frame status
        onwvf_ids = [ find(ad > minP(1,1) & ad < minP(2,1)) ; find(ad > minP(2,1))];
        
        temp=[sum(y(index1));sum(y(index2))];
        bleachInt=[bleachInt;temp];
        
        % Added* to include individual ON frame status
        for l = 1:size(temp,1)
            ObjbleachON(l,i) = temp(l);
        end
    end
    %plot figures for debugging
    if plotToggle==1
        clf;
        subplot(2,1,1);
        plot(data(:,i),'b.');hold all;
        subplot(2,1,1);
        plot(ad,'r-');
        title(temp);
        subplot(2,1,2);
        hist(ad,x);hold all;
        plot(maxP(:,1),maxP(:,2),'ro');
        %         pause;
    end
    % Added*
    SigDigitizedAllObjStatus_wvf(onwvf_ids,i) = 1;
end

% Added*
AllObjbleachON = sum(ObjbleachON,1);
%plot histogram of ON time
togglefig wvf_ONtime;
hist(AllObjbleachON,0:50:size(data,1));

% Added* to include individual ON frame status and Photon counts
ObjIntensityPeakT = ObjIntensityPeak';
ObjIntensitySumT = ObjIntensitySum';
stuniqueObjIntensityPeakFF = stuniqueObjIntensityPeak(1,:)'; %% First Frame for intensity distribution
stuniqueObjIntensitySumFF = stuniqueObjIntensitySum(1,:)';

use_alt_wvf=0;
if use_alt_wvf > 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Our Code from here %%%%
    %%%% The code below does similar wavelet based denoising as above but tries
    %%%% to threshold using simple method to show the thresholded ON status for
    %%%% each object in each frame . It maybe problematic for low
    %%%% signal-to-noise data
    %%%%
    
    FiltSigDigitizedAllObjStatus_wvf = zeros(size(AllObjIntensity));
    SigDigitizedAllObjStatus_wvf = zeros(size(AllObjIntensity));
    denoisedSignal = zeros(size(AllObjIntensity));
    %
    % input the max noise level for the data here, manually
    NoiseLevel = 13;
    %
    for k = 1:size(AllObjIntensity,2)
        
        % denoisedSignal(:,k) =  cmddenoise(AllObjIntensity(:,k),'haar',5,'s');
        % denoisedSignal(:,k) = wden(AllObjIntensity(:,k),'heursure','s','one',3,'sym8');
        
        % T = 20;                             % choose a threshold of 20
        % denoisedSignal(:,k) = doubledual_S1D(AllObjIntensity(:,k)',T);
        
        % decomposition.
        %lev = 5;
        %[c,l] = wavedec(AllObjIntensity(:,k),lev,'sym8');
        % threshold the decomposition structure [c,l].
        
        % denoisedSignal(:,k) = wden(AllObjIntensity(:,k),'heursure','s','one',lev,'sym8');
        s = AllObjIntensity(:,k);
        [C,L] = wavedec(s,3,'db1');
        [thr,sorh,keepapp] = ddencmp('den','wv',s);
        denoisedSignal(:,k) = wdencmp('gbl',C,L,'db1',1,thr,sorh,keepapp);
        
        dnsig = denoisedSignal(:,k);
        sthresh = 4000;  %% Sum object intensity threshold, change this
        
        msig = mean(AllObjIntensity(:,k)); %% dataset from different files here ,check values
        ssig = std(AllObjIntensity(:,k)); %% dataset from different files here ,check values
        
        thdnsig(k) = (max(dnsig) - min(dnsig))/3 + min(dnsig);
        thsig(k) =   (max(AllObjIntensity(:,k)) - min(AllObjIntensity(:,k)))/3 + min(AllObjIntensity(:,k));
              
        FiltSigDigitizedAllObjStatus_wvf(:,k) = denoisedSignal(:,k) > sthresh ;%min(sthresh,thdnsig(k));
        SigDigitizedAllObjStatus_wvf(:,k) = AllObjIntensity(:,k) > sthresh;%thsig(k);
        
    end  
end    


% Added*
% Show the original, denoised and thresholded signals
for k = 1:size(AllObjIntensity,2)
    
    % can use already calculated mean intensities of the objects instead
    % of reading the image stack again
    
    ObjInt(:,k) = AllObjIntensity(:,k);% can be a box so use 2d mean
    
    togglefig Obj_Trace_denoise
    subplot(4,1,1);
    plot(ObjInt(:,k));
    title(['Orginal Signal: Obj # ',num2str(k)],'Fontsize',fontSize*0.6);
    set(gca,'Fontsize',fontSize*0.6);
    
    subplot(4,1,2);
    stairs(denoisedata(:,k),'r');
    title(['Denoised: Obj # ',num2str(k)], 'Fontsize',fontSize*0.6);
    set(gca,'Fontsize',fontSize*0.6);
    
    subplot(4,1,3);
    stairs(SigDigitizedAllObjStatus_wvf(:,k),'k');
    ylim([0 1.5]);
    title(['Digitized: Obj # ',num2str(k)],'Fontsize',fontSize*0.6);
    set(gca,'Fontsize',fontSize*0.6);
    
    subplot(4,1,4);
    stairs(SigDigitizedAllObjStatus(:,k),'k');
    ylim([0 1.5]);
    title(['Digitized(SPTracking): Obj # ',num2str(k)],'Fontsize',fontSize*0.6);
    set(gca,'Fontsize',fontSize*0.6);
end