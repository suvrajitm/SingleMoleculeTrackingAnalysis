useWaveletFilter = 1;
AllObjIntensity = stuniqueObjIntensityMean;

% using the tracking information directly for ON/OFF Status
% calculations
SigDigitizedAllObjStatus = SigTrack;

if useWaveletFilter > 0
    % using mean intensity after tracking, you can use other as you wish
    SignalDigitization;
    %SigDigitizedAllObjStatus = SigDigitizedAllObjStatus_wvf;
end

% Calculate the ON/OFF frames 
SMolObjectOnOff;

% Perform the analysis of the object intensity profile
SMolObjectStatistics;

% Added* 
% Calculate total photon counts for each object before it photobleaches
PhotonBefPB = [];
for j=1:length(AllObjbleachON);
    PhotonBefPB(j) = sum(stuniqueObjIntensitySum(:,j))/camera_gain;
end
PhotonBefPB = PhotonBefPB';

togglefig Photon_beforebleach
hist(PhotonBefPB,20);
title([{'Photon Counts Distribution'},{'(Sum for each object before bleaching)'}] ,'FontSize',fontSize);
xlabel('Photon(a.u)','FontSize',fontSize);ylabel('Count','FontSize',fontSize);
hl = findobj(gca,'Type','patch');
set(hl,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0 0 0],'Linewidth',lwidth);
set(gca,'FontSize',fontSize,'Linewidth',lwidth);
print('-dtiff','-zbuffer',printResolution,['Photonbefbleach_dist.tif']);
 
