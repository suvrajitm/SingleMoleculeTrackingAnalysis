% Author       : Suvrajit Maji, CPCB PhD Student
% Email        : smaji@andrew.cmu.edu , suvrajit@gmail.com
% Affiliation  : Carnegie Mellon University
% Date Created : March 1, 2012
% Date Updated : April 25, 2013

global NumObjects;
NumObjects(frn) = NumShapeFiltObj;

togglefig Object_Detection;

map = [min(min(Iadj)) max(max(Iadj))];

if ~exist('folder','var')
    folder = 1;
end

subplot(1,2,1);
imshow(Iadj,map);
title({['Image Frame # ',num2str(frn)],['Sample # ' , num2str(folder)]});

subplot(1,2,2);
imshow(Iadj,map);
hold on;
clear Iadj;
if NumObjects(frn)>0
    plot(ObjectCentroids(:,1),ObjectCentroids(:,2),'oc','MarkerSize',6);
    hold off;
end
title({['Image Frame # ',num2str(frn)]; ['# Objects : ',num2str(NumObjects(frn))]});

impixelinfo;
drawnow;
