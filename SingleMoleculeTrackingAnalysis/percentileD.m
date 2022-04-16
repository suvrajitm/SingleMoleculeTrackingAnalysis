%percentiles calculations

load paper_Dall.mat

% percentiles
ps = [25; 50; 75; 100];

% percD columns are data # 
% rows are percentiles
percD = [];
for k =1 : size(Dall,2)
    D = Dall(:,k);
    D = D(D>0); % only positive D values
    percD = [percD prctile(D,ps)];
end