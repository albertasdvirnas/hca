% We test the method on 
% /media/albyback/My Passport/DATA/chromosomes/Data For Albertas/Files/2330P18/
% which we know should come from chromosome 22, roughly around 79100
% position.

% correctly placed using PCC
sum(cellfun(@(x) abs(x.pos(1)-79000) < 200,comparisonStructAll{1}))

% correctly placed using DTW with some Sakoe-Chiba band
sum(cellfun(@(x) abs(x.ucr.pos-79000) < 200,comparisonStructAll{1}))