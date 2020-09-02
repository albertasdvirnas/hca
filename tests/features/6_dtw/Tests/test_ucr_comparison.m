
addpath(genpath('Development/hca'));
[hcaStruct] = HCA_run_parallel('comparison_dtw.txt','theoryData/theories_2020-04-21_14_45_11_.txt' , '/home/albyback/workData/rawData/HumanDNAProject/BACs/Kymographs/P18/',10);


hcaStruct.comparisonStruct

coefs = cellfun(@(x) x.ucr.maxcoef, hcaStruct.comparisonStructAll{1});
pos = cellfun(@(x) x.ucr.pos, hcaStruct.comparisonStructAll{1});

figure,plot(coefs,'x')
figure,plot(pos,'x')

% hcaStruct.ccoefs = cellfun(@(x) x.ucr.maxcoef, hcaStruct.comparisonStructAll{1});
omparisonStructAll{1}{52}.ucr


[hcaStruct] = HCA_run_parallel('comparison_pcc.txt','theoryData/theories.txt' , '/home/albyback/workData/rawData/HumanDNAProject/BACs/Kymographs/P18/',1);
[hcaStruct] = HCA_run_parallel('comparison_dtw.txt','theoryData/theories.txt' , '/home/albyback/workData/rawData/HumanDNAProject/BACs/Kymographs/P18/',1);

% [hcaStruct] = HCA_run_parallel('comparison_pcc.txt','theoryData/theories_2020-04-22_13_10_30_.txt' , '/home/albyback/workData/rawData/HumanDNAProject/BACs/Kymographs/P18/',1);
