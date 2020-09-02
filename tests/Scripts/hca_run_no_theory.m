
% folder with kymos, later make sure this can go through all the folder in
% dir.
% fold = '/media/albyback/My Passport/DATA/Shared - Kymographs for Albertas testing/SareusEcoli_50,50,G+,G-/RawKymographs 190125 - 130nm - 0.249 nm_bp /*.tif';
% fold = '/media/albyback/My Passport/DATA/Shared - Kymographs for Albertas testing/Pmirabilis_DA62690,Urin/RawKymographs 190109 - 130nm - 0.243 nm_bp/*.tif';
fold = '/media/albyback/My Passport/DATA/Shared - Kymographs for Albertas testing/Ecoli_BL21/RawKymographs 180420 - 130nm - 0.250 nm_bp/*.tif';
setsFile = 'setsBacterial.txt';

sets = ini2struct(setsFile);

listing = dir(fold);

tifs = 'tifsBacterial.txt';

fd = fopen(tifs,'w');
for i=1:length(listing)
    fprintf(fd, '%s\n', fullfile(listing(i).folder,listing(i).name));
end
fclose(fd);


[consensusStructs, kymoStructs,barcodeGen] = hca_no_theory( tifs, sets );