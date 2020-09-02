function [barcodeGen] = gen_bac(data,sets)

listing = dir(fullfile(data,'*.tif'));

fd =fopen('tifs.txt','w');
for i=1:length(listing)
    fprintf(fd, '%s\n', fullfile(listing(i).folder,listing(i).name));
end
fclose(fd);

% load all user selected settings
import CBT.Hca.Settings.get_user_settings;
sets = get_user_settings(sets);

% add kymographs
import CBT.Hca.Import.add_kymographs_fun;
[kymoStructs] = add_kymographs_fun(sets);

%  put the kymographs into the structure
import CBT.Hca.Core.edit_kymographs_fun;
kymoStructs = edit_kymographs_fun(kymoStructs,sets.timeFramesNr);

% align kymos
import CBT.Hca.Core.align_kymos;
[kymoStructs] = align_kymos(sets,kymoStructs);

% generate barcodes
import CBT.Hca.Core.gen_barcodes;
barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoStructs, sets);

        
       
end

