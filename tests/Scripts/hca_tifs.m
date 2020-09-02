tifs = 'bars.txt';


try 
    fid = fopen(tifs); 
    tifNames = textscan(fid,'%s','delimiter','\n'); fclose(fid);
catch
    error('No valid tifs provided, please check the provided file');
end

kymoStructs = cell(1,length(tifNames{1}));

% loop over movie file folder
for idx = 1:length(tifNames{1})       
    kymoStructs{idx}.unalignedKymo = imread(tifNames{1}{idx});
    kymoStructs{idx}.name = tifNames{1}{idx};
end
