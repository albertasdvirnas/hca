%About correction factor for Euclidean distances
lengthS = 50:50:800;
lengthL = 500000;

Y= imgaussfilt(normrnd(0,1,1,lengthL),3);

n = 100;
distmax = zeros(length(lengthS),n);
edmax = zeros(length(lengthS),n);
for idy = 1:n;
    idy
    for idx=1:length(lengthS)
        X = imgaussfilt(normrnd(0,1,1,lengthS(idx)),3);
        comparisonFun = @(x,y,z) unmasked_MASS_PCC(y,x,z,2^(4+nextpow2(length(x))));

        [dist] = comparisonFun(X,Y,ones(1,length(X)));

        distmax(idx,idy) = max(dist(:));
        edmax(idx,idy) = sqrt(2*lengthS(idx)*(1-max(dist(:))));
    end
end

% % mean(distmax,2)
% figure,plot(lengthS,mean(distmax,2))
% curve = mean(edmax,2)';
% figure,plot(lengthS,curve)
% hold on
% plot(lengthS,curve./(lengthS))
% plot(lengthS,curve./(sqrt(lengthS)))
% legend({"ED","1/l normalized","1/sqrt(l) normalized" },'location', 'best')
% xlabel("l");
% ylabel("Distance")

%
% figure,plot(lengthS,mean(distmax,2))
curve = mean(edmax,2)';
curve2 = curve./(lengthS);
curve3 = curve./(sqrt(lengthS))
curve4 = mean(distmax,2)';
curve5 = mean(distmax,2)'.*lengthS.^(0.43);

% figure,plot(lengthS,curve)
% hold on
% plot(lengthS,curve2)
% plot(lengthS,curve3)
% legend({"ED","1/l normalized","1/sqrt(l) normalized" },'location', 'best')
% xlabel("l");
% ylabel("Distance")

curveN = curve./max(curve(:));
curve2N = curve2./max(curve2(:));
curve3N = curve3./max(curve3(:));
curve4N = curve4./max(curve4(:));
curve4N = curve4./max(curve4(:));
curve5N = curve5./max(curve5(:));

figure,plot(lengthS,curveN)
hold on
plot(lengthS,curve2N)
plot(lengthS,curve3N)
plot(lengthS,curve4N)
plot(lengthS,curve5N)

legend({"ED","1/l normalized","1/sqrt(l) normalized","PCC" ,'extra'},'location', 'best')
xlabel("l,length of the fragment");
ylabel("Normalized distance score")
title("Distance scores comparison for different normalized functions. # num samples = 100, N=500000");


%% instead run this for matrixprofile

% If use consensus
settings = 'C_settings.txt';

% load settings
sets = ini2struct( settings );


lenD = 5000;
lenQ = 5000;
import Pvalue.compute_hmm_p_value;
[queryfile, datafile] = compute_hmm_p_value(lenD,lenQ,[],sets.rand);
w = 50;
window = num2str(w);
pathToScript = '/home/albyback/git/sv/TEST/Paper_Scripts/A_data_plasmid_theories/scamp_script.sh';  % assumes script is in curent directory
 cmdStr       = [pathToScript ' ' datafile ' ' queryfile ' ' window];
 tic
system(cmdStr);
time = toc;
disp("Computed maximum scores");


  % compute pcc values
file1 =  strcat( strrep(queryfile,'pval','resultData'), 'index');
fileID = fopen(file1,'r'); formatSpec = '%f'; mpI = fscanf(fileID,formatSpec);  fclose(fileID);
file1 = strcat( strrep(queryfile,'pval','resultData'), 'column');
fileID = fopen(file1,'r'); formatSpec = '%f'; mp = fscanf(fileID,formatSpec);  fclose(fileID);
