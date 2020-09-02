% % chech if the sequences have to be reproducible
rng(1,'twister');

    
% alignment test simulated - we simulate the fragmented barcodes. Rather
% easy procedure as we just cut long barcode at random places
lengthN = 8000;
sigma = 2.3;

pcc = @(x,y) 1./length(x)*zscore(x,1)'*zscore(y,1);

% generate random vector for the chromosome
barN = normrnd(0,1,1,lengthN);
chromosomeBar = imgaussfilt(barN,sigma,'Padding','circular');
% want

fig=figure;
hAxis = subplot(1,1,1);
import CBT.Hca.UI.Helper.plot_single_concentric;
plot_single_concentric({chromosomeBar}, hAxis,'Chromosome example')