% To check that HCA theory is computed the same throughout the versions, we
% need to compare it on a few sequences of different lengths.

% Currently, compute_theory_barcode_old computes the theory

% It's best if all the tests are saved in one test file
testTheories = 'theoriestesting.txt';
import CBT.Hca.Import.import_names;
[names] = import_names(testTheories);


%% CBT HCA THEORY 4.0
% imports

% import settings
import CBT.Hca.Import.import_settings;
[sets] = import_settings('sets.txt');

% idx/can loop through these, or pass idx as user input
idx = 1;


tic;
% sets.method = sets.theoryGen.model;
import CBT.Hca.Core.Theory.compute_theory_barcode;
[theorySeq, header] = compute_theory_barcode(names{idx},sets);
toc;



%% CBT HCA THEORY 3.8.4


% This should be same as hca-v3.8.4
% [chr1,header] = create_memory_struct(sequence);
fasta = fastaread(names{idx});
disp('started generating theory barcode.. 0%')
tic
import CBT.Hca.Core.Theory.compute_hca_theory_barcode;
[ theorySeq2, bitmask, probSeq,theorSeq] = compute_hca_theory_barcode(fasta.Sequence,sets);
disp('done generating theory barcode... 100%')
toc
% note that the barcode's match, bar the edges (and might not be the same
% number of pixels. 
% figure,plot(theorySeq)
% hold on
% plot(theorySeq2)

figure,plot((theorySeq2-theorySeq))

%%
import CBT.Core.convert_bpRes_to_pxRes;
meanBpExt_pixels = sets.theoryGen.meanBpExt_nm / sets.theoryGen.pixelWidth_nm;
% x1=  convert_bpRes_to_pxRes(tempTheory',meanBpExt_pixels );
x2 =  convert_bpRes_to_pxRes(probSeq,meanBpExt_pixels );

fold ='/media/albyback/My Passport/DATA/chromosomes/all/*.fa';

listing = dir(fold);

settings = 'sets_theory.txt';
fastas = 'fastas.txt';
fd = fopen(fastas,'w');
for i=1:length(listing)
    fprintf(fd, '%s\n', fullfile(listing(i).folder,listing(i).name));
end
fclose(fd);


[theoryGen] = hca_theory_script( fastas,settings );


%% SIMPLE test to see how the FFT works on theory data

% note that we need to remove one point at the end
k = 2^10;
x = zeros(1,k);
y = zeros(1,k);

m = 2^6;
t = 100;
% x(m/2+1)=1;
x(100)=1;
x(200) = 1;
x(end)=1;
y(1:m) = images.internal.createGaussianKernel(3, m);
Y = fft(y);
X = fft(x);

% product
Z = X.*Y;

% inverser fourier
z = ifft(Z);

out = z(m:k);

z(t+m/2)
figure,plot(z)
