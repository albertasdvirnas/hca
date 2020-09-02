thyCurve_bpRes = rand(1,1000000);

tic
% m = 100;

h = images.internal.createGaussianKernel(1000, 2^14)';
   
k = 2^19;
n = length(thyCurve_bpRes);
m = length(h);
dist = zeros(1,n-m+1);

h = h(end:-1:1); %Reverse the query
h(m+1:k) = 0; %append zeros
    
    
for j = 1:k-m+1:n-k+1

    %The main trick of getting dot products in O(n log n) time
    X = fft(thyCurve_bpRes(j:j+k-1));
    Y = fft(h);

    Z = X.*Y;
    z = ifft(Z);
    dist(j:j+k-m) = z(m:k);
% 
%     Y = fft(y2);
%     Z = X.*Y;
%     z = ifft(Z)./m;
%     dist(2,j:j+k-m) = z(m:k)./(sigmax(m+j-1:j+k-1));

end
    
j = j+k-m;
k = n-j; % number of points left
    
if k >= m % if k < m, there are not enough points on long barcode to compute more PCC's

    %The main trick of getting dot products in O(n log n) time
    X = fft(thyCurve_bpRes(j+1:n));

    h(k+1:end)= [];

    Y = fft(h);
    Z = X.*Y;
    dist(j+1:n-m+1) = ifft(z(m:k));


end
    
toc