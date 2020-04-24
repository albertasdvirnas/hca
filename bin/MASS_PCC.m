% This code is created by Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan, Chetan Kumar Gupta and Eamonn Keogh.
% The overall time complexity of the code is O(n log n). The code is free to use for research purposes.
% The code may produce imaginary numbers due to numerical errors for long time series where batch processing on short segments can solve the problem.

%x is the long time series
%y is the query
%k is the size of pieces, preferably a power of two

function [dist] = MASS_PCC(x, y, k)
    % MASS_PCC - batch processing of PCC on short segments 
    %
    %   Args:
    %       x, y, k
    %
    %   Returns:
    %       dist, distance matrix
    
    %x is the data, y is the query
    m = length(y);
    n = length(x);
    dist = zeros(2,n-m+1);

    %compute y stats -- O(n)
    y = zscore(y,1); % zcore(y) if want to use std with /(m-1)
%     meany = mean(y);
%     sigmay = std(y,1);

    %compute x stats -- O(n)
%     meanx = movmean(x,[m-1 0]);
    sigmax = movstd(x,[m-1 0],1); % normalizes by m

    %k = 4096; %assume k > m
    %k = pow2(nextpow2(sqrt(n)));
    y2 = y;
    y = y(end:-1:1); %Reverse the query
    y(m+1:k) = 0; %append zeros
    y2(m+1:k) = 0; %append zeros

    Y = fft(y);
    Y2 = fft(y2);
    for j = 1:k-m+1:n-k+1
        %The main trick of getting dot products in O(n log n) time
        X = fft(x(j:j+k-1));

        Z = X.*Y;
        z = ifft(Z)./m;
        dist(1,j:j+k-m) = z(m:k)./(sigmax(m+j-1:j+k-1));
        
        Z = X.*Y2;
        z = ifft(Z)./m;
        dist(2,j:j+k-m) = z(m:k)./(sigmax(m+j-1:j+k-1));

    end

    if isempty(j)
        j = 0; % if nothing was computed
        k = n;
    else
        j = j+k-m;
        k = n-j; % number of points left
    end
    
    if k >= m % if k < m, there are not enough points on long barcode to compute more PCC's
        
        %The main trick of getting dot products in O(n log n) time
        X = fft(x(j+1:n));

        y(k+1:end)= [];

        Y = fft(y);
        
        Z = X.*Y;
        z = ifft(Z)./m;

        dist(1, j+1:n-m+1) = z(m:k)./(sigmax(j+m:n));
        
        y2(k+1:end)= [];

        Y = fft(y2);
        Z = X.*Y;
        z = ifft(Z)./m;
        
        dist(2,j+1:n-m+1) = z(m:k)./(sigmax(j+m:n));
    end
    
end