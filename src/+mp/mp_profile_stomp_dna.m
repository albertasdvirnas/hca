%     bTS = rand(1000,1);
%     aTS = rand(100,1);
%     r = 50;
%     kk =128;
% [mp,mpI] = mp_profile(aTS, bTS,r,kk);

function [mp,mpI] = mp_profile_stomp_dna(aTS, bTS,r,kk)

    mpLength = length(aTS)-r+1;
    mpColLength = length(bTS)-r+1;

	mp = zeros(mpLength,1); % instead of making twice the size, report a 
    mpI = zeros(mpLength,1); % a negative number in mpI? Then all procedure simplifies.
    mp2 = zeros(mpLength,1); % instead of making twice the size, report a 
    mpI2 = zeros(mpLength,1); % a negative number in mpI? Then all procedure simplifies.

    % forward
    [mp,mpI] = mp_profile_stomp(aTS, bTS,r,kk);
    
    % reverse. Here aTS is flipped, so mp(end) is actuall mp(1), mp(end-1)
    % mp(2), etc. Could pass mp as input ?
    [mp2(end:-1:1),mpI2(end:-1:1)] = mp_profile_stomp(aTS(end:-1:1), bTS,r,kk);
    
	updatePos(:,1) = mp(:) < mp2;
    mp(updatePos) = mp2(updatePos);
    mpI(updatePos) = mpI2(updatePos) + mpColLength;
    

    % combine

    
end

function [mp,mpI] = mp_profile_stomp(aTS, bTS,r,kk)
    % simple method to compute matrix profile for time series data comping
    % from DNA barcoding
    %
    % Args:
    %   aTS, bTS, r
    %
    % Returns:
    %   mp: matrix profile,
    %   mpI: matrix profile index
    
    mpLength = length(aTS)-r+1;
    mpColLength = length(bTS)-r+1;

    % we first loop through all subsequences of aTS of length r
%     import mp.unmasked_pcc_corr;
%     import mp.MASS_PCC;
    
    [X, sumx2, sumx, meanx, sigmax2, sigmax] = fastfindNNPre(bTS, r);
    [Y, sumy2, sumy, meany, sigmay2, sigmay] = fastfindNNPre(aTS, r);

%     sigmax'
%     [Zfft, sumz2, sumz, meanz, sigmaz2, sigmaz] = ...
%         fastfindNNPre(aTS, r);
    qRow = zeros(mpLength,1); % same as MP
    qCol = zeros(mpColLength,1); % columns distance profile

    % Start from
    dist = FAST_CC(bTS,aTS(1:r), kk );
    dist2 = FAST_CC(aTS,bTS(1:r), kk );

    qCol(:)= ((dist./r)-meanx*meany(1))./(sigmax.*sigmay(1));
%     qCol(end/2+1:end)= ((dist(end/2+1:end)./r)-meanx*meany(1))./(sigmax.*sigmay(1));

    qRow(:) = ((dist2./r)-meanx(1)*meany)./(sigmax(1).*sigmay);
%     qRow((end/2+1):end) = ((dist2((end/2+1):end)./r)-meanx(1)*meany)./(sigmax(1).*sigmay);

    
    mp = zeros(mpLength,1); % instead of making twice the size, report a 
    mpI = zeros(mpLength,1); % a negative number in mpI? Then all procedure simplifies.
    updatePos=false(mpLength,1);

    % evaluate initial matrix profile
    mp(:) = qRow;
    mpI(:) = 1;
    
    updatePos(:,1) = mp(:) < qRow;
    mp(updatePos) = qRow(updatePos);
    mpI(updatePos) = mpColLength+1;
    [ mp(1), mpI(1)] = max(qCol);
%     figure,plot(qCol)
    
%     [ mp(mpLength+1), mpI(mpLength+1)] = max(qCol(mpLength+1:end));

    for i = 2:mpLength
        % update values
%         if i == 3
%             dist(mpColLength+1)
%         end
        dist(2:mpColLength) = dist(1:mpColLength-1)- bTS(1:mpColLength-1).*aTS(i-1)+bTS(r+1:end).*aTS(i+r-1);

%         dist((mpColLength+2):end) = dist((mpColLength+1):(end-1))- bTS(1:(mpColLength-1)).*aTS(i+r-2)- bTS(r+1:end).*aTS(i)+bTS(r+1:end).*aTS(i+1)+bTS(2:(mpColLength)).*aTS(i+r-1);
        
        
%         dist(mpColLength+1) - bTS(r)*aTS(1)-bTS(1)*aTS(r)+bTS(2)*aTS(r+1)+bTS(r+1)*aTS(2)
%         sum(bTS(2:r+1).*flipud(aTS(2:r+1)))
%         dist(mpColLength+2)

        % update first values
        dist(1) = dist2(i);
%         dist(mpColLength+1) = dist2(mpLength+i);

        % test: first ellements are comp correct
%          sum(bTS(1:r).*flipud(aTS(2:r+1)))
%         dist(mpColLength+1)
%       sum(bTS(2:r+1).*aTS(2:r+1))
%         dist(2)
% correct
   %         dist(mpColLength+2)
        % correct

        
        % incorrect
%         sum(bTS(2:r+1).*flipud(aTS(2:r+1)))
%         dist(mpColLength+2)
     

        
%         sum(bTS(2:r+1).*flipud(aTS(2:r+1)))
%         dist(mpColLength+2)
        % ((dist(831)-bTS(831).*aTS(1)+bTS(832+r)*aTS(1+r) )/r-meanx(832)*meany(2))/(sigmax(832)*sigmay(2))

        qCol= ((dist./r)-meanx*meany(i))./(sigmax.*sigmay(i));
%         qCol((end/2+1):end)= ((dist((end/2+1):end)./r)-meanx*meany(i))./(sigmax.*sigmay(i));
        %

%         transpose(((dist((end/2+1):end)./r)-meanx*meany(i))./(sigmay(i)))
%         if i == 2
%             dist2(mpLength+i)
% 
%             figure,plot(qCol);hold on;
%         end
            
%         updatePos(:,1) = mp(:) < qRow(end/2+1:end);

%         [a, b] =  max(qCol);
%         if 

        [ mp(i), mpI(i)] = max(qCol);
%         [ mp(mpLength+i), mpI(mpLength+i)] = max(qCol(mpLength+1:end));       
    end
end


function [dist] = FAST_CC(x, y, k)
    % FAST CC- batch processing of CC on short segments 
    %
    %   Args:
    %       x, y, k
    %
    %   Returns:
    %       dist, distance matrix
    
    %x is the data, y is the query
    m = length(y);
    n = length(x);
    mpLen = n-m+1;
    dist = zeros(mpLen,1);

    %compute y stats -- O(n)
%     y = zscore(y,1); % zcore(y) if want to use std with /(m-1)
%     meany = mean(y);
%     sigmay = std(y,1);

    %compute x stats -- O(n)
%     meanx = movmean(x,[m-1 0]);
%     sigmax = movstd(x,[m-1 0],1); % normalizes by m

    %k = 4096; %assume k > m
    %k = pow2(nextpow2(sqrt(n)));
%     y2 = y;
    y = y(end:-1:1); %Reverse the query
    y(m+1:k) = 0; %append zeros
%     y2(m+1:k) = 0; %append zeros

    Y = fft(y);
%     Y2 = fft(y2);

    for j = 1:abs(k-m+1):n-k+1

        %The main trick of getting dot products in O(n log n) time
        X = fft(x(j:j+k-1));

        Z = X.*Y;
        z = ifft(Z);
        dist(j:j+k-m) = z(m:k);%./(sigmax(m+j-1:j+k-1));
        
%         Z = X.*Y2;
%         z = ifft(Z);
%         dist((mpLen+j):(mpLen+j+k-m)) = z(m:k);%./(sigmax(m+j-1:j+k-1));

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
        z = ifft(Z);

        dist(j+1:n-m+1) = z(m:k);%./(sigmax(j+m:n));
        
%         y2(k+1:end)= [];
% 
%         Y = fft(y2);
%         Z = X.*Y;
%         z = ifft(Z);
% 
%         dist((mpLen+j+1):(mpLen+n-m+1)) = z(m:k);%./(sigmax(j+m:n));

    end
%     dist = dist./m;
    
end


% m is winSize
function [X, sumx2, sumx, meanx, sigmax2, sigmax] = fastfindNNPre(x, m)
    n = length(x);
    x(n+1:2*n) = 0;
    X = fft(x);
    cum_sumx = cumsum(x);
    cum_sumx2 =  cumsum(x.^2);
    sumx2 = cum_sumx2(m:n)-[0;cum_sumx2(1:n-m)];
    sumx = cum_sumx(m:n)-[0;cum_sumx(1:n-m)];
    meanx = sumx./m;
    sigmax2 = (sumx2./m)-(meanx.^2);
    sigmax = sqrt(sigmax2);
end
