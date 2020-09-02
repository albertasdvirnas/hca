% time series
numPts = 200;
A = imgaussfilt(rand(1,numPts),2.3);
A =[ A A A];

numPts = length(A);

% embed time series into point could using Taken's embedding Theorem

d = 2;
% t = 2;

% critical bound
criticalBound = 2/sqrt(numPts);

% smallest time lag h where the sample autocorrelation function (ACF) ˆρh 
% becomes insignificant, i.e., smaller in absolute value than the critical bound
t = find(autocorr(A) < criticalBound,1,'first')-1;

% map to point cloud

[PC] = map_pc(A,t,d,numPts);
% 

figure,plot(PC(:,1),PC(:,2),'x')
% % x_i x_{i+tau}
% PC = zeros(numPts-t,d);
% i = 1:numPts-t;
% for j=i
%     PC(j,:) = [A(j) A(j+t)];
% end