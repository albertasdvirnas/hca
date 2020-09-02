len1 = 500;
len2 = 10000;

rng(1,'twister');

islinear = 1;
if islinear
    bar2 = imgaussfilt(normrnd(0,1,len2,1),2.3);
else
    
end

randPos = randi(len2-len1);

bar1 = bar2(randPos:randPos+len1-1);

barR = imgaussfilt(normrnd(0,1,len1,1),2.3);

% add noise
pcc = 0.7;
distFun = @(a,b) (zscore(a,1))'*(zscore(b,1))/length(b);

fun = @(x) ((zscore(bar1,1))'*(zscore((1-x)*zscore(bar1,1)+x*zscore(barR,1),1))/length(bar1))-pcc;
try
    sol = fzero(fun, [0 1]);
catch
    sol = 0;
end

% fun(sol)
bar1 = (1-sol)*zscore(bar1,1)+sol*zscore(barR,1);
% distFun(bar1, bar2(randPos:randPos+len1-1))

% distFun(bar1, bar2(randPos:randPos+len1-1))
bit1 = ones(len1,1);

r = 200;
k= 2^14;


import mp.mp_masked_profile_stomp_dna;
comparisonFun = @(x,y,z,w) mp_masked_profile_stomp_dna(x,y,z,w,2^(4+nextpow2(length(x))));

tic
[mp, mpI] = comparisonFun(bar1,bar2,bit1,r);
toc

max(mp)
% barBest = bar2(mpI:mpI+r-1);
% distFun(bar1,barBest)


isLinearTF =1;
mask = 10;
import CBT.Hca.UI.Helper.get_best_parameters_mp;
[ maxcoef,pos,or ] = get_best_parameters_mp( mp,mpI, 10, len1,len2,  isLinearTF, mask);

%could base the significance on how many of the top scores are 
