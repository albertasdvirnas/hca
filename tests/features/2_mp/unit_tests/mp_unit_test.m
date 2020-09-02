% script for unit tests showing how MP profile works for HMM
% should have at least one unit test for each feature that is tested.

% result = runtests('mp_unit_test');

% preconditions

%% Test 1: position
len2 = 1000;
bitmask = 0;
[pos] = mp_position_test(len2,bitmask,0);
assert(pos == len2+1)

%% Test 2 position with bitmask
len2 = 1000;
[pos] = mp_position_test(len2,5,0);
assert(pos == len2+1)

%% Test 3 mp true position loops around
len2 = 1000; len1 = 550;
[pos] = mp_position_test2(len2,5,len1,0);
assert(pos == len2+len1-101+2)

%% Test 4 inverse
len2 = 1002; len1 = 577;
[pos] = mp_position_test3(len2,8,len1,0);
assert(pos == len2+len1-101+2)

%% Test 5 long
len2 = 10000; len1 = 400;
[pos] = mp_position_test3(len2,8,len1,0);
assert(pos == len2+len1-101+2)

%% test52
len2 = 1002; len1 = 577;
[pos] = mp_position_test5(len2,8,len1,0);
assert(pos == len2+len1-101+2)

%% test6
len2 = 1002; len1 = 577;
[pos] = mp_position_test6(len2,8,len1,1);
assert(pos == -len1+2)

%% test 7
len2 = 1004; len1 = 540;
[pos] = mp_position_test6(len2,8,len1,0);
assert(pos == len2+3)

% bu tto test that it works we can also add pcc check!

%% testPCC
% want to test if pcc is in the correct place

len1 = 500;
len2 = 1000;
r = 100;
bar1 = imgaussfilt(normrnd(0,1,len1,1),2.3);
bar2 = imgaussfilt(normrnd(0,1,len2,1),2.3);
bar2(1:100) = bar1(2:101);
corr = mp_corr_test1(bar1,bar2,r,1);
assert(abs(corr-1)<0.0001)

%%
len1 = 500;
len2 = 1000;
r = 100;
bar1 = imgaussfilt(normrnd(0,1,len1,1),2.3);
bar2 = imgaussfilt(normrnd(0,1,len2,1),2.3);
bar2(1:100) = fliplr(bar1(2:101));
corr = mp_corr_test1(bar1,bar2,r,1);
assert(abs(corr-1)<0.0001)

%%
len1 = 500;
len2 = 1000;
r = 100;
bar1 = imgaussfilt(normrnd(0,1,len1,1),2.3);
bar2 = imgaussfilt(normrnd(0,1,len2,1),2.3);
bar2(301:400) = fliplr(bar1(2:101));
corr = mp_corr_test1(bar1,bar2,r,1);
assert(abs(corr-1)<0.0001)
%%
len1 = 500;
len2 = 1000;
r = 100;
bar1 = imgaussfilt(normrnd(0,1,len1,1),2.3);
bar2 = imgaussfilt(normrnd(0,1,len2,1),2.3);
bar2(301:400) = fliplr(bar1(2:101));
corr = mp_corr_test1(bar1,bar2,r,0);
assert(abs(corr-1)<0.0001)

%% check the correlation test..
len1 = 500;
len2 = 1000;
r = 100;
bar1 = imgaussfilt(normrnd(0,1,len1,1),2.3);
bar2 = imgaussfilt(normrnd(0,1,len2,1),2.3);
bit1 = ones(length(bar1),1);
% bit1(5:end-5) = 1;
bit2 = ones(length(bar2),1);
% bit2(7:end-7) = 1;

%%
% bar2(301:400) = fliplr(bar1(20:119));
% corr = mp_corr_test1(bar1,bar2,r,1);
% assert(abs(corr-1)<0.0001)

% 282&20
corr = mp_corr_test2(bar1,bar2,bit1, bit2, r,1);
assert(abs(corr-1)<0.0001)


bit1 = zeros(length(bar1),1);
bit1(5:end-1) = 1;
bit2 = ones(length(bar2),1);
corr = mp_corr_test2(bar1,bar2,bit1, bit2, r,1);
assert(abs(corr-1)<0.0001)
% seems ok
%%

bit1 = ones(length(bar1),1);
% bit1(5:end-1) = 1;
bit2 = zeros(length(bar2),1);
bit2(5:end-1) = 1;
corr = mp_corr_test2(bar1,bar2,bit1, bit2, r,1);
assert(abs(corr-1)<0.0001)
% bit2(7:end-7) = 1;


bit1 = zeros(length(bar1),1);
bit1(5:end-1) = 1;
bit2 = zeros(length(bar2),1);
bit2(5:end-1) = 1;
corr = mp_corr_test2(bar1,bar2,bit1, bit2, r,1);
assert(abs(corr-1)<0.0001)
