% script for unit tests showing how MP profile works for HMM
% should have at least one unit test for each feature that is tested.

% result = runtests('mp_unit_test');

% preconditions

%% Test 1: position
len2 = 1000;
[pos] = mp_position_test(len2);
assert(pos == len2+1)

%% Test 2 position with bitmask
len2 = 1000;
[pos] = mp_position_test(len2,5);
assert(pos == len2+1)

%% Test 3 mp true position loops around
len2 = 1000; len1 = 550;
[pos] = mp_position_test2(len2,5,len1);
assert(pos == len2+len1-101+2)

%% Test 4 inverse
len2 = 1002; len1 = 577;
[pos] = mp_position_test3(len2,8,len1);
assert(pos == len2+len1-101+2)

% 1451


