%% Here we test some suspicious data from theory test to track down the errors
%%
load('data1.mat');
load('data2.mat');
% figure,plot(data)
% hold on,plot(data2)

figure,plot(data-data2')

%%
load('tempt1.mat');
load('tempt2.mat');
n = 1000;
figure,plot(tempTheory(1:n)-xtraseq(1:n)')
