import SignalRegistration.unmasked_pcc_corr;
comparisonFun1 = @(x,y,z) unmasked_pcc_corr(x,y,z);


comparisonFun2 = @(x,y,z) unmasked_MASS_PCC(y,x,z,2^(4+nextpow2(length(x))));

N = 100;
bar = rand(1,N);
bit =zeros(1,N);
bit(2:end-4) = 1;

M = 10000001;
ref = rand(1,M);

[out] = comparisonFun1(bar,ref,bit);
[out2] = comparisonFun2(bar,ref,bit);

figure,plot((out(:)-out2(:)).^2)

xlabel('Position')
ylabel('Square difference (scoreNew-scoreOld)^2')
title('Comparison of unmasked_mass_pcc vs unmasked_pcc_corr','Interpreter','latex')
