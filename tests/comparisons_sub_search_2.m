
q = randi(10,1,10);
d = randi(10,1,20);

z = ones(1,10);




import SignalRegistration.unmasked_pcc_corr;
comparisonFun = @(x,y,z) unmasked_pcc_corr(x,y,z);

out = comparisonFun(q,d,z);

comparisonFun = @(x,y,z) unmasked_MASS_PCC(y,x,z,2^(4+nextpow2(length(x))));

out2 = comparisonFun(q,d,z);