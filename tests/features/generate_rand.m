function [bar1,bar2,randPos] = generate_rand(len1, len2, pcc, sigma, islinear)
    islinear = 1;

    rng(1,'twister');

    if islinear
        bar2 = imgaussfilt(normrnd(0,1,len2,1),2.3);
    else

    end

    randPos = randi(len2-len1);

    bar1 = bar2(randPos:randPos+len1-1);

    barR = imgaussfilt(normrnd(0,1,len1,1),sigma);

    % add noise
%     distFun = @(a,b) (zscore(a,1))'*(zscore(b,1))/length(b);

    fun = @(x) ((zscore(bar1,1))'*(zscore((1-x)*zscore(bar1,1)+x*zscore(barR,1),1))/length(bar1))-pcc;
    try
        sol = fzero(fun, [0 1]);
    catch
        sol = 0;
    end

    % fun(sol)
    bar1 = (1-sol)*zscore(bar1,1)+sol*zscore(barR,1);

end

