function [mpd,mpdI,mpdII,mpABBA,mpIABBA] = mp_dist_stomp(X1,X2, w, kk,par)

    %   Args:
    %           X1
    %           X2
    %           w
    %
    %   Returns:
    %       mpd
    n = length(X1);
    m = length(X2);

    import mp.mp_profile_stomp_dna;
    

    % P_AB
    [mpAB,mpIAB] = mp_profile_stomp_dna(X1, X2, w,kk);

    % P_BA
    [mpBA,mpIBA] = mp_profile_stomp_dna(X2, X1, w, kk);

    % P_ABBA
    mpABBA = [mpAB;mpBA];
    % I_AMMA
    mpIABBA = [mpIAB;mpIBA];

    % k, specific value.
    % Don't take the highest value, but the one close to the top values..
    k = par*(n+m);

    if length(mpABBA) > k
        [a,~] = sort(mpABBA,'desc');
        mpd = a(round(k));
        [~,mpdI] = max(mpABBA);
        mpdII = mpIABBA(mpdI);
    else
        [mpd,mpdI] = max(mpABBA); % max, gives the highest pcc as distance between the two, also would be nice to have position
        mpdII = mpIABBA(mpdI);
    end






end

