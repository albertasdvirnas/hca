function [zvals, cc] = compute_zval(thisBar, thisBit, theoryBar, theoryBit, barType, zeroModelParams, ccThresh, sets)
import SignalRegistration.masked_pcc_corr
import DL.Hca.beta_ev_cdf
import DL.Hca.trunc_normal_ev_cdf

lenBar = length(thisBar);
doBitmasking = sets.duallabel.doBitmask;
theoryNumPixels = sets.theory.totalNumBp / sets.theory.pixelWidth_nm * sets.theory.nmbp;
numStretchFactors = length(sets.theory.stretchFactors);
lambdaEff = zeroModelParams(1) * numStretchFactors * (theoryNumPixels - lenBar);

switch barType
  case 'dense'
    if doBitmasking
      cc = masked_pcc_corr(thisBar, theoryBar, thisBit, theoryBit);
      nuLenBar = sum(thisBit);
    else
      cc = MASS_PCC(theoryBar, thisBar, 2^(4 + nextpow2(lenBar)));
      nuLenBar = length(thisBar);
    end
    
    doZScore = cc > nanstd(cc(:)) * ccThresh;
    zvals = zeros(size(cc));
    zvals(doZScore) = -2 * log(vpa(1) - beta_ev_cdf(cc(doZScore), max(4, zeroModelParams(3) * nuLenBar), 1, lambdaEff));
    
    doIncreasedPrecision = isinf(zvals) & zvals > 0;
    zvals(doIncreasedPrecision) = -2 * log(vpa(1) - beta_ev_cdf(cc(doIncreasedPrecision), max(4, zeroModelParams(3) * nuLenBar), 1, lambdaEff, true));
  case 'sparse'
    cc = MASS_DOT_CC(theoryBar, thisBar, 2^(4 + nextpow2(lenBar)));
    
    doZScore = cc > nanmean(cc(:)) + nanstd(cc(:)) * ccThresh;
    zvals = zeros(size(cc));
    zvals(doZScore) = -2 * log(vpa(1) - trunc_normal_ev_cdf(cc(doZScore), norminv(linspace(0.001, 0.999, 31), zeroModelParams(3), 1 / sqrt(zeroModelParams(4) * lenBar)), 1 / sqrt(zeroModelParams(5) * lenBar), lambdaEff, 0, 1));
    
    doIncreasedPrecision = isinf(zvals) & zvals > 0;
    zvals(doIncreasedPrecision) = -2 * log(vpa(1) - trunc_normal_ev_cdf(cc(doIncreasedPrecision), norminv(linspace(0.001, 0.999, 51), zeroModelParams(3), 1 / sqrt(zeroModelParams(4) * lenBar)), 1 / sqrt(zeroModelParams(5) * lenBar), lambdaEff, 0, 1, true));
end

