function [x] = compute_theory_wrapper(ts, sets)
% theory wrapper.  Can compute multiple theory methods based on the
% input in the settings files

% Args:
%    ts : time series structures
%    sets : settings structure

import CBT.Hca.Core.Theory.run_simple_theory;
import CBT.Hca.Core.Theory.run_literature_theory;

%     if sum(ts > 4) > 0 % if there are nan's, don't compute the theory
%         x = nan(length(ts),1);
%     else
switch sets.theoryGen.method
  case "simple"
    x = run_simple_theory(ts);
  case "literature"
    x = run_literature_theory(ts',sets);
  case "TCGA"
    x = zeros(length(ts),1);
    dots = strfind(ts',sets.model.pattern);
    x(dots) = 1;
  case "enzyme"
    x = zeros(length(ts),1);
    [~, dots] = restrict(ts(:)', sets.theoryGen.pattern, 1);
    % Approximate normalization between 0 and 1
    x(dots(2:end)) = sets.theoryGen.pixelWidth_nm/sets.theoryGen.meanBpExt_nm;
  otherwise
    error('Incorrect method for theory generation selected' )
end
%     end


end

