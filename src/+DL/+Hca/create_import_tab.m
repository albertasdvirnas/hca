function [cache] = create_import_tab(~, tsHCC, tabTitle, cache)
% create_import_tab
if nargin < 4
  cache = containers.Map();
end

% create main tab for the analysis
hTabKymoImport = tsHCC.create_tab(strcat([tabTitle ' import tab']));
tsHCC.select_tab(hTabKymoImport);
hPanelKymoImport = uipanel(hTabKymoImport);

% import kymographs
import DL.Hca.launch_import_ui;
[~, cache] = launch_import_ui( ...
  hTabKymoImport, ...
  hPanelKymoImport, ...
  tsHCC, ...
  tabTitle, ...
  cache);
end

