function [lm, cache] = launch_import_ui( ...
  hMenuParent, ...
  hPanelKymoImport, ...
  tsHCA, ...
  tabTitle, ...
  cache)
if nargin < 5
  cache = containers.Map();
end

import Fancy.UI.FancyList.FancyListMgr;
lm = FancyListMgr();
lm.set_ui_parent(hPanelKymoImport);
lm.make_ui_items_listbox();

import Fancy.UI.FancyList.FancyListMgrBtnSet;
% flmbs1 = FancyListMgrBtnSet();
% flmbs1.NUM_BUTTON_COLS = 2;
% flmbs1.add_button(FancyListMgr.make_select_all_button_template());
% flmbs1.add_button(FancyListMgr.make_deselect_all_button_template());
flmbs2 = FancyListMgrBtnSet();
flmbs2.NUM_BUTTON_COLS = 2;
flmbs2.add_button(make_add_movies_directly_btn(tsHCA));
flmbs2.add_button(make_remove_movies_btn());

flmbs3 = FancyListMgrBtnSet();
flmbs3.NUM_BUTTON_COLS = 1;
flmbs3.add_button(extract_movies_from_list());

lm.add_button_sets(flmbs2,flmbs3);

  function [btnAddKymos] = make_add_movies_directly_btn(ts)
    import Fancy.UI.FancyList.FancyListMgrBtn;
    btnAddKymos = FancyListMgrBtn(...
      ['Add ', tabTitle], ...
      @(~, ~, lm) on_add_movies_directly(lm, ts));
    
    function [] = on_add_movies_directly(lm, ~)
      
      [itemFilenames, itemDirpath] = uigetfile( ...
        '*.*', ...
        strcat(['Select ' tabTitle ' file(s) to import']), ...
        pwd, ...
        'MultiSelect','on');
      
      if ~iscell(itemFilenames)
        itemFilenames = {itemFilenames};
      end
      
      lm.add_list_items(itemFilenames', repmat({itemDirpath},length(itemFilenames),1));
      
    end
  end

  function [btnRemoveKymos] = make_remove_movies_btn()
    import Fancy.UI.FancyList.FancyListMgrBtn;
    btnRemoveKymos = FancyListMgrBtn(...
      strcat(['Remove selected ' tabTitle]), ...
      @(~, ~, lm) on_remove_selected_movies(lm));
    function [] = on_remove_selected_movies(lm)
      lm.remove_selected_items();
    end
  end

  function [btnRemoveKymos] = extract_movies_from_list()
    import Fancy.UI.FancyList.FancyListMgrBtn;
    btnRemoveKymos = FancyListMgrBtn(...
      strcat(['Import ' tabTitle ' from list']), ...
      @(~, ~, lm) on_extract_movies_from_list(lm));
    function [] = on_extract_movies_from_list(lm)
      cache('selectedItems') = get_all_list_items(lm);
      delete(hMenuParent);
      uiresume(gcf);
    end
  end


end