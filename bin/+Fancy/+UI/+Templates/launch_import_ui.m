function [lm,cache] = launch_movie_import_ui(hMenuParent, hPanelKymoImport, tsHCA, tabTitle, cache)
    if nargin < 5
        cache = containers.Map();
    end
    
    % launch_movie_import_ui -
    %   adds a tab with list management UI/functionality  for
    %    movies

    import Fancy.UI.FancyList.FancyListMgr;
    lm = FancyListMgr();
    lm.set_ui_parent(hPanelKymoImport);
    lm.make_ui_items_listbox();
    
    import Fancy.UI.FancyList.FancyListMgrBtnSet;
    flmbs1 = FancyListMgrBtnSet();
    flmbs1.NUM_BUTTON_COLS = 2;
    flmbs1.add_button(FancyListMgr.make_select_all_button_template());
    flmbs1.add_button(FancyListMgr.make_deselect_all_button_template());
	flmbs2 = FancyListMgrBtnSet();
    flmbs2.NUM_BUTTON_COLS = 2;
    flmbs2.add_button(make_add_movies_directly_btn(tsHCA));
    flmbs2.add_button(make_remove_movies_btn());

    flmbs3 = FancyListMgrBtnSet();
    flmbs3.NUM_BUTTON_COLS = 1;
    flmbs3.add_button(extract_movies_from_list());

    lm.add_button_sets(flmbs1,flmbs2,flmbs3);

    function [btnAddKymos] = make_add_movies_directly_btn(ts)
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnAddKymos = FancyListMgrBtn(...
            strcat(['Add ' tabTitle ' directly']), ...
            @(~, ~, lm) on_add_movies_directly(lm, ts));
        
        function [] = on_add_movies_directly(lm, ts)

            [rawMovieFilenames, rawMovieDirpath] = uigetfile(...
            {'*.tif;*.tiff;*.mat;*.fasta;*.fa;*.txt'}, ...
            strcat(['Select raw ' tabTitle ' file(s) to import']), ...
            pwd, ...
            'MultiSelect','on');
            
            if ~iscell(rawMovieFilenames)
                rawMovieFilenames = {rawMovieFilenames};
            end
            
           lm.add_list_items(rawMovieFilenames', repmat({rawMovieDirpath},length(rawMovieFilenames),1));

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
            strcat(['Extract ' tabTitle ' from list']), ...
            @(~, ~, lm) on_extract_movies_from_list(lm));
        function [] = on_extract_movies_from_list(lm)
            [selectedItems, ~] = get_selected_list_items(lm);
            cache('selectedItems') = selectedItems;
            delete(hMenuParent);
            uiresume(gcf); 
        end
    end
    

end