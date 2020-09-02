function [] = run_dot_vs_theory_similarity_analysis(tsHCC)
    % run_dot_vs_theory_similarity_analysis
    % input tsHCC

    % create main tab for the analysis
    tabTitle = 'Dot vs theory similarity analysis';
    hTabKymoImport = tsHCC.create_tab(tabTitle);
    tsHCC.select_tab(hTabKymoImport);
    hPanelKymoImport = uipanel(hTabKymoImport);

    % import enzyme sequence
    import CBT.Hca.UI.launch_kymo_import_ui;
    lm = launch_kymo_import_ui(hPanelKymoImport, tsHCC);
    
    
    import CBT.Hca.UI.launch_kymo_import_ui;
    lm = launch_kymo_import_ui(hPanelKymoImport, tsHCC);

    % load settings and kymo structure
    import CBT.Hca.UI.load_settings_and_session_structure;
    cache = load_settings_and_session_structure(lm, tsHCC);

    % make consensus structure
    import CBT.Hca.UI.make_consensus_structure;
    cache = make_consensus_structure(lm, tsHCC, cache);

    % add theory and compare
    import CBT.Hca.UI.add_theory_and_compare;
    cache = add_theory_and_compare(lm, tsHCC,cache);

end