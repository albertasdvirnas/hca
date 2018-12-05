function [ theoryStruct,sets ] = get_user_theory( sets )
    % get_user_theory
    %
    % This function asks for all the required settings for the input of 
    % theory
    %
	%     Args:
    %         sets (struct): Input settings to the method
    % 
    %     Returns:
    %         theoryStruct: Return structure

    if sets.theory.askfortheory
        % loads figure window
        import Fancy.UI.Templates.create_figure_window;
        [hMenuParent, tsHCA] = create_figure_window('CB HCA theory selection tool','HCA');

        import Fancy.UI.Templates.create_import_tab;
        cache = create_import_tab(hMenuParent,tsHCA,'theory');
        uiwait(gcf);  
        
        dd = cache('selectedItems');
        sets.theoryFile = dd(1:end/2);
        sets.theoryFileFold = dd((end/2+1):end);
    end
    

    % now load theory
    import CBT.Hca.UI.Helper.load_theory;
    theoryStruct = load_theory(sets);
    
    if sets.theory.askfornmbp
        % here we add settings
        prompt = {'Nm/bp ratio for the experiments','Stretch factor', 'step'};
        title = 'Comparison to theory settings';
        dims = [1 35];
        definput = {'0.225','0.05','0.01'};
        answer = inputdlg(prompt,title,dims,definput);
        try
            sets.theory.nmbp  = str2double(answer{1}); % number of random cutouts from the input set
            sets.theory.stretchFactors = (1-str2double(answer{2})):str2double(answer{3}):(1+str2double(answer{2}));
        catch
            disp('Default nm/bp value is being assigned');
            sets.theory.nmbp = 0.225;
            sets.theory.stretchFactors = 0.95:0.01:1.05;
        end
    end

    tic
    import CBT.Hca.Core.Analysis.convert_nm_ratio;
    theoryStruct = convert_nm_ratio(sets.theory.nmbp, theoryStruct );
    toc
 
end

