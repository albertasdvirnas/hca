function [ sets ] = get_user_theory_sets( sets )
    

    if ~sets.promptfortheory    
        try 
            fid = fopen(sets.fastas); 
            fastaNames = textscan(fid,'%s','delimiter','\n'); fclose(fid);
            for i=1:length(fastaNames{1})
                [FILEPATH,NAME,EXT] = fileparts(fastaNames{1}{i});

                sets.theoryNames{i} = strcat(NAME,EXT);
                sets.theoryFold{i} = FILEPATH;
            end
        catch
            sets.promptfortheory  = 1;
        %         error('No valid fasta provided, please check the provided file');
        end
    end
    
    if sets.promptfortheory    
        % loads figure window
        import Fancy.UI.Templates.create_figure_window;
        [hMenuParent, tsHCA] = create_figure_window('Theory fasta(s) import tool','HCA');

        import Fancy.UI.Templates.create_import_tab;
        cache = create_import_tab(hMenuParent,tsHCA,'theory fasta(s)');
        uiwait(gcf);  
        
        dd = cache('selectedItems');
        sets.theoryNames = dd(1:end/2);
        sets.theoryFold = dd((end/2+1):end);
        delete(hMenuParent);

    end
    
    if sets.promptforsavetheory
       sets.resultsDir =  uigetdir(pwd,'Pick a directory where to save the theory results');
    end

end

