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
    
    askfornm = sets.theory.askfornmbp;

     if ~sets.theory.askfortheory
          try 
            fid = fopen(sets.theories); 
            fastaNames = textscan(fid,'%s','delimiter','\n'); fclose(fid);
            for i=1:length(fastaNames{1})
                [FILEPATH,NAME,EXT] = fileparts(fastaNames{1}{i});

                sets.theoryFile{i} = strcat(NAME,EXT);
                sets.theoryFileFold{i} = FILEPATH;
            end
        catch
            sets.theory.askfortheory = 1;
          end
          
     end
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
    
   [ fd,fr,fl ] = fileparts(sets.theoryFile{1});
    if ~isequal(fl,'.mat')&&~isequal(fl,'.txt')
        outdirpath = sets.output.matDirpath;

        % want this to use all the theory files and fols
        [t,matFilepathShort,theoryStruct, sets,theoryGen] = HCA_theory_parallel('',0.3,'theory_settings_parallel.txt',sets.theoryFile,sets.theoryFileFold);
         sets.output.matDirpath = outdirpath;
    else
        % now load theory
        import CBT.Hca.UI.Helper.load_theory;
        theoryStruct = load_theory(sets);
    end
    try delete(hMenuParent); catch; end;

    
    if askfornm
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
    theoryStruct = convert_nm_ratio(sets.theory.nmbp, theoryStruct,sets );
    toc
    
	datetime=datestr(now);
    datetime=strrep(datetime,' ','');%Replace space with underscore
    datetime=strrep(datetime,':','_');%Replace space with underscore
    
    if ispc
        name = strcat([sets.theoryFileFold{1} 'theoryStruct_' datetime '.mat']);
    else
        name = strcat([sets.theoryFileFold{1} '/theoryStruct_' datetime '.mat']);        
    end
    
    save(name, '-v7.3','theoryStruct');
    
    
end

