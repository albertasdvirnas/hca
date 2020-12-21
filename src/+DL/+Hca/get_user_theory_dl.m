function [ theoryStruct,sets ] = get_user_theory_dl( sets )
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

tmpTheoryStruct = cell(1,2);
btype = {'CB', 'Dots'};

for i=1:2
  if not(sets.theory.askfortheory)
    try
      switch i
        case 1
          fid = fopen(sets.theories);
        case 2
          fid = fopen(sets.sparseTheories);
      end
      thrNames = textscan(fid,'%s','delimiter','\n');
      fclose(fid);
      for j=1:length(thrNames{1})
        [FILEPATH,NAME,EXT] = fileparts(thrNames{1}{j});
        sets.theoryFile{j} = strcat(NAME,EXT);
        sets.theoryFileFold{j} = FILEPATH;
      end
    catch
      sets.theory.askfortheory = 1;
    end
  end
  if sets.theory.askfortheory
    % loads figure window
    import Fancy.UI.Templates.create_figure_window;
    [hMenuParent, tsHCA] = create_figure_window([btype{i} ' HCA theory selection tool'],'HCA');
    
    import Fancy.UI.Templates.create_import_tab;
    cache = create_import_tab(hMenuParent, tsHCA, [btype{i} '-theory']);
    uiwait(gcf);
    
    dd = cache('selectedItems');
    sets.theoryFile = dd(1:end/2);
    sets.theoryFileFold = dd((end/2+1):end);
    
    % Sort by name, hopefully both barcode types have the same order....
    [sets.theoryFile, sortedIds] = sort(sets.theoryFile);
    sets.theoryFileFold = sets.theoryFileFold(sortedIds);
    
    delete(hMenuParent);
    
    if i==2
      answer = questdlg('What is the second type of theory?', ...
        'Choose label type', ...
        'Intensity profile', ...
        'Dots in pixels', ...
        'Intensity profile');
      if strcmp(answer, 'Dots in pixels')
        tmpTheoryStruct{i} = cell(1, length(sets.theoryFile));
        for j=1:length(sets.theoryFile)
          tmpTheoryStruct{i}{j}.filename = fullfile(sets.theoryFileFold{j}, ...
            sets.theoryFile{j});
        end
        break
      end
    end
  end
  
  % now load theory
  import CBT.Hca.UI.Helper.load_theory;
  tmpTheoryStruct{i} = load_theory(sets);
end

% assignin('base', 'tmpTheoryStruct', tmpTheoryStruct)

if not(length(tmpTheoryStruct{1}) == length(tmpTheoryStruct{2}))
  throw(MException('theoryImport:empty', ...
    'Theory pair member is empty, aborting theory import.'));
end

if sets.theory.askfortheory && strcmp(answer, 'Dots in pixels')
  import Microscopy.Simulate.Core.apply_point_spread_function
  for j=1:length(tmpTheoryStruct{1})
    thisDots = importdata(tmpTheoryStruct{2}{j}.filename);
    dotbar = zeros(1, tmpTheoryStruct{1}{j}.length);
    dotbar(ceil(thisDots)) = sets.bitmasking.prestretchPixelWidth_nm/sets.bitmasking.prestretchPixelWidth_nm;
    dotbar = apply_point_spread_function( ...
      dotbar, ...
      sets.bitmasking.psfSigmaWidth_nm/sets.bitmasking.prestretchPixelWidth_nm, ...
      1);
    [folder, filename, ext] = fileparts(tmpTheoryStruct{2}{j}.filename);
    newPath = fullfile(folder, ['dotbar_' filename ext]);
    fd = fopen(newPath,'w');
    fprintf(fd, strcat([' %5.' num2str(sets.theory.precision) 'f ']), dotbar);
    fclose(fd);
    tmpTheoryStruct{2}{j}.filename = newPath;
  end
end
% import CBT.Hca.Core.Analysis.convert_nm_ratio;
% tmpTheoryStruct{1} = convert_nm_ratio(sets.theory.nmbp, tmpTheoryStruct{1}, sets);
% tmpTheoryStruct{2} = convert_nm_ratio(sets.theory.nmbp, tmpTheoryStruct{2}, sets);
for j=1:length(tmpTheoryStruct{1})
  theoryStruct{j} = tmpTheoryStruct{1}{j};
  theoryStruct{j}.filename2 = tmpTheoryStruct{2}{j}.filename;
end