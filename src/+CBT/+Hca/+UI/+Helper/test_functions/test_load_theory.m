function [] = test_load_theory()
    % tests to show that theory loading works

    sets.theories = '/home/albyback/git/hca/tests/theoryData/theories_mat_2020-04-10_23_20_18_.txt';
    sets.theories = '/home/albyback/git/hca/tests/theoryData/theories_2020-04-10_23_42_14_.txt';
    sets.theory.precision = 5;
    sets.output.matDirpath = 'output/';
    newNmBp = 0.2;
% 
%     fd = fopen(sets.theories,'w');
%     fprintf(fd,'%s \n', '/home/albyback/git/hca/tests/data/theoryGen_2019-01-15_14_33_24_session.mat');
%     fclose(fd);
    fid = fopen(sets.theories); 
    fastaNames = textscan(fid,'%s','delimiter','\n'); fclose(fid);
    for i=1:length(fastaNames{1})
        [FILEPATH,NAME,EXT] = fileparts(fastaNames{1}{i});
        sets.theoryFile{i} = strcat(NAME,EXT);
        sets.theoryFileFold{i} = FILEPATH;
    end

           
    % case 1: 1 .mat file
    import CBT.Hca.UI.Helper.load_theory;
    [ theoryStruct ] = load_theory( sets );

    % convert
    import CBT.Hca.Core.Analysis.convert_nm_ratio;
    [ theoryStruct ] = convert_nm_ratio( newNmBp, theoryStruct,sets)
end

