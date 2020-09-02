function [theoryGen] = gen_theory(sequence,sets)
    %
    import CBT.Hca.Core.Theory.choose_cb_model;
    [sets.model ] = choose_cb_model(sets.theoryGen.method);

    import CBT.Hca.Core.Theory.compute_free_conc;
    sets = compute_free_conc(sets);

    [ theory, header] = gen_hca_theory( sequence,sets);
 
    theoryGen.theoryBarcodes{1} = theory;
    theoryGen.theoryNames{1}= header;
    theoryGen.theoryIdx{1} = 1;
    theoryGen.bpNm{1} = sets.theoryGen.meanBpExt_nm/sets.theoryGen.psfSigmaWidth_nm;
    theoryGen.sets = sets.theoryGen;
%     theoryGen.grayscaleVideoRescaled = grayscaleVideoRescaled;
%     theoryGen.molExtracted = molExtracted;
%     theoryGen.moldata = moldata;
    % save sets
    matFilename = strcat(['example_session.mat']);
    matFilepath = fullfile(sets.resultsDir, matFilename);
    save(matFilepath, 'theoryGen');

        
end

