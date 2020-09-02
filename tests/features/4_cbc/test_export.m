load('cbc_test.mat');

% load some result of HCA which was saved before for testing

consensusStructs.consensusStruct{1} = consensusStruct;

sets.output = sets.output.matDirpath;

% if ~exists(sets.output)
mkdir(sets.output);
% end

import CBT.Hca.Export.export_cbc_compatible_consensus;
export_cbc_compatible_consensus(consensusStructs, barcodeGen,kymoStructs,sets);
