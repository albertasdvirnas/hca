function [] = plot_best_bar_dl(fig1, barcodeGen, comparisonStruct, theoryStruct, userDefinedSeqCushion, ii, scoreType, typeInd)
% plot_best_bar

% plots best barcode vs theory in case barcode is always larger than
% theory
if nargin < 7
  userDefinedSeqCushion = 0;
end

len1=size(barcodeGen, 2);

thisCompStruct = subsref(comparisonStruct{ii}, substruct('.', scoreType));

pos = thisCompStruct.pos(1);
or = thisCompStruct.or(1);
stretch = thisCompStruct.bestBarStretch;
indcoef = thisCompStruct.indcoef(typeInd,1);
idx = thisCompStruct.idx;

switch typeInd
  case 1
    btype = 'cb';
    fname = theoryStruct{idx}.filename;
  case 2
    btype = 'dots';
    fname = theoryStruct{idx}.filename2;
end
tname = theoryStruct{idx}.name;
% theory length
thrLen = theoryStruct{idx}.length;

% load theory file
fileID = fopen(fname,'r');
formatSpec = '%f';
theorBar = fscanf(fileID,formatSpec);
fclose(fileID);

niceName = tname;
pl = [strfind(niceName,'NC') strfind(niceName,'NZ')];
niceName = niceName(pl:end);
pl = [strfind(niceName,'|') strfind(niceName,' ')];
niceName = strrep(niceName(1:(min(pl)-1)),'_','\_');
if isempty(niceName)
  niceName = strrep(tname,'_','\_');
end



% bitmask. In case of linear barcode, would like to modify this
%     theorBit = ones(1,thrLen);

% load either theory barcode or the consensus barcode
% try
  expBar = barcodeGen{1,ii}.rawBarcode;
  expBit = barcodeGen{1,ii}.rawBitmask;
% catch
%   try
%     expBar = consensusStruct.rawBarcode;
%     expBit = consensusStruct.rawBitmask;
%   catch
%     expBar = consensusStruct{ii-length(barcodeGen)}.rawBarcode;
%     expBit = consensusStruct{ii-length(barcodeGen)}.rawBitmask;
%   end
%   
% end

expLen = length(expBar);

% interpolate to the length which gave best CC value
expBar = interp1(expBar, linspace(1,expLen,expLen*stretch));
expBit = expBit(round(linspace(1,expLen,expLen*stretch)));
expBar(~expBit)= nan;


% expBar with expanded cushion
expBar = [ repmat(nan,1,userDefinedSeqCushion) expBar repmat(nan,1,userDefinedSeqCushion)];

%     numSt = min(userDefinedSeqCushion,userDefinedSeqCushion-comparisonStruct{ii}.pos(1));
% don't allow theory start to loop over
theoryStart = pos-userDefinedSeqCushion;
theoryEnd = pos-userDefinedSeqCushion+length(expBar)-1;

if theoryStart < 1
  % in circular case, these should be taken from the end of theorBar
  theorBar = [ repmat(nan,abs(theoryStart)+1,1); theorBar];
  theoryEnd = theoryEnd + abs(theoryStart)+1;
  theoryStart = 1;
  thrLen = thrLen+abs(theoryStart)+1;
end

if theoryEnd > thrLen % now split this into linear and nonlinear case..
  theorBar = [ theorBar; theorBar(1:theoryEnd-thrLen)];
end

theoryIndBp = ceil((theoryStart:theoryEnd)*theoryStruct{idx}.pixelWidth_nm/theoryStruct{idx}.meanBpExt_nm);

barStruct.bar1 = expBar;
barStruct.bar2= theorBar;

if or == 1
  barStruct.matchTable = [1 length(expBar) theoryStart theoryEnd or];
else
  barStruct.matchTable = [1 length(expBar) theoryEnd theoryStart or];
end

import CBT.Hca.UI.Helper.create_full_table;
[temp_table,barfragq, barfragr] = create_full_table(barStruct.matchTable, barStruct.bar1,barStruct.bar2,1);

switch typeInd
  case 1
    barfragq{1}(~isnan(barfragq{1})) = zscore(barfragq{1}(~isnan(barfragq{1})));
    barfragr{1}(~isnan(barfragr{1})) = zscore(barfragr{1}(~isnan(barfragr{1})));
  case 2
    barfragq{1}(~isnan(barfragq{1})) = barfragq{1}(~isnan(barfragq{1}))/mean(barfragq{1}(~isnan(barfragq{1})))*mean(barfragr{1}(~isnan(barfragr{1})));
end

%     figure,
plot(theoryIndBp, barfragq{1})
hold on
plot(theoryIndBp, barfragr{1})



xlabel('Position along the sequence cushion (bp)','Interpreter','latex')
ylabel('Z-scored','Interpreter','latex')
if ii <= len1
  name = num2str(ii);
else
  name = 'consensus';
end

title(strcat(['Experimental barcode vs theory ']),'Interpreter','latex');
%
legend({strcat(['$\hat Z_{\rm ' name '}^{' btype '}=$' num2str(indcoef,'%0.2f')]), niceName},'Interpreter','latex')



end

