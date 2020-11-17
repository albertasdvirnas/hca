function [] = plot_best_bar_dl(fig1, barcodeGen, consensusStruct, comparisonStruct,theoryStruct ,userDefinedSeqCushion, ii, scoreType, typeInd)
% plot_best_bar

% plots best barcode vs theory in case barcode is always larger than
% theory
if nargin < 7
  userDefinedSeqCushion = 0;
end

len1=size(barcodeGen, 2);

switch scoreType
  case 'dual'
    pos = comparisonStruct{ii}.dual.pos(1);
    or = comparisonStruct{ii}.dual.or(1);
    stretch = comparisonStruct{ii}.dual.bestBarStretch;
    indcoef = cellfun(@(x) x.dual.indcoef,comparisonStruct,'UniformOutput',false)';
    indcoef = indcoef{ii}(typeInd,1);
  case 'dense'
    pos = comparisonStruct{ii}.dense.pos(1);
    or = comparisonStruct{ii}.dense.or(1);
    stretch = comparisonStruct{ii}.dense.bestBarStretch;
    indcoef = cellfun(@(x) x.dense.indcoef,comparisonStruct,'UniformOutput',false)';
    indcoef = indcoef{ii}(typeInd,1);
  case 'sparse'
    pos = comparisonStruct{ii}.sparse.pos(1);
    or = comparisonStruct{ii}.sparse.or(1);
    stretch = comparisonStruct{ii}.sparse.bestBarStretch;
    indcoef = cellfun(@(x) x.sparse.indcoef,comparisonStruct,'UniformOutput',false)';
    indcoef = indcoef{ii}(typeInd,1);
end

switch typeInd
  case 1
    btype = 'cb';
  case 2
    btype = 'dots';
end

% load theory file
fileID = fopen(theoryStruct{comparisonStruct{ii}.idx}.filename,'r');
formatSpec = '%f';
theorBar = fscanf(fileID,formatSpec);
fclose(fileID);

niceName = theoryStruct{comparisonStruct{ii}.idx}.name;
pl = [strfind(niceName,'NC') strfind(niceName,'NZ')];
niceName = niceName(pl:end);
pl = [strfind(niceName,'|') strfind(niceName,' ')];
niceName = strrep(niceName(1:(min(pl)-1)),'_','\_');
if isempty(niceName)
  niceName = strrep(theoryStruct{comparisonStruct{ii}.idx}.name,'_','\_');
end


% theory length
thrLen = theoryStruct{comparisonStruct{ii}.idx}.length;

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
plot(barfragq{1})
hold on
plot(barfragr{1})



xlabel('Position along the sequence cushion (px)','Interpreter','latex')
ylabel('Z-scored','Interpreter','latex')
if ii <= len1
  name = num2str(ii);
else
  name = 'consensus';
end

title(strcat(['Experimental barcode vs theory ']),'Interpreter','latex');
%
legend({strcat(['$\hat C_{\rm ' name '}^{' btype '}=$' num2str(indcoef,'%0.2f')]), niceName},'Interpreter','latex')



end

