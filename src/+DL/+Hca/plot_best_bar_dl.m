function [] = plot_best_bar_dl(fig1, barcodeGen, comparisonStruct, theoryStruct, userDefinedSeqCushion, ii, scoreType, typeInd)
% plot_best_bar

% plots best barcode vs theory in case barcode is always larger than
% theory
if nargin < 7
  userDefinedSeqCushion = 0;
end

len1=size(barcodeGen, 2);

thisCompStruct = subsref(comparisonStruct{ii}, substruct('.', scoreType));

or = thisCompStruct.or(1);
indcoef = thisCompStruct.maxcoefPartsCC(typeInd,1);
idx = thisCompStruct.idx;

switch typeInd
  case 1
    btype = 'dense';
    fname = theoryStruct{idx}.filename;
    try
      pos = thisCompStruct.optPos(1);
      stretch = thisCompStruct.optStr(1);
    catch
      pos = thisCompStruct.pos(1);
      stretch = thisCompStruct.bestBarStretch;
    end
  case 2
    btype = 'sparse';
    fname = theoryStruct{idx}.filename2;
    pos = thisCompStruct.pos(1);
    stretch = thisCompStruct.bestBarStretch;
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
v = linspace(1,expLen,expLen*stretch);
expBar = interp1(expBar, v);
expBit = logical(expBit(round(v)));
expBar(~expBit)= nan;


% expBar with expanded cushion
expBar = [nan(1, userDefinedSeqCushion) expBar nan(1, userDefinedSeqCushion)];
expBit = [false(1, userDefinedSeqCushion) expBit false(1, userDefinedSeqCushion)];

%     numSt = min(userDefinedSeqCushion,userDefinedSeqCushion-comparisonStruct{ii}.pos(1));
% don't allow theory start to loop over
theoryStart = pos-userDefinedSeqCushion;
theoryEnd = pos-userDefinedSeqCushion+length(expBar)-1;    

if theoryStart < 1
  % in circular case, these should be taken from the end of theorBar
  theorBar = [nan(abs(theoryStart)+1,1); theorBar];
  theoryEnd = theoryEnd + abs(theoryStart)+1;
  theoryStart = 1;
  thrLen = thrLen+abs(theoryStart)+1;
end

if theoryEnd > thrLen % now split this into linear and nonlinear case..
  theorBar = [ theorBar; theorBar(1:theoryEnd-thrLen)];
end

theoryIndBp = ceil((theoryStart:theoryEnd)*theoryStruct{idx}.pixelWidth_nm/theoryStruct{idx}.meanBpExt_nm);

barfragq = expBar;
barfragr = theorBar(theoryStart:theoryEnd);

switch typeInd
  case 1
    barfragq = (barfragq - nanmean(barfragq(expBit)))/nanstd(barfragq(expBit));
    barfragr = (barfragr - nanmean(barfragr(expBit)))/nanstd(barfragr(expBit));
  case 2
    barfragq = barfragq/nanmean(barfragq)*nanmean(barfragr);
end

if or == 2
  barfragq = fliplr(barfragq);
end

%     figure,
plot(theoryIndBp, barfragq)
hold on
plot(theoryIndBp, barfragr)

xlim([theoryIndBp(1) theoryIndBp(end)])


xlabel('Position along the sequence cushion (bp)','Interpreter','latex')
ylabel('Z-scored','Interpreter','latex')
if ii <= len1
  name = barcodeGen{ii}.name(1:end-4);
else
  name = 'consensus';
end

title(strcat(['Experimental barcode vs theory ']),'Interpreter','latex');
%
legend({strcat(['$\hat C_{\rm ' name '}^{' btype '}=$' num2str(indcoef,'%0.2f')]), niceName},'Interpreter','latex')



end

