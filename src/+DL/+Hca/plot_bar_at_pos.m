function [] = plot_bar_at_pos(barcodeGen, ii, theoryStruct, pos, or, stretch, idx, score, userDefinedSeqCushion, typeInd)
% plot_best_bar

% plots best barcode vs theory in case barcode is always larger than
% theory

switch typeInd
  case 1
    fname = theoryStruct{idx}.filename;
  case 2
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

expBar = barcodeGen{ii}.rawBarcode;
expBar = interp1(expBar, linspace(1,length(expBar),length(expBar)*stretch));
name = barcodeGen{ii}.name(1:end-4);
% expBar with expanded cushion
expBar = [nan(1,userDefinedSeqCushion) expBar nan(1,userDefinedSeqCushion)];

%     numSt = min(userDefinedSeqCushion,userDefinedSeqCushion-comparisonStruct{ii}.pos(1));
% don't allow theory start to loop over
theoryStart = pos-userDefinedSeqCushion;
theoryEnd = pos-userDefinedSeqCushion+length(expBar)-1;

if theoryStart < 1
  % in circular case, these should be taken from the end of theorBar
  theorBar = [ nan(abs(theoryStart)+1,1); theorBar];
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

title(strcat(['Experimental barcode vs theory ']),'Interpreter','latex');
%
legend({strcat(['$\hat C_{\rm{' name '}}^{\rm{bionano}}=$' num2str(score,'%0.2f')]), niceName},'Interpreter','latex')



end

