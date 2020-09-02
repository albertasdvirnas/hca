function plot_concetric_with_root(hAxis,alignedBar,bit,posVec,sets)
    % plot concentric
    
    if nargin < 5
         sets = struct;
         sets.A = 'experiment';
         sets.B = 'theory';
        sets.BG_COLOR = [1.0, 1.0, 1.0];
        sets.COLORBAR_COLOR = [0 0 0 ];
        sets.COLOR_MAP = 1/255*[200,200,200];
        sets.INNER_BAND_WIDTH = 8; %initial radius
        sets.FINAL_GAP_WIDTH = 2; % between barcodes
        sets.BARCODE_WIDTH = 8; % barcodewidth
    end
    
 
%     sets.NONCIRCADD = 20; % add pixels to the end to make sure it looks non-circular
%     numBarcodes = 2;

    % this table should be a bit different, in case comparison is linear.
    % When we do not allow loop-arounds..
%     import CBT.Hca.UI.Helper.create_full_table;
%     [temp_table,barfragq, barfragr] = create_full_table(barStruct.matchTable, barStruct.bar1,barStruct.bar2,1);

%     import Plot.shift_barcode_data;
%     [barStruct] = shift_barcode_data(barStruct);

    import Barcoding.Visualizing.colormap_kry;

    % barcodes and bitmasks to be plotted
    bars = alignedBar(posVec(:,1));
    bits = bit(posVec(:,1));
     
    % barcode lengths
    lT = cellfun(@(x) length(x),bars);

    % circle radius'es - depend on how many barcodes there are
%     radii = cumsum([0 sets.INNER_BAND_WIDTH sets.BARCODE_WIDTH sets.FINAL_GAP_WIDTH sets.BARCODE_WIDTH ]');
        % maximum length  
    radii =   cumsum([0 sets.INNER_BAND_WIDTH sets.BARCODE_WIDTH repmat([sets.FINAL_GAP_WIDTH sets.BARCODE_WIDTH],1,length(bars)-1)]');
  
%     repmat([sets.FINAL_GAP_WIDTH sets.BARCODE_WIDTH],1,length(bars)-1)
%         if sets.theory.isLinearTF == 1
%             maxLen = 2*max(lT,lN);%+sets.NONCIRCADD;
%         else
%             maxLen = max(lT,lN);
%         end
    
    % since we don't know the size, this is determined by posVec later, now
    % a quick hack to make it work
    maxLen = 3000;

    % angles. The last point is the same as the first point. We start at
    % -pi
    theta = pi:-(2*pi/maxLen):-pi;
    % xcoord
    Xs = radii*cos(theta);
    % ycoord
    Ys = radii*sin(theta);

    % create Cs structure. Keep 1st and 3rd 5th rows empty
    Cs = NaN(2*length(bars) + 1, maxLen+1);
    % if barcodes are linear, we add some extra nan's. Last point not
    % really necessary
    for i=1:length(bars)
        bars{i}(~bits{i}) = nan;
%         Cs(2*i,:) = [bars{i}; NaN(maxLen - length(bars{i}),1); bars{i}(1)];
        bars{i} = (bars{i}-nanmean(bars{i}))/nanstd(bars{i});
        if ~posVec(i,3)
             bars{i} = flipud( bars{i});
        end

        if posVec(i,2) > 0
        	Cs(2*i,posVec(i,2):posVec(i,2)+lT(i)-1) = bars{i};
        else
            Cs(2*i,end+posVec(i,2):end) = bars{i}(1:-posVec(i,2)+1);
            Cs(2*i,1:lT(i)+posVec(i,2)-1) = bars{i}(-posVec(i,2)+2:lT(i));
        end

    end

    hGoSurface = pcolor(hAxis, Xs, Ys, Cs);
    % set(hAxis, 'color', sets.BG_COLOR, 'xtick', [], 'ytick', [], 'box', 'off');
    [~] = colorbar(hAxis, 'Color', sets.COLORBAR_COLOR, 'FontSize', 12);

    % now can add also colorbar 
%         text(radii(3)*sin(6*pi/4),radii(3)*cos(6*pi/4),sets.A,'FontSize', 6,'Interpreter','latex');
%         text(radii(5)*sin(6*pi/4),radii(5)*cos(6*pi/4),sets.B,'FontSize', 6,'Interpreter','latex');

%         set(hAxis.Parent, 'BackgroundColor', sets.BG_COLOR);
%         legend({'one','two'},'Interpreter','latex')
%         set(hAxis,'Position')
% set(ax,'YTickLabel',a)
%         colormap gray;
        axis square;
        axis off;
        set(hAxis, 'XTickLabel', [], 'YTickLabel',[]); 
        grid off;
        shading flat; %shading interp;
end

