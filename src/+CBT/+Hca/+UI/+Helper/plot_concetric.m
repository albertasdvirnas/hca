function plot_concetric(hAxis,barStruct,sets)
    % plot concentric
    
    if nargin < 3
         sets = struct;
         sets.A = 'experiment';
         sets.B = 'theory';
    end
    
    sets.BG_COLOR = [1.0, 1.0, 1.0];
    sets.COLORBAR_COLOR = [0 0 0 ];
    sets.COLOR_MAP = 1/255*[200,200,200];
    sets.INNER_BAND_WIDTH = 8; %initial radius
    sets.FINAL_GAP_WIDTH = 2; % between barcodes
    sets.BARCODE_WIDTH = 8; % barcodewidth
%     sets.NONCIRCADD = 20; % add pixels to the end to make sure it looks non-circular
%     numBarcodes = 2;

    % this table should be a bit different, in case comparison is linear.
    % When we do not allow loop-arounds..
    import CBT.Hca.UI.Helper.create_full_table;
    [temp_table,barfragq, barfragr] = create_full_table(barStruct.matchTable, barStruct.bar1,barStruct.bar2,1);

%     import Plot.shift_barcode_data;
%     [barStruct] = shift_barcode_data(barStruct);

    import Barcoding.Visualizing.colormap_kry;
    for i=1:length(barfragq)

        bar1 = barfragq{i};
        bar2 = barfragr{i}';

        lT = length(bar1);
        lN = length(bar2);

        % circle radius'es
        radii = cumsum([0 sets.INNER_BAND_WIDTH sets.BARCODE_WIDTH sets.FINAL_GAP_WIDTH sets.BARCODE_WIDTH ]');
        % maximum length  
        
        if sets.theory.isLinearTF == 1
            maxLen = 2*max(lT,lN);%+sets.NONCIRCADD;
        else
            maxLen = max(lT,lN);
        end

        % angles. The last point is the same as the first point. We start at
        % -pi
        theta = pi:-(2*pi/maxLen):-pi;
        % xcoord
        Xs = radii*cos(theta);
        % ycoord
        Ys = radii*sin(theta);

        % create Cs structure. Keep 1st and 3rd 5th rows empty
        Cs = NaN(2*2 + 1, maxLen+1);
        % if barcodes are linear, we add some extra nan's. Last point not
        % really necessary
        Cs(2,:) = [bar1 NaN(1, maxLen - length(bar1)) bar1(1)];
        try
            Cs(4,:) = [bar2 NaN(1, maxLen - length(bar2)) bar2(1)];
        catch
            Cs(4,:) = [bar2; NaN(maxLen - length(bar2),1); bar2(1)];
        end

        % create figure 
%         f = figure;
        % what if it is subfigure
%         hPanelUnaligned = uipanel('Parent', fig1);
%         hAxis = axes(...
%             'Units', 'normal', ...
%             'Position', [0, 0.1, 0.9, 0.9], ...
%             'Parent', hPanelUnaligned ...
%             );

        hGoSurface = pcolor(hAxis, Xs, Ys, Cs);
%         set(hAxis, 'color', sets.BG_COLOR, 'xtick', [], 'ytick', [], 'box', 'off');
        [~] = colorbar(hAxis, 'Color', sets.COLORBAR_COLOR, 'FontSize', 12);
        text(radii(3)*sin(6*pi/4),radii(3)*cos(6*pi/4),sets.A,'FontSize', 6,'Interpreter','latex');
        text(radii(5)*sin(6*pi/4),radii(5)*cos(6*pi/4),sets.B,'FontSize', 6,'Interpreter','latex');

%         set(hAxis.Parent, 'BackgroundColor', sets.BG_COLOR);
%         legend({'one','two'},'Interpreter','latex')
%         set(hAxis,'Position')
% set(ax,'YTickLabel',a)
    %     colormap gray;
        axis square;
        axis off;
        set(hAxis, 'XTickLabel', [], 'YTickLabel',[]); 
        grid off;
        shading flat; %shading interp;
    end
end

