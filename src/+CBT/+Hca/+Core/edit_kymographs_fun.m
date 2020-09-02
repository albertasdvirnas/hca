function [ kymoStructs ] = edit_kymographs_fun(kymoStructs, timeFramesNr )
    % edit_kymographs_fun
    % Function for putting the kymographs into a session file
    %
    %     Args:
    %         kymoStructs: Kymograph structure
    %         timeFramesNr: number of timeframes
    % 
    %     Returns:
    %         kymoStructs: kymoStructs
    %

    % If the setting for number of timeframes is non-zero, we have to
    % remove some rows from the kymographs, and remove the kymographs with
    % insuficient amount of rows altogether
    if timeFramesNr ~= 0
        % number of timegrames
        timeframeTotal = cellfun(@(x) size(x.unalignedKymo, 1), kymoStructs);
        % remove kymo's that have not enough timeframes
        kymoStructs(timeframeTotal < timeFramesNr) = [];
        % dispplay which kymographs were removed
        disp(strcat(['Kymographs nr ' num2str(find(timeframeTotal < timeFramesNr)) ' were removed because they do not have enough time-frames'])); 
        % save only the preselected number of rows
        for i=1:length(kymoStructs)
%             try
                kymoStructs{i}.unalignedKymo = kymoStructs{i}.unalignedKymo(1:timeFramesNr,:);
%             catch
%                 kymoStructs{i} = [];
%             end
        end
    end    
    
end

