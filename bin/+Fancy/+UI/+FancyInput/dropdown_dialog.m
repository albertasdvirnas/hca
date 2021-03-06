function [choice, choiceIdx] = dropdown_dialog(dlgTitleStr, txtPromptMsg, choiceStrs, defaultChoiceIdx)
    if nargin < 4
        defaultChoiceIdx = 1;
    end

    hDialog = dialog(...
        'Units', 'normalized', ...
        'OuterPosition',[0.2 0.2 0.6 0.6], ...
        'Name', dlgTitleStr);
    [~] = uicontrol(...
        'Parent',hDialog,...
        'Style','text',...
        'Units', 'normalized', ...
        'OuterPosition',[0.2 0.6 0.6 0.1], ...
        'FontSize', 20, ...
        'String', txtPromptMsg);
       
    [~] = uicontrol(...
        'Parent',hDialog,...
        'Style','popup',...
        'Units', 'normalized', ...
        'OuterPosition',[0.2 0.4 0.6 0.1], ...
        'String', choiceStrs,...
        'Callback',@popup_callback);
       
    [~] = uicontrol(...
        'Parent',hDialog,...
        'Units', 'normalized', ...
        'OuterPosition',[0.2 0.2 0.6 0.1], ...
        'String','Continue',...
        'Callback','delete(gcf)');
       
    choiceIdx = defaultChoiceIdx;
    choice = choiceStrs{choiceIdx};
       
    % Wait for d to close before running to completion
    uiwait(hDialog);
   
    function popup_callback(popup, ~)
        idx = popup.Value;
        popup_items = popup.String;
        % This code uses dot notation to get properties.
        % Dot notation runs in R2014b and later.
        % For R2014a and earlier:
        % idx = get(popup,'Value');
        % popup_items = get(popup,'String');
        choiceIdx = idx;
        choice = char(popup_items(idx,:));
    end
end