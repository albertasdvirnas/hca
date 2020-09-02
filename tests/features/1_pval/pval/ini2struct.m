function outStruct = ini2struct( iniFilepath )
    % INI2STRUCT - Reads and parses an INI file to return a structure with
    %    section names and keys as fields.
    %
    % Inputs:
    %   iniFilepath
    %     the filepath for the INI
    %
    % Outputs
    %   outStruct
    %     a scalar struct with potentially nested scalar structs for field
    %     values where values are pulled from ini-file
    %
    % A slight modification of ini2struct.m:
    %   http://www.mathworks.com/matlabcentral/fileexchange/45725-ini2struct/content//ini2struct.m
    % Original BSD license here:
    %  http://www.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=45725
    % 
    % Authors:
    %   Saair Quaderi
    %     minor cleanup, moved from CMN_HelperFunctions
    %   Tobias Ambjörnsson (2015)

    outStruct = struct;
    inputFile = fopen(iniFilepath,'r'); % open file
    
    if inputFile == -1
        error(['Failed to open file: ''', iniFilepath, '''']);
    end
    
    sectionName = '';
    while ~feof(inputFile) % read until it ends
        currentLine = strtrim(fgetl(inputFile)); % remove leading/trailing spaces
        if isempty(currentLine) || currentLine(1)=='%' || currentLine(1)=='#' % skip empty & comments lines
            continue
        end
        if currentLine(1)=='[' % section header
            sectionName = strtok(currentLine(2:end), ']');
            sectionName = matlab.lang.makeValidName(sectionName);
            outStruct.(sectionName) = struct; % create field
            continue
        end
        
        [keyName, val] = strtok(currentLine, '='); % Expected line format: "Key = Value ; comment"
        keyName = matlab.lang.makeValidName(keyName);
        
        val = strtrim(val(2:end)); % remove spaces after "="
        
        if isempty(val) || val(1)=='%' || val(1)=='#' % empty entry
            val = [];
        elseif val(1)=='"' % double-quoted string
            val = strtok(val, '"');
        elseif val(1)=='''' % single-quoted string
            val = strtok(val, '''');
        else
            val = strtok(val, '%'); % remove inline comment
            val = strtok(val, '#'); % remove inline comment
            val = strtrim(val); % remove spaces before comment

            % convert string to number(s)
            [numVal, successfulConversion] = str2num(val); %#ok<ST2NM>
            if successfulConversion
                val = numVal;
            end
        end
        
        if isempty(sectionName)
            outStruct.(keyName) = val;
        else
            outStruct.(sectionName).(keyName) = val;
        end
    end
    fclose(inputFile);
end
        