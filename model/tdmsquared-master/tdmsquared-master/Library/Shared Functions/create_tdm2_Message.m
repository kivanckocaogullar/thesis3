function message = create_tdm2_Message(msgType, timestamp, blockName, text)
    % Force all inputs to row char arrays
    blockName = char(blockName(:).');   % row vector
    text = char(text(:).');             % row vector
    timeStr = ['t = ', num2str(timestamp)];

    % Extract display name
    if contains(blockName, '/')
        index = strfind(blockName, '/'); 
        blockNameShort = blockName(1:(index(end)-1));
    else
        blockNameShort = blockName;
    end

    linkStr = ['<a href="matlab:hilite_system(''', blockNameShort, ''')">', blockNameShort, '</a>: '];

    % Choose message type
    switch msgType
        case 'info'
            prefix = 'tDMSquared Info @ ';
        case 'warning'
            prefix = 'tDMSquared Warning @ ';
        case 'error'
            prefix = 'tDMSquared Error @ ';
        otherwise
            prefix = 'tDMSquared Message @ ';
    end

    % Concatenate everything (all char row vectors)
    message = [prefix, timeStr, ' caused by ', linkStr, text];
end
