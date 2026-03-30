function label = getModalityLabel(specs)
% target space modality
% e.g. 'UAV | Uniform Grid' or 'UGV | Normal Grid'

    % Target type
    if specs.TargetType == 1
        typeStr = 'UAV (Full Volume)';
    elseif specs.TargetType == 2
        typeStr = 'UGV (Floor Slab)';
    else
        typeStr = 'Unknown Target Type';
    end

    % Discretisation mode
    if specs.TargetMode == 1
        modeStr = 'Uniform Grid';
    elseif specs.TargetMode == 2
        modeStr = 'Normal Grid';
    else
        modeStr = 'Unknown Grid Mode';
    end

    label = sprintf('%s | %s', typeStr, modeStr);
end