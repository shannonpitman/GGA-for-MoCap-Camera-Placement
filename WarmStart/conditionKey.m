function key = conditionKey(numCams, costFunc, targetType, gridMode, spacing)
%CONDITIONKEY  Label for one experimental condition, e.g. '7C_CF3_TT1_GM1_sp1.00'.
%   Same format as the ConditionKey in batchRunGA schedules, so the seed
%   archive and the batch logs use one name for a condition.
    key = sprintf('%dC_CF%d_TT%d_GM%d_sp%.2f', numCams, costFunc, targetType, gridMode, spacing);
end
