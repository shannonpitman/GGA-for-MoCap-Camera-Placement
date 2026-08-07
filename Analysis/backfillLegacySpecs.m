function specs = backfillLegacySpecs(specs)
% BACKFILLLEGACYSPECS  Add any fields that older saved `specs` structs are
% missing, so re-evaluation under the corrected cost functions works.
%
% Older runs (pre-FocalWide / pre-PreComputed.maxCameraRangeWide saves)
% only contain a subset of the hardware/cost parameters that the current
% cost functions rely on. Calling resUncertaintyCost / dynamicOcclusionCost
% on those structs throws "Unrecognized field name" errors. This helper
% fills in any missing fields with the current defaults from
% setupHardwareSpecs and setupCostParams, leaving anything already present
% untouched.

    if ~isfield(specs, 'Cams') || isempty(specs.Cams)
        error('backfillLegacySpecs:NoCams', ...
            'specs.Cams is missing — cannot infer the hardware defaults.');
    end

    % --- Hardware-level defaults from the current setupHardwareSpecs ---
    defaults = setupHardwareSpecs(specs.Cams);
    fns = fieldnames(defaults);
    for i = 1:numel(fns)
        if ~isfield(specs, fns{i}) || isempty(specs.(fns{i}))
            specs.(fns{i}) = defaults.(fns{i});
        end
    end

    % --- PreComputed cost-param defaults: only re-derive missing pieces ---
    if ~isfield(specs, 'PreComputed') || isempty(specs.PreComputed)
        specs = setupCostParams(specs);
        return;
    end

    % Build a reference PreComputed struct using setupCostParams against a
    % temporary copy, then merge in only the fields that are missing from
    % the loaded one.
    tmp = setupCostParams(specs);
    refPC = tmp.PreComputed;
    pcFns = fieldnames(refPC);
    for i = 1:numel(pcFns)
        if ~isfield(specs.PreComputed, pcFns{i}) || ...
                isempty(specs.PreComputed.(pcFns{i}))
            specs.PreComputed.(pcFns{i}) = refPC.(pcFns{i});
        end
    end

    % Make sure derived target fields exist for older saves.
    if ~isfield(specs, 'TargetHomogeneous') || isempty(specs.TargetHomogeneous)
        if isfield(specs, 'Target') && ~isempty(specs.Target) ...
                && isfield(specs, 'NumPoints')
            specs.TargetHomogeneous = [specs.Target'; ones(1, specs.NumPoints)];
        end
    end
end
