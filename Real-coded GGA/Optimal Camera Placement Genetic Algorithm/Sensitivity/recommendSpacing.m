function recs = recommendSpacing(sweep, tolerancePct, verbose)
%RECOMMENDSPACING  Coarsest spacing meeting a tolerance, per (config, CF).
%
%   recs = recommendSpacing(sweep, tolerancePct) returns a struct array
%   with the largest grid spacing whose absolute relative deviation from
%   the finest-grid (0.1 m) reference cost stays under tolerancePct (in
%   percent, e.g. 2 for 2%).
%
%   recs = recommendSpacing(sweep, tolerancePct, verbose) suppresses the
%   console table when verbose is false (default true).
%
%   Each element of recs has fields:
%       config             - configuration name
%       cf                 - cost-function tag (e.g. 'CF3_combined')
%       referenceSpacing   - finest spacing used as reference
%       referenceCost      - cost at that reference spacing
%       recommendedSpacing - coarsest spacing within tolerance
%       deviationPctAtRec  - absolute relative deviation [%] at that spacing
%       tolerancePct       - tolerance threshold used
%
%   Notes
%     * Visual inspection of the cost-vs-spacing curves remains the
%       primary check; this helper just makes the threshold explicit so a
%       written justification can cite a number rather than "looks flat".
%     * If even the second spacing already exceeds tolerance, the
%       reference (finest) spacing is recommended and a warning is logged.
%     * The recommendation is across both raw CFs (CF1, CF2) and the
%       combined CF3. Use the strictest (smallest) recommended spacing
%       across all three for a conservative pick.

    if nargin < 3, verbose = true; end

    results = sweep.results;
    cfNames = sweep.cfNames;
    configNames = {sweep.configs.name};

    recs = repmat(struct('config','', 'cf','', ...
                         'referenceSpacing',NaN, 'referenceCost',NaN, ...
                         'recommendedSpacing',NaN, 'deviationPctAtRec',NaN, ...
                         'tolerancePct',tolerancePct), ...
                  numel(configNames)*numel(cfNames), 1);
    k = 0;

    for c = 1:numel(configNames)
        for f = 1:numel(cfNames)
            sub = filterRows(results, configNames{c}, cfNames{f});
            ref = sub.cost(1);
            if ref == 0
                relPct = abs(sub.cost - ref) * 100;     % avoid div-by-zero
            else
                relPct = abs((sub.cost - ref) / ref) * 100;
            end

            within = relPct <= tolerancePct;
            % Coarsest spacing where this AND all finer spacings stay within
            % tolerance — guards against a false-positive crossing back inside.
            cumWithin = cummin(double(within));
            okIdx = find(cumWithin > 0, 1, 'last');
            if isempty(okIdx) || okIdx == 1
                recommended = sub.spacing(1);
                devAtRec    = relPct(1);
                if verbose && (isempty(okIdx) || okIdx == 1) && relPct(2) > tolerancePct
                    warning('recommendSpacing:NoRoom', ...
                        ['No spacing coarser than the %.3g m reference stays ' ...
                         'within %.2f%% for %s / %s. Falling back to ref.'], ...
                        sub.spacing(1), tolerancePct, configNames{c}, cfNames{f});
                end
            else
                recommended = sub.spacing(okIdx);
                devAtRec    = relPct(okIdx);
            end

            k = k + 1;
            recs(k).config             = configNames{c};
            recs(k).cf                 = cfNames{f};
            recs(k).referenceSpacing   = sub.spacing(1);
            recs(k).referenceCost      = ref;
            recs(k).recommendedSpacing = recommended;
            recs(k).deviationPctAtRec  = devAtRec;
            recs(k).tolerancePct       = tolerancePct;
        end
    end

    if verbose
        printRecommendations(recs);
    end
end


function printRecommendations(recs)
    fprintf('\n=================================================================\n');
    fprintf('  Recommended grid spacing (tolerance = %.2f%%)\n', recs(1).tolerancePct);
    fprintf('=================================================================\n');
    fprintf('  %-32s  %-15s  %12s  %10s\n', 'Config', 'Cost function', 'Recommended', 'Dev [%]');
    fprintf('  %s\n', repmat('-', 1, 75));
    for i = 1:numel(recs)
        fprintf('  %-32s  %-15s  %10.3g m  %9.3f\n', ...
            truncate(recs(i).config, 32), recs(i).cf, ...
            recs(i).recommendedSpacing, recs(i).deviationPctAtRec);
    end
    fprintf('-----------------------------------------------------------------\n');
    % Most conservative spacing per config across all CFs
    cfgs = unique({recs.config}, 'stable');
    for c = 1:numel(cfgs)
        sel = strcmp({recs.config}, cfgs{c});
        worst = min([recs(sel).recommendedSpacing]);
        fprintf('  Conservative pick for %-25s : %.3g m\n', truncate(cfgs{c},25), worst);
    end
    fprintf('=================================================================\n\n');
end


function s = truncate(s, n)
    if numel(s) > n, s = [s(1:n-1) '~']; end
end


function sub = filterRows(results, configName, cfn)
    rows = results( strcmp({results.config}, configName) & ...
                    strcmp({results.cf},     cfn) );
    [sp, ord] = sort([rows.spacing]);
    sub.spacing   = sp(:);
    sub.cost      = reshape([rows(ord).cost], [], 1);
end
