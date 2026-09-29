function pretty = prettifyFilterDesc(desc)
%PRETTIFYFILTERDESC  Turn a loadGARuns filter tag into a readable heading.
%
%   prettifyFilterDesc('7C_CF3_TT1_GM1_sp100cm')
%     -> '7 cameras  .  CF3 (Combined)  .  UAV  .  Uniform grid  .  1.00 m spacing'
%
%   loadGARuns returns a compact tag naming the instance a figure was
%   built from. Every per-instance figure needs that instance spelled out
%   in its heading, otherwise the exported files are indistinguishable
%   from one another. Shared so convergence, diversity and anything added
%   later all phrase the instance identically.
%
%   Unrecognised tokens pass through unchanged.
%
%   See also loadGARuns, thesisTitle.
    tokens = strsplit(desc, '_');
    parts  = strings(1, 0);
    cfNames = {'Resolution Uncertainty', 'Dynamic Occlusion', 'Combined'};
    ttNames = {'UAV', 'UGV'};
    gmNames = {'Uniform grid', 'Normal grid'};
    for k = 1:numel(tokens)
        t = tokens{k};
        if endsWith(t, 'C') && all(isstrprop(t(1:end-1), 'digit'))
            n = str2double(t(1:end-1));
            parts(end+1) = sprintf('%d cameras', n);                         %#ok<AGROW>
        elseif startsWith(t, 'CF')
            n = str2double(t(3:end));
            if n>=1 && n<=numel(cfNames)
                parts(end+1) = sprintf('CF%d (%s)', n, cfNames{n});          %#ok<AGROW>
            else
                parts(end+1) = t;                                           %#ok<AGROW>
            end
        elseif startsWith(t, 'TT')
            n = str2double(t(3:end));
            if n>=1 && n<=numel(ttNames)
                parts(end+1) = ttNames{n};                                  %#ok<AGROW>
            else
                parts(end+1) = t;                                           %#ok<AGROW>
            end
        elseif startsWith(t, 'GM')
            n = str2double(t(3:end));
            if n>=1 && n<=numel(gmNames)
                parts(end+1) = gmNames{n};                                  %#ok<AGROW>
            else
                parts(end+1) = t;                                           %#ok<AGROW>
            end
        elseif startsWith(t, 'sp') && endsWith(t, 'cm')
            num = str2double(t(3:end-2));
            parts(end+1) = sprintf('%.2f m spacing', num/100);              %#ok<AGROW>
        else
            parts(end+1) = string(t);                                       %#ok<AGROW>
        end
    end
    pretty = char(strjoin(parts, '  ·  '));
end
