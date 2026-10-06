function clim = autoCaxis(data, prcLow, prcHigh)
% autoCaxis  Robust percentile-based color scaling
% Usage:
%   clim = autoCaxis(data, 2, 98)

    v = data(:);
    v = v(isfinite(v));

    if isempty(v)
        clim = [-1 1];
        return;
    end

    clim = prctile(v, [prcLow prcHigh]);

    % Safety in case map is flat
    if clim(1) == clim(2)
        clim = clim + [-eps eps];
    end
end