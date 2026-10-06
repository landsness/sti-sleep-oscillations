function ROI = DrawROI_arn_v2(data, clim)
% DrawROI_arn_v2 — draw polygon ROI on a jet map with confirm/redraw
% IMPORTANT: does NOT close all figures (so it won't kill your main fig)

ROI = [];

while true
    % ----- Display map for drawing -----
    figDraw = figure('Name', 'Draw ROI', 'Color', 'w');
    imagesc(data);
    axis image off;
    colormap(jet); colorbar;

    if nargin >= 2 && ~isempty(clim)
        caxis(clim);
    end

    title('Click vertices; double-click to close polygon', ...
          'FontSize', 12, 'FontWeight', 'bold');

    % ----- Draw -----
    disp('👉 Click polygon vertices; double-click to close polygon');
    h = drawpolygon('LineWidth', 1.5, 'Color', 'k');
    ROI = createMask(h);
    ROI = logical(ROI);

    % ----- Show quick metrics -----
    roiPix = nnz(ROI);
    rp = regionprops(ROI, 'Centroid');
    fprintf('ROI pixels = %d\n', roiPix);
    if ~isempty(rp)
        fprintf('Centroid = [%.2f, %.2f]\n', rp.Centroid(1), rp.Centroid(2));
    end

    % ----- Optional overlay figure (separate) -----
    figOverlay = figure('Name', 'ROI Overlay', 'Color', 'w');
    imagesc(data);
    axis image off;
    colormap(jet); colorbar;
    if nargin >= 2 && ~isempty(clim)
        caxis(clim);
    end
    hold on;
    contour(ROI, [1 1], 'k', 'LineWidth', 1.5);
    hold off;
    title('ROI overlay', 'FontSize', 12, 'FontWeight', 'bold');

    % ----- Confirm -----
    resp = input('Happy with ROI? (0 = yes, 1 = no/redraw): ');
    if resp == 0
        % Close only the figures we created here
        if isgraphics(figDraw, "figure"), close(figDraw); end
        if isgraphics(figOverlay, "figure"), close(figOverlay); end
        break;
    else
        % Close only these and loop again
        if isgraphics(figDraw, "figure"), close(figDraw); end
        if isgraphics(figOverlay, "figure"), close(figOverlay); end
        fprintf('🔁 Redraw ROI...\n');
        ROI = [];
        continue;
    end
end

% Safety
if isempty(ROI)
    warning('⚠️ No ROI drawn — returning blank mask.');
    ROI = false(size(data));
end
end