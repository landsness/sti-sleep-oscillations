%% ===============================================================
%  SO_ROI_Quantification_Mirrored_v2.m
%  Author: Arnav Ajay Jadav | Landsness Lab, WashU
%
%  Purpose:
%   - Load raw SO maps (bsl/acute/wk1) from Mouse_##_SO_RawMaps.mat
%   - User draws ROI1 and ROI2 ON ONE chosen map (gold standard)
%   - Mirror ROI across midline (contra)
%   - Apply same ROIs to all maps to compute RAW ipsi/contra means
%   - Compute ratios at end (acute/bsl, wk1/bsl) for convenience
%   - Preview values before append
%   - Append one row per mouse into ROI_Summary_V2.xlsx
%   - Verify write by reading back the written row
% ===============================================================

clear; clc; close all;

%% ---------------- USER INPUTS ----------------
baseDir   = "Z:\Results\";
excelPath = fullfile(baseDir, "ROI_Summary_V2.xlsx");

mouseNum = input('Enter mouse number (e.g., 11, 15, 24): ', 's');
mouseNumN = str2double(mouseNum);


pattern = sprintf("Mouse_%s_SO_RawMaps.mat", mouseNum);
fprintf("🔍 Searching for %s under %s ...\n", pattern, baseDir);

fileList = dir(fullfile(baseDir, "**", pattern));
if isempty(fileList)
    error("❌ No SO_RawMaps.mat found for mouse %s", mouseNum);
end

% choose newest
[~, idx] = max([fileList.datenum]);
rawMapsPath = fullfile(fileList(idx).folder, fileList(idx).name);
fprintf("✅ Using %s\n", rawMapsPath);

S = load(rawMapsPath);
req = ["so_map_bsl","so_map_acute","so_map_wk1"];
for r = req
    if ~isfield(S, r)
        error("❌ Raw maps file missing variable: %s", r);
    end
end

bsl   = S.so_map_bsl;
acute = S.so_map_acute;
wk1   = S.so_map_wk1;

% optional recovery map (only for drawing choice)
recovery = wk1 - acute;

%% ---------------- CHOOSE DRAW MAP ----------------
fprintf("\nChoose which map to draw ROIs on FIRST (gold standard):\n");
fprintf("  1 = Acute\n");
fprintf("  2 = Baseline (BSL)\n");
fprintf("  3 = Week 1 (WK1)\n");
fprintf("  4 = Recovery (WK1 - Acute)\n");
drawChoice = input("Enter 1/2/3/4: ");

switch drawChoice
    case 1
        drawMap = acute; drawLabel = "ACUTE";
    case 2
        drawMap = bsl;   drawLabel = "BSL";
    case 3
        drawMap = wk1;   drawLabel = "WK1";
    case 4
        drawMap = recovery; drawLabel = "RECOVERY";
    otherwise
        error("Invalid choice.");
end

%% ---------------- DISPLAY DRAW MAP ----------------
fig = figure('Name', sprintf('Mouse %s — Draw ROIs on %s', mouseNum, drawLabel), ...
             'Color', 'w');
imagesc(drawMap);
axis image off;
colormap(jet); colorbar;

% robust scaling (needs your autoCaxis.m function on path)
clim = autoCaxis(drawMap, 2, 98);
caxis(clim);

title(sprintf('Mouse %s — DRAW ON %s (ROI1 then ROI2)', mouseNum, drawLabel), ...
      'FontSize', 13, 'FontWeight', 'bold');
drawnow;

%% ---------------- ROI DRAW + METRICS ----------------
SO_values = struct();

roiMean = @(mask, data) mean(data(mask), 'omitnan');

for roiIdx = 1:2
    fprintf('\n==============================\n');
    fprintf('✏️  Draw ROI %d on %s now...\n', roiIdx, drawLabel);

    % Draw on chosen map
    roi = DrawROI_arn_v2(drawMap, clim);  % logical mask
    if ~islogical(roi), roi = logical(roi); end

    % Mirror across midline
    mirrored_roi = fliplr(roi);

    % Overlay verification
    figure(fig); hold on;
    contour(roi, [1 1], 'y', 'LineWidth', 1.4);
    contour(mirrored_roi, [1 1], 'w--', 'LineWidth', 1.2);
    hold off; drawnow;

    % Geometry
    roi_area    = nnz(roi);
    contra_area = nnz(mirrored_roi);

    c1 = regionprops(roi,'Centroid');
    c2 = regionprops(mirrored_roi,'Centroid');

    if isempty(c1), centroid = [NaN NaN]; else, centroid = c1.Centroid; end
    if isempty(c2), centroid_contra = [NaN NaN]; else, centroid_contra = c2.Centroid; end

    SO_values.(sprintf("ROI%d_Pixels", roiIdx))           = roi_area;
    SO_values.(sprintf("ROI%d_ContraPixels", roiIdx))     = contra_area;
    SO_values.(sprintf("ROI%d_CentroidX", roiIdx))        = centroid(1);
    SO_values.(sprintf("ROI%d_CentroidY", roiIdx))        = centroid(2);
    SO_values.(sprintf("ROI%d_ContraCentroidX", roiIdx))  = centroid_contra(1);
    SO_values.(sprintf("ROI%d_ContraCentroidY", roiIdx))  = centroid_contra(2);

    % RAW means across maps
    SO_values.(sprintf("ROI%d_bsl_Ipsi", roiIdx))     = roiMean(roi, bsl);
    SO_values.(sprintf("ROI%d_bsl_Contra", roiIdx))   = roiMean(mirrored_roi, bsl);

    SO_values.(sprintf("ROI%d_acute_Ipsi", roiIdx))   = roiMean(roi, acute);
    SO_values.(sprintf("ROI%d_acute_Contra", roiIdx)) = roiMean(mirrored_roi, acute);

    SO_values.(sprintf("ROI%d_wk1_Ipsi", roiIdx))     = roiMean(roi, wk1);
    SO_values.(sprintf("ROI%d_wk1_Contra", roiIdx))   = roiMean(mirrored_roi, wk1);

    % Convenience ratios (Eric can recompute)
    SO_values.(sprintf("ROI%d_ACUToverBSL_Ipsi", roiIdx)) = ...
        SO_values.(sprintf("ROI%d_acute_Ipsi", roiIdx)) ./ SO_values.(sprintf("ROI%d_bsl_Ipsi", roiIdx));
    SO_values.(sprintf("ROI%d_ACUToverBSL_Contra", roiIdx)) = ...
        SO_values.(sprintf("ROI%d_acute_Contra", roiIdx)) ./ SO_values.(sprintf("ROI%d_bsl_Contra", roiIdx));

    SO_values.(sprintf("ROI%d_WK1overBSL_Ipsi", roiIdx)) = ...
        SO_values.(sprintf("ROI%d_wk1_Ipsi", roiIdx)) ./ SO_values.(sprintf("ROI%d_bsl_Ipsi", roiIdx));
    SO_values.(sprintf("ROI%d_WK1overBSL_Contra", roiIdx)) = ...
        SO_values.(sprintf("ROI%d_wk1_Contra", roiIdx)) ./ SO_values.(sprintf("ROI%d_bsl_Contra", roiIdx));

    % Debug line (quick sanity)
    fprintf("ROI%d nnz=%d | acute ipsi=%.4f | bsl ipsi=%.4f\n", roiIdx, nnz(roi), ...
        SO_values.(sprintf("ROI%d_acute_Ipsi", roiIdx)), SO_values.(sprintf("ROI%d_bsl_Ipsi", roiIdx)));
end

%% ---------------- HEADERS + FIELD LIST ----------------
headers = {'Mouse', ...
    'ROI1_bsl_Ipsi','ROI1_bsl_Contra','ROI1_acute_Ipsi','ROI1_acute_Contra','ROI1_wk1_Ipsi','ROI1_wk1_Contra', ...
    'ROI1_Pixels','ROI1_ContraPixels','ROI1_CentroidX','ROI1_CentroidY','ROI1_ContraCentroidX','ROI1_ContraCentroidY', ...
    'ROI1_ACUToverBSL_Ipsi','ROI1_ACUToverBSL_Contra','ROI1_WK1overBSL_Ipsi','ROI1_WK1overBSL_Contra', ...
    'ROI2_bsl_Ipsi','ROI2_bsl_Contra','ROI2_acute_Ipsi','ROI2_acute_Contra','ROI2_wk1_Ipsi','ROI2_wk1_Contra', ...
    'ROI2_Pixels','ROI2_ContraPixels','ROI2_CentroidX','ROI2_CentroidY','ROI2_ContraCentroidX','ROI2_ContraCentroidY', ...
    'ROI2_ACUToverBSL_Ipsi','ROI2_ACUToverBSL_Contra','ROI2_WK1overBSL_Ipsi','ROI2_WK1overBSL_Contra'};

fields = headers(2:end);  % everything except Mouse corresponds to SO_values fields

%% ---------------- BUILD ROWDATA (NaN-safe) ----------------
rowData = cell(1, numel(headers));
rowData{1} = mouseNumN;

for k = 1:numel(fields)
    rowData{k+1} = getFieldOrNaN(SO_values, fields{k});
end

% sanity check
if numel(rowData) ~= numel(headers)
    error("Header/row length mismatch: headers=%d rowData=%d", numel(headers), numel(rowData));
end

%% ---------------- PREVIEW TABLE ----------------
fprintf("\n================ PREVIEW: DATA TO APPEND =================\n");

T = cell2table(rowData, 'VariableNames', matlab.lang.makeValidName(headers));
disp(T);

if any(isnan(table2array(T(:,2:end))), 'all')
    warning("⚠️ Some numeric values are NaN. If unexpected, debug ROI masks / maps.");
end
fprintf("----------------------------------------------------------\n");

resp = input("❓ Append these values to Excel? (1 = YES, 0 = NO): ");
if resp ~= 1
    fprintf("🚫 Not writing to Excel. Exiting.\n");
    return;
end

%% ---------------- WRITE/APPEND TO EXCEL ----------------
if ~isfile(excelPath)
    writecell(headers, excelPath, 'Sheet', 1, 'Range', 'A1');
    nextRow = 2;
else
    [~, ~, raw] = xlsread(excelPath);
    nextRow = size(raw, 1) + 1;

    existingHeaders = raw(1,:);
    if numel(existingHeaders) ~= numel(headers) || any(~strcmp(existingHeaders, headers))
        error("Excel headers do not match expected format. Delete ROI_Summary_V2.xlsx or fix header row.");
    end
end

writecell(rowData, excelPath, 'Sheet', 1, 'Range', sprintf('A%d', nextRow));
fprintf("\n💾 Results appended to %s (row %d)\n", excelPath, nextRow);

%% ---------------- VERIFY WRITE (read back) ----------------
[~,~,raw2] = xlsread(excelPath);
if size(raw2,1) >= nextRow
    writtenRow = raw2(nextRow, 1:numel(headers));
    fprintf("\n✅ Verified row %d read-back (first 6 cols):\n", nextRow);
    disp(writtenRow(1:min(6,numel(headers))));
else
    warning("⚠️ Could not verify write-back. Excel may be locked/open or network delay.");
end


%% ==================== helper function ====================
function v = getFieldOrNaN(S, f)
    f = char(f);
    if isfield(S, f) && ~isempty(S.(f))
        v = S.(f);
        if ~isscalar(v)
            v = v(1);
        end
        if ~isfinite(v)
            v = NaN;
        end
    else
        v = NaN;
    end
end