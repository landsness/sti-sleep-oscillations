%% ===============================================================
%  SO_GroupTopoplots_AverageFromRawMaps_SavePNGs_v2_arn.m
% Purpose:
%   - Ask user for a list of mouse numbers (e.g., "11 12 13 14 15 16 17")
%   - Auto-find each Mouse_##_SO_RawMaps.mat anywhere under Z:\Results\
%   - Load so_map_bsl / so_map_acute / so_map_wk1 for each mouse
%   - Compute GROUP MEAN maps for BSL, Acute, WK1 across selected mice
%   - Display ONLY the 3 group-average topoplots (no per-mouse saving)
%  - Accept mouse IDs as strings (supports: 11, 2DBSI, 4DBSI, etc.)
%  - Auto-find each Mouse_<ID>_SO_RawMaps.mat under Z:\Results\
%  - Compute group mean maps (BSL/Acute/WK1)
%  - SAVE 3 separate PNGs + VERIFY save
%  - NEW: set 0-valued background to WHITE using custom colormap
% ===============================================================

clear; clc; close all;

%% ---------------- USER INPUTS ----------------
resultsRoot = "Z:\Results\";

mouseStr = input("Enter mouse IDs (space/comma separated, e.g., 11 12 2DBSI 4DBSI): ", "s");
mouseIDs = parseIDList(mouseStr);

if isempty(mouseIDs)
    error("No mouse IDs parsed. Example: 11 12 2DBSI 4DBSI");
end

fprintf("\nMouse IDs to include (%d): %s\n", numel(mouseIDs), strjoin(mouseIDs, ", "));

% Output folder (auto-created)
stamp = datestr(now, "yyyymmdd_HHMMSS");
outRoot = fullfile(resultsRoot, "SO_Maps_GroupPlots");
outDir  = fullfile(outRoot, "GroupMean_" + stamp);
if ~isfolder(outDir), mkdir(outDir); end

fprintf("📁 Output folder:\n  %s\n", outDir);

%% ---------------- FIND + LOAD + ACCUMULATE ----------------
reqVars = ["so_map_bsl","so_map_acute","so_map_wk1"];

bslStack   = [];
acuteStack = [];
wk1Stack   = [];

usedPaths = strings(0);
usedIDs   = strings(0);

for i = 1:numel(mouseIDs)
    id = mouseIDs(i);
    pattern = "Mouse_" + id + "_SO_RawMaps.mat";

    fprintf("\n🔍 Searching for %s under %s ...\n", pattern, resultsRoot);
    fileList = dir(fullfile(resultsRoot, "**", pattern));

    if isempty(fileList)
        warning("❌ Not found: %s (skipping %s)", pattern, id);
        continue;
    end

    % choose newest match
    [~, idx] = max([fileList.datenum]);
    rawMapsPath = fullfile(fileList(idx).folder, fileList(idx).name);
    fprintf("✅ Using %s\n", rawMapsPath);

    S = load(rawMapsPath);

    % confirm required vars exist
    ok = true;
    for r = reqVars
        if ~isfield(S, r)
            warning("❌ %s missing var %s (skipping %s)", rawMapsPath, r, id);
            ok = false; break;
        end
    end
    if ~ok, continue; end

    bsl   = S.so_map_bsl;
    acute = S.so_map_acute;
    wk1   = S.so_map_wk1;

    % sanity checks
    if ~ismatrix(bsl) || ~ismatrix(acute) || ~ismatrix(wk1)
        warning("❌ Non-2D map in %s (skipping %s)", rawMapsPath, id);
        continue;
    end
    if any(size(bsl) ~= size(acute)) || any(size(bsl) ~= size(wk1))
        warning("❌ Size mismatch in %s (skipping %s)", rawMapsPath, id);
        continue;
    end

    % init stacks
    if isempty(bslStack)
        bslStack   = zeros([size(bsl), 0], "like", bsl);
        acuteStack = zeros([size(acute), 0], "like", acute);
        wk1Stack   = zeros([size(wk1), 0], "like", wk1);
    else
        if any(size(bsl) ~= size(bslStack(:,:,1)))
            warning("❌ Map dim mismatch vs earlier mice in %s (skipping %s)", rawMapsPath, id);
            continue;
        end
    end

    % append
    bslStack(:,:,end+1)   = bsl;
    acuteStack(:,:,end+1) = acute;
    wk1Stack(:,:,end+1)   = wk1;

    usedPaths(end+1) = string(rawMapsPath);
    usedIDs(end+1)   = id;

    fprintf("➕ Added %s | N=%d\n", id, size(bslStack,3));
end

if isempty(bslStack) || size(bslStack,3) < 1
    error("No mice successfully loaded. Check resultsRoot / naming / permissions.");
end

N = size(bslStack,3);

fprintf("\n================ SUMMARY ================\n");
fprintf("Included IDs (N=%d): %s\n", numel(usedIDs), strjoin(usedIDs, ", "));
fprintf("Stack size: %dx%dx%d\n", size(bslStack,1), size(bslStack,2), N);
fprintf("=========================================\n\n");

%% ---------------- GROUP MEANS ----------------
mean_bsl   = mean(bslStack,   3, "omitnan");
mean_acute = mean(acuteStack, 3, "omitnan");
mean_wk1   = mean(wk1Stack,   3, "omitnan");

% shared scaling for comparability
allVals = [mean_bsl(:); mean_acute(:); mean_wk1(:)];
allVals = allVals(isfinite(allVals));
clim = prctile(allVals, [2 98]);

% NEW: custom colormap where 0 -> white (first color)
cmap = jet(256);
cmap(1,:) = [1 1 1];

% NEW: force lower bound at 0 so zeros map to cmap(1,:)
clim0 = [0 clim(2)];

%% ---------------- DISPLAY FIGURE ----------------
fig = figure("Color","w", "Name","Group Mean SO Topoplots");
tiledlayout(1,3,"Padding","compact","TileSpacing","compact");

nexttile; imagesc(mean_bsl); axis image off; colormap(cmap); colorbar; caxis(clim0);
title(sprintf("Group Mean BSL (N=%d)", N), "FontWeight","bold");

nexttile; imagesc(mean_acute); axis image off; colormap(cmap); colorbar; caxis(clim0);
title(sprintf("Group Mean Acute (N=%d)", N), "FontWeight","bold");

nexttile; imagesc(mean_wk1); axis image off; colormap(cmap); colorbar; caxis(clim0);
title(sprintf("Group Mean WK1 (N=%d)", N), "FontWeight","bold");

drawnow;

%% ---------------- SAVE 3 SEPARATE PNGs + VERIFY ----------------
% We save each panel by re-rendering into a clean single-axis figure so the PNG
% is exactly what Eric expects (no multi-panel confusion).

savePanel(mean_bsl,   clim0, fullfile(outDir, sprintf("GroupMean_BSL_N%d.png",   N)), sprintf("Group Mean BSL (N=%d)", N));
savePanel(mean_acute, clim0, fullfile(outDir, sprintf("GroupMean_Acute_N%d.png", N)), sprintf("Group Mean Acute (N=%d)", N));
savePanel(mean_wk1,   clim0, fullfile(outDir, sprintf("GroupMean_WK1_N%d.png",   N)), sprintf("Group Mean WK1 (N=%d)", N));

% also save a text log of included mice + file paths (super helpful for audit)
logPath = fullfile(outDir, "IncludedMice_Log.txt");
fid = fopen(logPath, "w");
fprintf(fid, "Group mean generated: %s\n", datestr(now));
fprintf(fid, "Included IDs (N=%d): %s\n\n", numel(usedIDs), strjoin(usedIDs, ", "));
fprintf(fid, "Files used:\n");
for k = 1:numel(usedPaths)
    fprintf(fid, "  %s : %s\n", usedIDs(k), usedPaths(k));
end
fclose(fid);

fprintf("\n✅ Saved outputs:\n");
fprintf("  %s\n", outDir);
fprintf("  %s\n", logPath);

if isfile(logPath)
    fprintf("✅ Verified log exists.\n");
else
    warning("⚠️ Log did not verify (permissions/network delay?)");
end

%% ==================== helpers ====================

function ids = parseIDList(s)
    s = strrep(s, ",", " ");
    s = strrep(s, ";", " ");
    parts = strsplit(strtrim(s));
    ids = strings(0);
    for i = 1:numel(parts)
        p = strtrim(parts{i});
        if strlength(p)==0, continue; end
        ids(end+1) = string(p); %#ok<AGROW>
    end
    % keep stable order, drop duplicates
    [~, ia] = unique(ids, "stable");
    ids = ids(sort(ia));
end

function savePanel(map, clim, outPath, ttl)
    f = figure("Color","w","Visible","off");
    ax = axes(f);

    imagesc(ax, map);
    axis(ax, "image"); axis(ax, "off");

    % NEW: custom colormap where 0 -> white
    cmap = jet(256);
    cmap(1,:) = [1 1 1];
    colormap(ax, cmap);

    caxis(ax, clim);
    colorbar(ax);
    title(ax, ttl, "FontWeight","bold");

    % High-res export
    exportgraphics(f, outPath, "Resolution", 300);

    close(f);

    % verification print
    if isfile(outPath)
        fprintf("✅ Saved: %s\n", outPath);
    else
        warning("❌ Save failed: %s", outPath);
    end
end
