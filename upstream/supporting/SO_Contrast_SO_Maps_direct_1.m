%% ===============================================================
%  SO_Contrast_SO_Maps_direct.m
%  Extract 0.1-1 Hz SO power maps for ALL contrasts (HbO/HbR/HbT/Calcium)
%  by loading the *-Power.mat files DIRECTLY from the folder tree.
%  No Excel database / parseRunsJM required.
%
%  Applies the brain mask (xform_isbrain) from the LandmarksAndMask file
%  (runInfo.saveMaskFile) so non-brain edge artifacts are removed (set NaN).
%
%  Expected layout:
%    BASE_DIR / <mouse> / {bsl, acute, oneweek} /
%        <date>-<mouse>-fc1-Power.mat
%        <date>-<mouse>-LandmarksAndMask.mat
%
%  whole_spectra_map(:,:,freq,contrast) contrast order (LeeWrapper 301-305):
%     1 = HbO   2 = HbR   3 = HbT   4 = Calcium (hemodynamically corrected)
% ===============================================================
clear; clc; close all;

%% ---------------- USER INPUT ----------------
BASE_DIR = 'C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\PreProcessed_Data';
OUT_DIR  = 'C:\Users\landsness\Box\Manuscripts\STI_SO\Revisions\SO_ContrastMaps';

SO_LOW = 0.1;  SO_HIGH = 1.0;
contrastNames = {'HbO','HbR','HbT','Calcium'};     % channels 1..4
sessions = {'bsl',{'bsl'}; 'acute',{'acute'}; 'wk1',{'oneweek','weekone','wk1','week1'}};
ratio_caxis = [0 2];

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

%% ---------------- FIND MOUSE FOLDERS ----------------
d = dir(BASE_DIR);
mouseDirs = {d([d.isdir] & ~ismember({d.name},{'.','..'})).name};

for mi = 1:numel(mouseDirs)
    mouseID = mouseDirs{mi};
    fprintf('\n== Mouse %s ==\n', mouseID);
    SO = struct();
    gotAny = false;

    for si = 1:size(sessions,1)
        sessName = sessions{si,1};
        folder = '';
        for f = sessions{si,2}
            cand = fullfile(BASE_DIR, mouseID, f{1});
            if exist(cand,'dir'), folder = cand; break; end
        end
        if isempty(folder), continue; end

        pf = dir(fullfile(folder, '*Power.mat'));
        if isempty(pf), fprintf('  [%s] no Power.mat\n', sessName); continue; end
        powerFile = fullfile(folder, pf(1).name);

        S = load(powerFile, 'whole_spectra_map', 'hz');
        if ~isfield(S,'whole_spectra_map') || ~isfield(S,'hz')
            fprintf('  [%s] file missing whole_spectra_map/hz\n', sessName); continue;
        end
        hz = S.hz(:);
        SO_Hz = find(hz > SO_LOW & hz < SO_HIGH);
        nC = size(S.whole_spectra_map, 4);

        % --- brain mask from LandmarksAndMask (saveMaskFile) ---
        mask = [];
        mf = dir(fullfile(folder, '*LandmarksAndMask.mat'));
        if ~isempty(mf)
            M = load(fullfile(folder, mf(1).name), 'xform_isbrain');
            if isfield(M,'xform_isbrain')
                mask = ~isnan(M.xform_isbrain) & (M.xform_isbrain > 0.5);
            end
        end
        if isempty(mask)
            fprintf('  [%s] WARNING: no LandmarksAndMask/xform_isbrain; maps left unmasked\n', sessName);
        end

        for c = 1:min(nC, numel(contrastNames))
            sm = mean(S.whole_spectra_map(:,:,SO_Hz,c), 3);   % linear SO power
            if ~isempty(mask), sm(~mask) = NaN; end            % apply brain mask
            SO.(contrastNames{c}).(sessName) = sm;
        end
        gotAny = true;
        fprintf('  [%s] ok  (%d freq bins, %d contrasts, mask %s)\n', ...
            sessName, numel(SO_Hz), nC, string(~isempty(mask)));
    end

    if ~gotAny, continue; end

    save(fullfile(OUT_DIR, sprintf('Mouse_%s_SO_AllContrasts.mat', mouseID)), 'SO');

    % ---- ratio-map comparison figure (rows=contrast, cols=acute/bsl, wk1/bsl) ----
    fig = figure('Visible','off','Position',[100 100 620 240*numel(contrastNames)]);
    tl = tiledlayout(numel(contrastNames), 2, 'TileSpacing','compact','Padding','compact');
    for r = 1:numel(contrastNames)
        cn = contrastNames{r};
        for k = 1:2
            st = {'acute','wk1'}; st = st{k};
            nexttile;
            if isfield(SO,cn) && isfield(SO.(cn),'bsl') && isfield(SO.(cn),st)
                b = SO.(cn).bsl; b(b==0) = NaN;
                m = SO.(cn).(st); m(m==0) = NaN;
                ratio = m./b;
                imagesc(ratio, 'AlphaData', ~isnan(ratio)); caxis(ratio_caxis);
            end
            axis image off; colormap(jet);
            title(sprintf('%s  %s/bsl', cn, st), 'Interpreter','none');
        end
    end
    cb = colorbar; cb.Layout.Tile = 'east';
    title(tl, sprintf('Mouse %s - SO (0.1-1 Hz) ratio by contrast (brain-masked)', mouseID), 'Interpreter','none');
    saveas(fig, fullfile(OUT_DIR, sprintf('Mouse_%s_SO_ContrastRatios.png', mouseID)));
    close(fig);
end

fprintf('\nDone. Outputs in:\n  %s\n', OUT_DIR);
disp('Compare HbT/HbO/HbR ratio maps against Calcium at the lesion.');
