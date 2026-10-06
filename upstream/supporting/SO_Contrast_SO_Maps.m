%% ===============================================================
%  SO_Contrast_SO_Maps.m
%  Extract 0.1-1 Hz SO power maps for ALL contrasts (HbO/HbR/HbT/Calcium)
%  from the existing *-Power.mat files, to test whether the lesional SO
%  "suppression" seen in calcium also appears in the hemodynamic channels.
%
%  Reuses the same Power.mat / whole_spectra_map the calcium SO maps came
%  from, so the method is identical by construction. No reprocessing of
%  dataFluor / dataHb needed.
%
%  whole_spectra_map(:,:,freq,contrast) contrasts (LeeWrapper lines 301-305):
%     1 = HbO   2 = HbR   3 = HbT   4 = Calcium (hemodynamically corrected)
% ===============================================================
clear; clc; close all;

%% ---------------- USER INPUT (same as SO_PowerAnalysis_Arnav_v2) ----------------
excelFile = "Z:\Databases\DataBase_DBSI_arn.xlsx";
excelRows.bsl   = [172 173 174 175 176 177];
excelRows.acute = [179 180 181 182 183 184];
excelRows.wk1   = [186 187 188 189 190 191];

SO_LOW = 0.1;  SO_HIGH = 1.0;
CONTRASTS = struct('idx',{1,2,3,4},'name',{'HbO','HbR','HbT','Calcium'});
ratio_caxis = [0 2];

%% ---------------- LOAD RUN INFO ----------------
runsInfo = struct();
runsInfo.bsl   = parseRunsJM(excelFile, excelRows.bsl);
runsInfo.acute = parseRunsJM(excelFile, excelRows.acute);
runsInfo.wk1   = parseRunsJM(excelFile, excelRows.wk1);
mouseList = unique({runsInfo.bsl.mouseName});
sessionTypes = ["bsl","acute","wk1"];

for i = 1:numel(mouseList)
    mouseID = mouseList{i};
    bslRun = runsInfo.bsl(strcmp({runsInfo.bsl.mouseName}, mouseID)); bslRun = bslRun(1);
    outDir = fullfile(bslRun.saveFolder, 'ROI_Data');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    fprintf('\n== Mouse %s ==\n', mouseID);

    SO = struct();   % SO.(contrast).(session) = 128x128 linear SO power
    for st = sessionTypes
        runs = runsInfo.(st);
        r = runs(strcmp({runs.mouseName}, mouseID));
        if isempty(r), fprintf('  no %s run\n',st); continue; end
        r = r(1);
        pf = fullfile(r.saveFolder, sprintf('%s-%s-fc1-Power.mat', r.recDate, r.mouseName));
        if ~exist(pf,'file'), fprintf('  missing %s\n',pf); continue; end
        S = load(pf,'whole_spectra_map','hz');
        SO_Hz = find(S.hz > SO_LOW & S.hz < SO_HIGH);
        for c = CONTRASTS
            SO.(c.name).(st) = mean(S.whole_spectra_map(:,:,SO_Hz,c.idx), 3);  % linear SO power
        end
        fprintf('  %s done\n', st);
    end

    % --- save raw maps (all contrasts) ---
    save(fullfile(outDir, sprintf('Mouse_%s_SO_AllContrasts.mat', mouseID)), 'SO');

    % --- ratio maps + side-by-side figure (acute/bsl and wk1/bsl per contrast) ---
    fig = figure('Visible','off','Position',[100 100 1100 620]);
    tl = tiledlayout(2,numel(CONTRASTS),'TileSpacing','compact','Padding','compact');
    for c = CONTRASTS
        for k = 1:2
            st = {'acute','wk1'}; st = st{k};
            nexttile;
            if isfield(SO,c.name) && isfield(SO.(c.name),'bsl') && isfield(SO.(c.name),st)
                b = SO.(c.name).bsl; b(b==0)=NaN;
                m = SO.(c.name).(st); m(m==0)=NaN;
                imagesc(m./b,'AlphaData',~isnan(m./b)); caxis(ratio_caxis);
            end
            axis image off; colormap(jet);
            title(sprintf('%s  %s/bsl', c.name, st),'Interpreter','none');
        end
    end
    cb=colorbar; cb.Layout.Tile='east';
    title(tl, sprintf('Mouse %s - SO (0.1-1 Hz) ratio maps by contrast', mouseID),'Interpreter','none');
    saveas(fig, fullfile(outDir, sprintf('Mouse_%s_SO_ContrastRatios.png', mouseID)));
    close(fig);
end
disp('Done. Compare the HbT/HbO/HbR ratio maps against Calcium at the lesion.');
