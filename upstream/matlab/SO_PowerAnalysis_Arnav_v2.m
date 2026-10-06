%% ===============================================================
%  SO_PowerAnalysis_Arnav_v2.m
%  Author: Arnav Ajay Jadav | Landsness Lab, WashU
%
%  Purpose:
%   1) Load fc1 Power.mat for Baseline / Acute / Week1 from database
%   2) Compute slow oscillation power maps (0.1–1 Hz)
%   3) Save:
%       - PNGs + GIF in  <bsl saveFolder>\SO_Images\
%       - Numeric MATs in <bsl saveFolder>\ROI_Data\
%         * Mouse_<id>_SO_RawMaps.mat  (bsl/acute/wk1 maps)
%         * Mouse_<id>_SO_Recovery.mat (wk1/bsl - acute/bsl)
%   4) Keep your old visual outputs, but *prioritize numeric maps* for ROI
% ===============================================================

clear; clc; close all;

%% ---------------- USER INPUT ----------------
excelFile = "Z:\Databases\DataBase_DBSI_arn.xlsx";

% Put [] if you want to temporarily disable any session
excelRows.bsl   = [172, 173, 174, 175, 176, 177];
excelRows.acute = [179, 180, 181, 182, 183, 184];
excelRows.wk1   = [186, 187, 188, 189, 190, 191];

% Plot scaling
caxis_range = [-60 -40];   % dB power map range for display
ratio_caxis = [0 2];       % ratio image range

% Frequency band
SO_LOW = 0.1;
SO_HIGH = 1.0;
SPECTRA_CHAN = 4;          % you used index 4 previously; keep same

%% ---------------- LOAD RUN INFO ----------------
runsInfo = struct();
if ~isempty(excelRows.bsl),   runsInfo.bsl   = parseRunsJM(excelFile, excelRows.bsl); end
if ~isempty(excelRows.acute), runsInfo.acute = parseRunsJM(excelFile, excelRows.acute); end
if ~isempty(excelRows.wk1),   runsInfo.wk1   = parseRunsJM(excelFile, excelRows.wk1); end

if ~isfield(runsInfo, 'bsl') || isempty(runsInfo.bsl)
    error('Baseline rows are empty / missing. Need baseline to anchor save directory & mouse list.');
end

% Unique mice anchored from baseline entries
mouseList = unique({runsInfo.bsl.mouseName});

%% ---------------- MAIN LOOP ----------------
for i = 1:numel(mouseList)

    mouseID = mouseList{i};
    fprintf('\n==============================\n');
    fprintf('🐭 Processing mouse: %s\n', mouseID);

    % --- Grab baseline run (anchors folder organization)
    bslRun = runsInfo.bsl(strcmp({runsInfo.bsl.mouseName}, mouseID));
    if isempty(bslRun)
        fprintf('⚠️ No baseline run found for mouse %s. Skipping.\n', mouseID);
        continue;
    end
    bslRun = bslRun(1); % in case multiple, take first

    % Output folders (inside baseline date folder)
    soImagesDir = fullfile(bslRun.saveFolder, 'SO_Images');
    roiDataDir  = fullfile(bslRun.saveFolder, 'ROI_Data');

    if ~exist(soImagesDir, 'dir'), mkdir(soImagesDir); end
    if ~exist(roiDataDir, 'dir'),  mkdir(roiDataDir);  end

    % Store raw numeric SO maps here
    powerData = struct();        % holds linear (not dB) mean SO power maps
    meta = struct();
    meta.mouseID = mouseID;
    meta.SO_band_Hz = [SO_LOW SO_HIGH];

    % --- Pull each session run and compute SO map
    sessionTypes = ["bsl","acute","wk1"];

    for sessionType = sessionTypes

        if ~isfield(runsInfo, sessionType) || isempty(runsInfo.(sessionType))
            continue;
        end

        sessionRuns = runsInfo.(sessionType);

        % Find this mouse in that session
        sessionRun = [];
        for j = 1:numel(sessionRuns)
            if strcmp(sessionRuns(j).mouseName, mouseID)
                sessionRun = sessionRuns(j);
                break;
            end
        end

        if isempty(sessionRun)
            fprintf('⚠️ No %s run found for mouse %s.\n', upper(sessionType), mouseID);
            continue;
        end

        % Build power file path
        powerFile = fullfile(sessionRun.saveFolder, ...
            sprintf('%s-%s-fc1-Power.mat', sessionRun.recDate, sessionRun.mouseName));

        if ~exist(powerFile, 'file')
            fprintf('❌ Missing power file: %s\n', powerFile);
            continue;
        end

        % Load spectra
        S = load(powerFile, 'whole_spectra_map', 'hz');
        if ~isfield(S, 'whole_spectra_map') || ~isfield(S, 'hz')
            fprintf('❌ File missing expected variables whole_spectra_map/hz: %s\n', powerFile);
            continue;
        end

        hz = S.hz;
        whole_spectra_map = S.whole_spectra_map;

        SO_Hz = find(hz > SO_LOW & hz < SO_HIGH);
        if isempty(SO_Hz)
            fprintf('❌ No frequencies in SO band for %s (%s). Skipping.\n', mouseID, sessionType);
            continue;
        end

        % Mean SO power map (linear power)
        so_map = mean(whole_spectra_map(:,:,SO_Hz,SPECTRA_CHAN), 3);

        powerData.(sessionType) = so_map;

        % Metadata
        meta.(sessionType).recDate    = sessionRun.recDate;
        meta.(sessionType).saveFolder = sessionRun.saveFolder;
        meta.(sessionType).powerFile  = powerFile;

        % ---- Save PNG (display in dB)
        fig = figure('Visible','off');
        imagesc(10*log10(so_map));
        colormap([1 1 1; jet]); colorbar;
        caxis(caxis_range);
        axis image off;
        title(sprintf('Mouse %s - %s Slow Oscillation Power', mouseID, upper(sessionType)));

        pngName = fullfile(soImagesDir, sprintf('Mouse_%s_%s_SO.png', mouseID, sessionType));
        saveas(fig, pngName);
        close(fig);

        fprintf('✅ %s map done.\n', upper(sessionType));
    end

    % --- Save RAW MAPS MAT (this is what ROI quant script should load)
    rawMapsPath = fullfile(roiDataDir, sprintf('Mouse_%s_SO_RawMaps.mat', mouseID));

    so_map_bsl   = [];
    so_map_acute = [];
    so_map_wk1   = [];

    if isfield(powerData,'bsl'),   so_map_bsl   = powerData.bsl; end
    if isfield(powerData,'acute'), so_map_acute = powerData.acute; end
    if isfield(powerData,'wk1'),   so_map_wk1   = powerData.wk1; end

    save(rawMapsPath, 'so_map_bsl','so_map_acute','so_map_wk1','meta');
    fprintf('💾 Saved raw maps MAT: %s\n', rawMapsPath);

    % --- Ratios + Recovery (only if baseline + acute + wk1 exist)
    hasBSL = isfield(powerData,'bsl')   && ~isempty(powerData.bsl);
    hasAc  = isfield(powerData,'acute') && ~isempty(powerData.acute);
    hasW1  = isfield(powerData,'wk1')   && ~isempty(powerData.wk1);

    if hasBSL

        % protect baseline zeros
        bsl_no_zeros = powerData.bsl;
        bsl_no_zeros(bsl_no_zeros == 0) = NaN;

        if hasAc
            acute_no_zeros = powerData.acute;
            acute_no_zeros(acute_no_zeros == 0) = NaN;

            acute_bsl = acute_no_zeros ./ bsl_no_zeros;

            fig = figure('Visible','off');
            imagesc(acute_bsl); colormap([1 1 1; jet]); colorbar;
            caxis(ratio_caxis); axis image off;
            title(sprintf('Mouse %s - Acute/Baseline Ratio', mouseID));

            pngName = fullfile(soImagesDir, sprintf('Mouse_%s_Acute_Baseline_Ratio.png', mouseID));
            saveas(fig, pngName);
            close(fig);
        else
            acute_bsl = [];
        end

        if hasW1
            wk1_no_zeros = powerData.wk1;
            wk1_no_zeros(wk1_no_zeros == 0) = NaN;

            wk1_bsl = wk1_no_zeros ./ bsl_no_zeros;

            fig = figure('Visible','off');
            imagesc(wk1_bsl); colormap([1 1 1; jet]); colorbar;
            caxis(ratio_caxis); axis image off;
            title(sprintf('Mouse %s - Week1/Baseline Ratio', mouseID));

            pngName = fullfile(soImagesDir, sprintf('Mouse_%s_Week1_Baseline_Ratio.png', mouseID));
            saveas(fig, pngName);
            close(fig);
        else
            wk1_bsl = [];
        end

        % --- Recovery map: (wk1/bsl) - (acute/bsl)
        if ~isempty(wk1_bsl) && ~isempty(acute_bsl)
            so_recovery_map = wk1_bsl - acute_bsl;

            fig = figure('Visible','off');
            imagesc(so_recovery_map);
            colormap(jet); colorbar;
            caxis([-0.5 0.5]);
            axis image off;
            title(sprintf('Mouse %s - SO Recovery Map (Wk1 - Acute)', mouseID));

            pngName = fullfile(soImagesDir, sprintf('Mouse_%s_SO_Recovery.png', mouseID));
            saveas(fig, pngName);
            close(fig);

            % Save numeric recovery map MAT into ROI_Data
            recPath = fullfile(roiDataDir, sprintf('Mouse_%s_SO_Recovery.mat', mouseID));
            save(recPath, 'so_recovery_map','meta');
            fprintf('💾 Saved recovery MAT: %s\n', recPath);
        end
    end

    % --- GIF evolution (if any maps exist)
    gifFile = fullfile(soImagesDir, sprintf('Mouse_%s_SO_Evolution.gif', mouseID));
    fig = figure('Visible','off');

    wroteFirst = false;
    for t = sessionTypes
        if isfield(powerData, t) && ~isempty(powerData.(t))
            imagesc(10*log10(powerData.(t)));
            colormap(jet); colorbar; caxis(caxis_range);
            axis image off;
            title(sprintf('Mouse %s - %s SO Power', mouseID, upper(t)));

            frame = getframe(fig);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);

            if ~wroteFirst
                imwrite(imind, cm, gifFile, 'gif', 'Loopcount', inf, 'DelayTime', 1);
                wroteFirst = true;
            else
                imwrite(imind, cm, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
            end
        end
    end
    close(fig);

    fprintf('🎞️ Saved GIF (if maps existed): %s\n', gifFile);
    fprintf('✅ Done mouse %s\n', mouseID);

end

disp('✅ Processing Complete. Raw maps + images saved.');