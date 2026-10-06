%% ===============================================================
%  Extract_ROI_Spectra.m
%  Export ROI-averaged POWER SPECTRA (PSD vs frequency) for the
%  lesional (ROI1), perilesional (ROI2), and contralateral-homolog ROIs,
%  at baseline / acute / week-1, for every mouse.
%
%  Answers Reviewer's request for "examples of power spectra so readers
%  can evaluate the extracted SO power metrics."  Reuses the SAME
%  whole_spectra_map the SO metrics come from, so it is identical to the
%  pipeline by construction (no reprocessing).
%
%  Output is TINY (a few ROIs x freq x mouse) -> easy to stage back so the
%  figure can be built/iterated in Python.
%
%  whole_spectra_map(:,:,freq,contrast):  4 = Calcium (hemodynamically
%  corrected) -- the channel the paper's SO metric uses.
%
%  Layout expected:
%    BASE_DIR/<mouse>/{bsl,acute,oneweek}/<date>-<mouse>-fc1-Power.mat
%    BASE_DIR/<mouse>/<session>/<date>-<mouse>-LandmarksAndMask.mat
% ===============================================================
clear; clc;

%% ---------------- USER INPUT ----------------
BASE_DIR = 'C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\PreProcessed_Data';
ROI_XLSX = 'C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\SO_Analysis\ROI_Summary_V2.xlsx';
OUT_DIR  = 'C:\Users\landsness\Box\Manuscripts\STI_SO\STI_SO_figures\data\derived';  % read by figures/panels/fig2c_roi_spectra.py
CALC_CHAN = 4;                 % calcium (corrected) channel in whole_spectra_map
SO_BAND  = [0.1 1.0];          % for reference lines only; full spectrum is exported
sessions = {'bsl',{'bsl'}; 'acute',{'acute'}; 'wk1',{'oneweek','weekone','wk1','week1'}};
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

%% ---------------- ROI centroids / sizes ----------------
% Mouse IDs are read as TEXT: the table mixes numeric WFCI IDs (1..24) with the
% two pilot animals '2DBSI' and '4DBSI'. (Rev. 2026-10-04: the earlier numeric
% read silently dropped both pilot animals, giving n = 23 instead of 25.)
opts = detectImportOptions(ROI_XLSX);
opts = setvartype(opts, 'Mouse', 'string');
T = readtable(ROI_XLSX, opts);
T.Mouse = strtrim(T.Mouse);
% expected columns: Mouse, ROI1_CentroidX, ROI1_CentroidY, ROI1_Pixels,
%                          ROI2_CentroidX, ROI2_CentroidY, ROI2_Pixels
getrow = @(mid) T(T.Mouse==string(mid),:);

%% ---------------- walk mice ----------------
d = dir(BASE_DIR);
mouseDirs = {d([d.isdir] & ~ismember({d.name},{'.','..'})).name};

hz = [];                       % filled from first file
EXPORT = struct('mouse',{},'session',{},'roi',{},'psd',{});
rowk = 0;

for mi = 1:numel(mouseDirs)
    mouseID = mouseDirs{mi};
    % folder -> ROI-table ID:  '7' -> '7',  'Mouse_2DBSI' -> '2DBSI'
    mid = regexprep(mouseID, '^Mouse_', '');
    r = getrow(mid);
    if isempty(r), fprintf('Folder %s: no ROI row, skipping\n', mouseID); continue; end
    if height(r) > 1, error('Folder %s: %d ROI rows match ID %s', mouseID, height(r), mid); end

    cx1 = r.ROI1_CentroidX(1); cy1 = r.ROI1_CentroidY(1); n1 = r.ROI1_Pixels(1);
    cx2 = r.ROI2_CentroidX(1); cy2 = r.ROI2_CentroidY(1); n2 = r.ROI2_Pixels(1);
    rad1 = sqrt(n1/pi); rad2 = sqrt(n2/pi);
    % contralateral homolog = mirror ROI1 across midline (128-grid, bregma-centered)
    cx1c = 129 - cx1; cy1c = cy1;

    for si = 1:size(sessions,1)
        sess = sessions{si,1}; folder = '';
        for f = sessions{si,2}
            cand = fullfile(BASE_DIR, mouseID, f{1});
            if exist(cand,'dir'), folder = cand; break; end
        end
        if isempty(folder), continue; end
        pf = dir(fullfile(folder,'*Power.mat'));
        if isempty(pf), fprintf('Mouse %s [%s]: no Power.mat\n',mouseID,sess); continue; end

        S = load(fullfile(folder,pf(1).name),'whole_spectra_map','hz');
        if ~isfield(S,'whole_spectra_map') || ~isfield(S,'hz'), continue; end
        if isempty(hz), hz = S.hz(:); end
        wsm = squeeze(S.whole_spectra_map(:,:,:,CALC_CHAN));   % 128x128xF

        % brain mask (optional)
        mask = true(128,128);
        mf = dir(fullfile(folder,'*LandmarksAndMask.mat'));
        if ~isempty(mf)
            M = load(fullfile(folder,mf(1).name),'xform_isbrain');
            if isfield(M,'xform_isbrain'), mask = ~isnan(M.xform_isbrain) & (M.xform_isbrain>0.5); end
        end

        rois = { 'lesional',   cx1,  cy1,  rad1;
                 'perilesional',cx2, cy2,  rad2;
                 'contra',     cx1c, cy1c, rad1 };
        [X,Y] = meshgrid(1:128,1:128);
        for ri = 1:size(rois,1)
            nm = rois{ri,1}; cx = rois{ri,2}; cy = rois{ri,3}; rad = rois{ri,4};
            disk = ((X-cx).^2 + (Y-cy).^2) <= rad^2 & mask;
            if nnz(disk)==0, continue; end
            psd = zeros(numel(hz),1);
            for k = 1:numel(hz)
                sl = wsm(:,:,k);
                psd(k) = mean(sl(disk),'omitnan');
            end
            rowk = rowk + 1;
            EXPORT(rowk).mouse   = char(mid);   % text ID, e.g. '7' or '2DBSI'
            EXPORT(rowk).session = sess;
            EXPORT(rowk).roi     = nm;
            EXPORT(rowk).psd     = psd(:)';
        end
        fprintf('Mouse %s [%s]: ok (%d freq bins)\n', mouseID, sess, numel(hz));
        clear S wsm
    end
end

%% ---------------- save tiny export ----------------
outfile = fullfile(OUT_DIR,'ROI_Spectra_export.mat');
save(outfile,'EXPORT','hz','SO_BAND','CALC_CHAN','-v7');   % -v7 = small, python-readable via scipy
fprintf('\nSaved %s  (%d rows, %d freq bins, %d mice)\n', outfile, numel(EXPORT), numel(hz), ...
        numel(unique({EXPORT.mouse})));   % expect 25 mice

% also a flat CSV (mouse,session,roi,f1,f2,...) for convenience
csvfile = fullfile(OUT_DIR,'ROI_Spectra_export.csv');
fid = fopen(csvfile,'w');
fprintf(fid,'mouse,session,roi'); fprintf(fid,',f%.5g',hz); fprintf(fid,'\n');
for i=1:numel(EXPORT)
    fprintf(fid,'%s,%s,%s',EXPORT(i).mouse,EXPORT(i).session,EXPORT(i).roi);
    fprintf(fid,',%.6g',EXPORT(i).psd); fprintf(fid,'\n');
end
fclose(fid);
fprintf('Saved %s\n', csvfile);
disp('Done. Stage ROI_Spectra_export.mat (or .csv) back so the figure can be built.');
