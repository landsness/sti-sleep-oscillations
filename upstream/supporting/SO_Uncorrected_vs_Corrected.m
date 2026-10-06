%% ===============================================================
%  SO_Uncorrected_vs_Corrected.m
%  Extract SO (0.1-1 Hz) power maps from UNCORRECTED (xform_datafluor) and
%  CORRECTED (xform_datafluorCorr) GCaMP, to show the hemodynamic correction
%  removes the PT "bright artifact" contribution at the ischemic core.
%
%  Also saves a baseline-brightness proxy (temporal mean magnitude) to
%  visualize the bright core in the uncorrected signal.
%
%  Uses the lab's own PowerAnalysis.m so processing is identical to the
%  pipeline. Loads the big dataFluor.mat only for the mice you list.
%
%  Layout: BASE_DIR/<mouse>/{bsl,acute,oneweek}/
%     <date>-<mouse>-fc1-dataFluor.mat   (xform_datafluor, xform_datafluorCorr, runInfo)
%     <date>-<mouse>-LandmarksAndMask.mat (xform_isbrain)
% ===============================================================
clear; clc; close all;

%% ---------------- USER INPUT ----------------
BASE_DIR = 'C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\PreProcessed_Data';
OUT_DIR  = 'C:\Users\landsness\Box\Manuscripts\STI_SO\Revisions\SO_ContrastMaps';
MICE     = {'11','13'};      % start with large-infarct mice where artifact is most visible
% add the LeeWrapper2Cam Scripts folder (for PowerAnalysis.m) to the path:
addpath('C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\SO_Analysis\MATLAB_Code_Used\LeeWrapper2Cam Scripts');

SO_LOW = 0.1; SO_HIGH = 1.0;
sessions = {'bsl',{'bsl'}; 'acute',{'acute'}; 'wk1',{'oneweek','weekone','wk1','week1'}};
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

for mi = 1:numel(MICE)
    m = MICE{mi};
    fprintf('\n== Mouse %s ==\n', m);
    U = struct(); Cc = struct(); B = struct();   % uncorrected SO / corrected SO / brightness
    for si = 1:size(sessions,1)
        sess = sessions{si,1};
        folder = '';
        for f = sessions{si,2}
            cand = fullfile(BASE_DIR,m,f{1}); if exist(cand,'dir'), folder=cand; break; end
        end
        if isempty(folder), continue; end
        ff = dir(fullfile(folder,'*dataFluor.mat'));
        if isempty(ff), fprintf('  [%s] no dataFluor.mat\n',sess); continue; end
        mf = dir(fullfile(folder,'*LandmarksAndMask.mat'));

        S = load(fullfile(folder,ff(1).name),'xform_datafluor','xform_datafluorCorr','runInfo');
        fs = S.runInfo.samplingRate;
        mask = true(128,128);
        if ~isempty(mf)
            M = load(fullfile(folder,mf(1).name),'xform_isbrain');
            if isfield(M,'xform_isbrain'), mask = ~isnan(M.xform_isbrain) & (M.xform_isbrain>0.5); end
        end

        % --- SO power maps via the lab's PowerAnalysis (0.1-1 Hz mean) ---
        [wsU,~,~,hz] = PowerAnalysis(double(S.xform_datafluor),     fs, mask);
        [wsC,~,~,~ ] = PowerAnalysis(double(S.xform_datafluorCorr), fs, mask);
        band = hz>SO_LOW & hz<SO_HIGH;
        soU = mean(wsU(:,:,band),3); soU(~mask)=NaN;
        soC = mean(wsC(:,:,band),3); soC(~mask)=NaN;

        % --- brightness proxy: temporal mean of |uncorrected| (shows bright core) ---
        bright = mean(abs(double(S.xform_datafluor)),3,'omitnan'); bright(~mask)=NaN;

        U.(sess)=soU; Cc.(sess)=soC; B.(sess)=bright;
        fprintf('  [%s] fs=%.2f Hz  done\n', sess, fs);
        clear S wsU wsC
    end
    save(fullfile(OUT_DIR, sprintf('Mouse_%s_SO_UncorrVsCorr.mat', m)), 'U','Cc','B');
    fprintf('  saved Mouse_%s_SO_UncorrVsCorr.mat\n', m);
end
disp('Done. Stage the *_SO_UncorrVsCorr.mat files back for the composite figure.');
