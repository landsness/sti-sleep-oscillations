%% ===============================================================
%  Extract_ROI_Traces.m
%  Export ROI-averaged TIME TRACES of uncorrected vs hemodynamically
%  corrected GCaMP (delta F/F), to answer the Reviewer's request to
%  "show hemodynamically corrected fluorescence traces."
%
%  For representative mouse(s), extracts a short time window of:
%     - uncorrected fluorescence  (xform_datafluor)
%     - corrected fluorescence    (xform_datafluorCorr)   [hemodynamic corr.]
%  averaged over lesional (ROI1), perilesional (ROI2), and contralateral
%  homolog ROIs, at baseline and acute.
%
%  Output is TINY (a few ROIs x ~60 s x 2 corrections) -> stage back and
%  the figure gets built/iterated in Python.
%
%  Layout expected:
%    BASE_DIR/<mouse>/{bsl,acute,oneweek}/<date>-<mouse>-fc1-dataFluor.mat
%       (xform_datafluor, xform_datafluorCorr, runInfo)
%    BASE_DIR/<mouse>/<session>/<date>-<mouse>-LandmarksAndMask.mat
% ===============================================================
clear; clc;

%% ---------------- USER INPUT ----------------
BASE_DIR = 'C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\PreProcessed_Data';
ROI_XLSX = 'C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\SO_Analysis\ROI_Summary_V2.xlsx';
OUT_DIR  = 'C:\Users\landsness\Box\Manuscripts\STI_SO\Revisions\SO_ContrastMaps';
MICE     = {'8','9','14','24'};  % near-median lesional acute suppression; pick best SNR from these
WINDOW_SEC = 60;                % length of exported trace segment
START_SEC  = 30;                % skip first N s (settling) before the window
sessions = {'bsl',{'bsl'}; 'acute',{'acute'}; 'wk1',{'oneweek','weekone','wk1','week1'}};   % bsl + acute + week 1
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

T = readtable(ROI_XLSX);
getrow = @(mnum) T(T.Mouse==mnum,:);

EXPORT = struct('mouse',{},'session',{},'roi',{},'corr',{},'fs',{},'trace',{});
rowk = 0;

for mi = 1:numel(MICE)
    mouseID = MICE{mi}; mnum = str2double(mouseID);
    r = getrow(mnum);
    if isempty(r), fprintf('Mouse %s: no ROI row\n',mouseID); continue; end
    cx1=r.ROI1_CentroidX(1); cy1=r.ROI1_CentroidY(1); n1=r.ROI1_Pixels(1);
    cx2=r.ROI2_CentroidX(1); cy2=r.ROI2_CentroidY(1); n2=r.ROI2_Pixels(1);
    rad1=sqrt(n1/pi); rad2=sqrt(n2/pi);
    cx1c=129-cx1; cy1c=cy1;

    for si=1:size(sessions,1)
        sess=sessions{si,1}; folder='';
        for f=sessions{si,2}
            cand=fullfile(BASE_DIR,mouseID,f{1}); if exist(cand,'dir'),folder=cand;break; end
        end
        if isempty(folder), continue; end
        ff=dir(fullfile(folder,'*dataFluor.mat'));
        if isempty(ff), fprintf('Mouse %s [%s]: no dataFluor\n',mouseID,sess); continue; end
        S=load(fullfile(folder,ff(1).name),'xform_datafluor','xform_datafluorCorr','runInfo');
        fs=S.runInfo.samplingRate;

        mask=true(128,128);
        mf=dir(fullfile(folder,'*LandmarksAndMask.mat'));
        if ~isempty(mf)
            M=load(fullfile(folder,mf(1).name),'xform_isbrain');
            if isfield(M,'xform_isbrain'), mask=~isnan(M.xform_isbrain)&(M.xform_isbrain>0.5); end
        end

        i0=max(1,round(START_SEC*fs)); i1=min(size(S.xform_datafluor,3), i0+round(WINDOW_SEC*fs)-1);
        idx=i0:i1;
        [X,Y]=meshgrid(1:128,1:128);
        rois={'lesional',cx1,cy1,rad1; 'perilesional',cx2,cy2,rad2; 'contra',cx1c,cy1c,rad1};
        srcs={'uncorr',S.xform_datafluor; 'corr',S.xform_datafluorCorr};

        for ri=1:size(rois,1)
            nm=rois{ri,1}; cx=rois{ri,2}; cy=rois{ri,3}; rad=rois{ri,4};
            disk=((X-cx).^2+(Y-cy).^2)<=rad^2 & mask;
            if nnz(disk)==0, continue; end
            for ci=1:size(srcs,1)
                data=srcs{ci,2};
                tr=zeros(numel(idx),1);
                for kk=1:numel(idx)
                    sl=data(:,:,idx(kk)); tr(kk)=mean(sl(disk),'omitnan');
                end
                rowk=rowk+1;
                EXPORT(rowk).mouse=mnum; EXPORT(rowk).session=sess;
                EXPORT(rowk).roi=nm; EXPORT(rowk).corr=srcs{ci,1};
                EXPORT(rowk).fs=fs; EXPORT(rowk).trace=tr(:)';
            end
        end
        fprintf('Mouse %s [%s]: ok (fs=%.2f, %d samples)\n',mouseID,sess,fs,numel(idx));
        clear S
    end
end

outfile=fullfile(OUT_DIR,'ROI_Traces_export.mat');
save(outfile,'EXPORT','WINDOW_SEC','START_SEC','-v7');
fprintf('\nSaved %s  (%d rows)\n',outfile,numel(EXPORT));
disp('Done. Stage ROI_Traces_export.mat back so the figure can be built.');
