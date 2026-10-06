%% ===============================================================
%  SO_Artifact_BandCheck.m
%  Test whether the PT "bright GCaMP artifact" lives OUTSIDE the SO band.
%
%  From the UNCORRECTED GCaMP (xform_datafluor), compute frequency-band
%  power topoplots:
%     - LOW  (0 - 0.1 Hz)  : where a slow/bright artifact should live
%     - SO   (0.1 - 1 Hz)  : the signal of interest
%     - FULL (0 - fs/2)    : whole-band reference
%  plus a brightness proxy and the power spectrum at the lesion core vs the
%  contralateral homolog (the clinching plot).
%
%  If the artifact blob appears in LOW/FULL but NOT in SO -> the bright
%  artifact cannot contaminate the SO measurements.
%
%  Uses the lab's PowerAnalysis.m. Loads the big dataFluor.mat only for the
%  mice listed. Outputs are small -> easy to share.
% ===============================================================
clear; clc; close all;

%% ---------------- USER INPUT ----------------
BASE_DIR  = 'C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\PreProcessed_Data';
OUT_DIR   = 'C:\Users\landsness\Box\Manuscripts\STI_SO\Revisions\SO_ContrastMaps';
ROI_XLSX  = 'C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\SO_Analysis\ROI_Summary_V2.xlsx';
MICE      = {'11','13'};     % representative large-infarct mice; add more as desired
addpath('C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\SO_Analysis\MATLAB_Code_Used\LeeWrapper2Cam Scripts');

SO_LOW=0.1; SO_HIGH=1.0; LOW_HI=0.1; CORER=4;   % core disk radius (px)
sessions = {'bsl',{'bsl'}; 'acute',{'acute'}; 'wk1',{'oneweek','weekone','wk1','week1'}};
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

% ROI1 centroids
T = readtable(ROI_XLSX);
getC = @(mnum) deal(T.ROI1_CentroidX(T.Mouse==mnum), T.ROI1_CentroidY(T.Mouse==mnum));

for mi = 1:numel(MICE)
    m = MICE{mi}; mnum = str2double(m);
    [cx,cy] = getC(mnum);
    if isempty(cx), warning('no ROI1 centroid for %s',m); cx=64;cy=64; else, cx=cx(1);cy=cy(1); end
    fprintf('\n== Mouse %s  (ROI1 @ %.0f,%.0f) ==\n', m, cx, cy);
    R = struct();
    for si=1:size(sessions,1)
        sess=sessions{si,1}; folder='';
        for f=sessions{si,2}, c=fullfile(BASE_DIR,m,f{1}); if exist(c,'dir'),folder=c;break; end, end
        if isempty(folder), continue; end
        ff=dir(fullfile(folder,'*dataFluor.mat')); if isempty(ff), fprintf('  [%s] no dataFluor\n',sess); continue; end
        mf=dir(fullfile(folder,'*LandmarksAndMask.mat'));
        S=load(fullfile(folder,ff(1).name),'xform_datafluor','runInfo'); fs=S.runInfo.samplingRate;
        mask=true(128,128);
        if ~isempty(mf), M=load(fullfile(folder,mf(1).name),'xform_isbrain');
            if isfield(M,'xform_isbrain'), mask=~isnan(M.xform_isbrain)&(M.xform_isbrain>0.5); end, end

        [wsm,~,~,hz] = PowerAnalysis(double(S.xform_datafluor), fs, mask);
        bLow = hz>0 & hz<LOW_HI;  bSO = hz>SO_LOW & hz<SO_HIGH;
        mLow=mean(wsm(:,:,bLow),3); mSO=mean(wsm(:,:,bSO),3); mFull=mean(wsm,3);
        mLow(~mask)=NaN; mSO(~mask)=NaN; mFull(~mask)=NaN;
        bright=mean(abs(double(S.xform_datafluor)),3,'omitnan'); bright(~mask)=NaN;

        % core (ROI1 disk) and contralateral-mirror spectra
        [X,Y]=meshgrid(1:128,1:128);
        core=((X-cx).^2+(Y-cy).^2)<=CORER^2 & mask;
        contra=((X-(129-cx)).^2+(Y-cy).^2)<=CORER^2 & mask;
        spec_core   = squeeze(mean(mean(wsm.*core,1,'omitnan'),2,'omitnan'))/mean(core(:));
        spec_contra = squeeze(mean(mean(wsm.*contra,1,'omitnan'),2,'omitnan'))/mean(contra(:));

        R.(sess)=struct('mLow',mLow,'mSO',mSO,'mFull',mFull,'bright',bright, ...
                        'spec_core',spec_core,'spec_contra',spec_contra,'hz',hz, ...
                        'cx',cx,'cy',cy,'fs',fs);
        fprintf('  [%s] fs=%.2f Hz done\n',sess,fs);
        clear S wsm
    end
    save(fullfile(OUT_DIR,sprintf('Mouse_%s_ArtifactBandCheck.mat',m)),'R');
    fprintf('  saved Mouse_%s_ArtifactBandCheck.mat\n',m);
end
disp('Done. Stage the *_ArtifactBandCheck.mat files back for the figure.');
