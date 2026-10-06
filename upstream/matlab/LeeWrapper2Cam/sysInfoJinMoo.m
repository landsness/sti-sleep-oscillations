function info = sysInfoJinMoo(systemType)
%session2procInfo Making OIS preprocess info from session type
%   Input:
%       sessiontype = char array showing type of systemtype ('fcOIS1','fcOIS2')
%   Output:
%       info = struct with info such as rgb indices
%           rgb = 1x3 vector specifying indices for red, green, and blue
%           numLEDs = number of leds
%           LEDFiles = string array containing name of text files showing
%           LED spectra
%           readFcn = function handle for reading from raw file (tiff, dat)
%           invalidFrameInd = any temporal frame index that should be
%           removed prior to processing (these are any indices that are not
%           dark frames yet still need to be removed)
%           gbox = gaussian filter box size for smoothing image
%           gsigma = gaussian filter sigma for smoothing image

switch lower(systemType)
    case {'lee1_gcamp','gcamp'} %jpc220719 added
        info.rgb = [3 2 NaN];
        info.numLEDs = 3;
        info.binFactor = [4 4 4];
        info.LEDFiles = [
            "473Excite_bgSubtractJM.txt",...
            "532OIS_bgsubtractJM.txt", ...
            "625OIS_bgsubtractJM.txt"];
        info.fluorFiles = {'gcamp6f_emission.txt'};
        info.FADFiles = [];
        info.chHb = [2 3];
        info.chFluor = 1; %this means channel calcium!!
        info.chFAD = [];
        info.chSpeckle = [];
        info.chLaser = [];
        info.numInvalidFrames = 1;
        info.validThr = [1 inf];
        info.gbox = 5;
        info.gsigma = 1.2;
        info.Contrasts={'HbO','HbR','HbT','Calcium'};
        info.colors = {'r','b','k','g'};
        info.cam1Chan = [1 2]; %Only sees green
        info.cam2Chan = 3;
    case {'lee1_fad_rgeco','lee1_rgeco','rgeco'}
        info.rgb = [4 2 NaN];
        info.numLEDs = 4;
        info.binFactor = [4 4 4 4];
        info.LEDFiles = [
            "473Excite_bgSubtractJM.txt",...
            "532OIS_bgsubtractJM.txt", ...
            "532Excite_bgsubtractJM.txt",...
            "625OIS_bgsubtractJM.txt"];
        info.fluorFiles = {'jrgeco1a_emission.txt'};
        info.FADFiles = {'fad_emission.txt'};
        info.chHb = [2 4];
        info.chFluor = 3; %this means channel calcium!!
        info.chFAD = 1;
        info.chSpeckle = [];
        info.chLaser = [];
        info.numInvalidFrames = 1;
        info.validThr = [1 inf];
        info.gbox = 5;
        info.gsigma = 1.2;
        info.Contrasts={'HbO','HbR','HbT','Calcium','FAD'};
        info.colors = {'r','b','k','m','g'};
        info.cam1Chan = [1 2];
        info.cam2Chan = [3 4];
    case {'both','lee1_gcamp_rgeco','gcamp and rgeco','rgeco and gcamp'}
        info.rgb = [4 2 NaN];
        info.numLEDs = 4;
        info.binFactor = [4 4 4 4];
        info.LEDFiles = [
            "473Excite_bgSubtractJM.txt",...
            "532OIS_bgsubtractJM.txt", ...
            "532Excite_bgsubtractJM.txt",...
            "625OIS_bgsubtractJM.txt"];
        info.fluorFiles = {'gcamp6f_emission.txt','jrgeco1a_emission.txt'};
        info.FADFiles = [];
        info.chHb = [2 4];
        info.chFluor = [1 3]; %this means channel calcium!!
        info.chFAD = [];
        info.chSpeckle = [];
        info.chLaser = [];
        info.numInvalidFrames = 1;
        info.validThr = [1 inf];
        info.gbox = 5;
        info.gsigma = 1.2;
        info.Contrasts={'HbO','HbR','HbT','GCaMP','RGECO'};
        info.colors = {'r','b','k','m','g'};
        info.cam1Chan = [1 2];
        info.cam2Chan = [3 4];
end

paramPath = what('bauerParams'); % path to bauerParams module
if numel(paramPath)>1
    paramPath=paramPath(1);
end
sourceSpectraLoc = fullfile(paramPath.path,'ledSpectra'); % path to led spectra text files
for ledInd = 1:numel(info.LEDFiles)
    info.LEDFiles{ledInd} = fullfile(sourceSpectraLoc,info.LEDFiles{ledInd});
end
fluorSpectraLoc = fullfile(paramPath.path,'probeSpectra'); % path to fluor spectra text files
for chInd = 1:numel(info.fluorFiles)
    info.fluorFiles{chInd} = fullfile(fluorSpectraLoc,info.fluorFiles{chInd});
end


if ~isempty(info.FADFiles)
    for chInd = 1:numel(info.FADFiles)
        info.FADFiles{chInd} = fullfile(fluorSpectraLoc,info.FADFiles{chInd});
    end
end

