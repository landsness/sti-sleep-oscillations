function [rawFig, displacement] = qcRaw(runInfo,time,raw,representativeCh,isbrain,system)
%qcRawFig Generates figure with quality control for raw data
% Inputs:
%   time = time values in seconds. Cell array with one cell for one channel
%   raw = raw data. Cell array with one cell for one channel
%   representativeCh = which channel should be used for representative
%   subplots
load(runInfo.saveMaskFile) %error!

for rch = 1:size(representativeCh)
    [displacement(rch,:)] = checkMovement(raw(:,:,representativeCh(rch),:));
end

% apply the mask to raw data
for i=1:size(raw,3)
    rawMask = raw(:,:,i,:).*isbrain;
    rawMask(rawMask==0) = NaN;
    raw(:,:,i,:) = rawMask;
end
rawSpatialAvg = nan(size(raw,4),size(raw,3));
fsTemp = diff(time);
T = fsTemp(1);
fs= 1./T;


for ch = 1:size(raw,3)
    rawSpatialAvg(:,ch) = squeeze(nanmean(nanmean(raw(:,:,ch,:),1),2));
    rawStd(ch) = std(rawSpatialAvg(:,ch));
    rawStdPer(ch) = 100*rawStd(ch)./mean(rawSpatialAvg(:,ch));
    rawTrend(:,ch) = bsxfun(@rdivide,rawSpatialAvg(:,ch),mean(rawSpatialAvg(:,ch),1));

       
%     rawTrend = double(bsxfun(@rdivide,rawSpatialAvg(:,ch),mean(rawSpatialAvg(:,ch))))';
    tNum(:,ch) = numel(rawTrend(:,ch));
    fdata_tmp(:,ch) = abs(fft(-log(rawTrend(:,ch)),size(rawTrend(:,ch),1)));
    hz_tmp(:,ch) = linspace(0,fs,tNum(:,ch));
    fdata(:,ch) = fdata_tmp(1:end/2,ch);
    hz(:,ch) =hz_tmp(1:end/2,ch);
end

% diffF2F = zeros(size(representativeCh),2);
% diff12F = zeros(size(representativeCh),2);

% difference from frame to frame
for rch = 1:size(representativeCh,2)
    diffF2Fraw = diff(raw(:,:,representativeCh(rch),:),[],4);
    diffF2Fmag = squeeze(nansum(nansum(abs(diffF2Fraw),1),2))';
    diffF2F(rch,:) = [0 diffF2Fmag];
    
    % difference from first frame to frame
    diff12Fraw = bsxfun(@minus,raw(:,:,representativeCh(rch),:),raw(:,:,representativeCh(rch),1)); %all frames minus the first frame
    diff12F(rch,:) = squeeze(squeeze(nansum(nansum(abs(diff12Fraw),1),2)))';
end

saveLoc = [runInfo.saveFolder filesep runInfo.recDate '-' runInfo.mouseName '-'...
runInfo.session '' num2str(runInfo.run) '-rawQCData.mat'];



%%Plotting
rawFig = figure('Position',[100 100 1000 700]);

numCh = size(rawSpatialAvg,2);
chStr = cell([1,numCh]);
for ch = 1:numel(chStr)
    chStr{ch} = ['Ch ' num2str(ch)];
end
%add title to subplot
mouseInfo = [runInfo.recDate, '-', runInfo.mouseName, '-', runInfo.session, char(string(runInfo.run))];
currVer = convertCharsToStrings(version('-release'));
if currVer == "2018a"
    suptitle(mouseInfo)
else
    sgtitle(mouseInfo)
end

% plot average abs time course
s1 = subplot('Position',[0.05 0.65 0.35 0.25]);
if numel(runInfo.fluorChInd) == 1
    if contains(lower(system),'gcamp')
        Colors = [0 0 1; 0 1 0; 1 0.5 0;];
        TickLabels = {'B','G','R' };
        legendName = {'Blue fluoro', 'green LED',  'red LED'}; 
    elseif contains(lower(system),'rgeco')
        Colors=[0 0 1; 0 1 0; 1 0 1; 1 0 0];
        TickLabels={'FAD', 'G', 'RGECO', 'R'};
        legendName = {'FAD', 'green LED', 'RGECO' 'red LED'};
    elseif strcmp(system, 'NewBiotech') %rmb added 241016
        Colors = [0 0 1; 0 1 0; 1 0.5 0;];
        TickLabels = {'B','G','R' };
        legendName = {'Blue fluoro', 'green LED',  'red LED'};   
    elseif strcmp(system, 'EastOIS1') %jpc220719 added
        Colors = [0 0 1; 0 1 0; 1 0.5 0;];
        TickLabels = {'B','G','R' };
        legendName = {'Blue fluoro', 'green LED',  'red LED'}; 
    end
elseif numel(runInfo.fluorChInd) == 2
    Colors=[0 0 1; 0 1 0; 1 0 1; 1 0 0];
    TickLabels={'GCaMP', 'G', 'RGECO', 'R'};
    legendName = {'GCaMP', 'green LED', 'RGECO' 'red LED'}; 
end

for ch = 1:numCh
    plot(s1,time,rawSpatialAvg(:,ch),'Color',Colors(ch,:));
    hold on;
end
title(s1,'Raw Data');
legend(legendName,'Location','northeastoutside');
xlabel('Time(sec)'); ylabel('Counts');

% plot average time trend
s2 = subplot('Position',[0.5 0.65 0.4 0.25]);
for ch = 1:numCh
    plot(s2,time,rawTrend(:,ch),'Color',Colors(ch,:));
    hold on;
end
title(s2,'Normalized Raw Data');
legend(legendName);
xlabel('Time(sec)');
ylabel('Normalized Counts');

% plot std dev
s3 = subplot('Position',[0.05 0.36 0.35 0.2]);
plot(s3,rawStdPer);
hold on;
scatter(s3,1:numel(rawStdPer),rawStdPer,'filled','k');
hold off;
ylabel('% Deviation');
title(s3,'Std Deviation');
xticks(1:numCh);
xticklabels(TickLabels)

% plot diff per frame
s4 = subplot('Position',[0.5 0.36 0.4 0.2]);
hold on;
xlabel('Time(sec)');
yyaxis(s4,'left');
plot(s4,time,diffF2F(1,:),'g');
hold on
if size(diffF2F,1) == 2
    plot(s4,time,diffF2F(2,:),'r');
    hold on
end
ylabel({'Sum Abs Diff Frame to Frame'; '(Counts)'});
yyaxis(s4,'right');
plot(s4,time,diff12F(1,:),'c');
hold on
if size(diff12F,1) == 2
    plot(s4,time,diff12F(2,:),'m');
    hold on
end
ylabel({'Sum Abs Diff Compared to First Frame'; '(Counts)'});
if numel(runInfo.fluorChInd) == 1
    legend({'Frame to Frame Change','Change Compared to First Frame'},'Location','northwest');
elseif numel(runInfo.fluorChInd) == 2
    legend({'GCaMP Frame to Frame Change','GCaMP Change Compared to First Frame','RGECO Frame to Frame Change','RGECO Change Compared to First Frame'},'Location','northwest');
end

% plot fft
s5 = subplot('Position',[0.05 0.06 0.35 0.2]);
for ch = 1:numCh
    loglog(s5,hz(:,ch),fdata(:,ch));
    hold on;
end
title('FFT of Log-meaned Raw Data');
xlim([min(hz(:,ch)) max(hz(:,ch))]);
legend(legendName,'Location','southwest');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% plot offset
s6 = subplot('Position',[0.5 0.06 0.4 0.2]);
plot(s6,time,displacement); title('Absolute offset in pixels');
xlabel('Time (sec)'); ylabel('Offset (pixels)');



save(saveLoc, 'runInfo','time','rawSpatialAvg','time','diffF2F',...
'diff12F','displacement','hz','fdata','fs','rawStd');
end


