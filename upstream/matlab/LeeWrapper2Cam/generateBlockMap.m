function [fh,peakMaps] = generateBlockMap(data,runInfo,numBlocks,peakRange,xform_isbrain)

% pres of each block during stim period
parampath=what('bauerparams');
load(which('noVasculatureMask.mat')) 
mask_new = logical(mask_new);
numRows = length(runInfo.Contrasts);
fh = figure('units','normalized','outerposition',[0 0 1 1]);
colormap jet
peakMaps = nan(128,128,numBlocks+1,numRows);    %JPC added 220805  
for jj = 1:size(data,5) %changed from length(runInfo.Contrasts) 221031
     peakMap = squeeze(mean(data(:,:,:,:,jj),4,'omitnan'));
     peakMap = squeeze(mean(peakMap(:,:,(round(peakRange{jj}(1)*runInfo.samplingRate)):(round(peakRange{jj}(end)*runInfo.samplingRate))),3,'omitnan'));
     maxVal = prctile(abs(peakMap(mask_new)),90,'all')*2.5;
     maxVal = round(maxVal,2,'significant');
    for ii = 1:numBlocks+1
        p = subplot(numRows,numBlocks+1,(numBlocks+1)*(jj-1)+ii);
        if ii == numBlocks+1
            imagesc(peakMap,'AlphaData',xform_isbrain)
            peakMaps(:,:,ii,jj) = peakMap;            %JPC added 220805     
            axis image
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
            title('Averaged')
            Pos = get(p,'Position');
            cb = colorbar;
            try
            caxis([-maxVal maxVal])
            end
            set(p,'Position',Pos)
            try
            set(cb,'YTick',[-maxVal,0,maxVal]);
            end

            if sum(contains({'HbO','HbR','HbT'},runInfo.Contrasts{jj}))>0
            set(get(cb,'label'),'string','Hb(\Delta\muM)');
            elseif sum(contains({'Calcium','FAD'},runInfo.Contrasts{jj}))>0
                set(get(cb,'label'),'string','Fluorescence(\DeltaF/F%)');
            end
        else
            pre = squeeze(mean(data(:,:,(round(peakRange{jj}(1)*runInfo.samplingRate)):(round(peakRange{jj}(end)*runInfo.samplingRate)),ii,jj),3,'omitnan'));
            imagesc(pre,'AlphaData',xform_isbrain);
            peakMaps(:,:,ii,jj) = pre;    %JPC added 220805  
            
            cb = colorbar;
            try
            caxis([-maxVal maxVal])
            end
            try
            set(cb,'YTick',[-maxVal,0,maxVal]);
            end            
            axis image
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
            if jj==1
                title(strcat('Pres',{' '},num2str(ii)))
            end
            if ii == 1
                ylabel(runInfo.Contrasts{jj})
            end
        end
    end
end
