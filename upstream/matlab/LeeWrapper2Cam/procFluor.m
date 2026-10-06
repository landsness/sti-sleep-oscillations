function fluorData = procFluor(rawData,baseline)
%procFluor Preprocessing for fluorescence data
%   detrendBool = if true, detrends using info's highpass filter frequency
fluorData = double(rawData);
clear rawData
 if ~exist('baseline','var')
    baseline = nanmean(fluorData,length(size(fluorData)));
 else
     baseline = double(baseline);
 end

% --- DEBUG: Check baseline and ratio in procFluor ---
baseline_mean = mean(baseline(:));
fprintf('DEBUG (procFluor): Baseline mean = %g\n', baseline_mean);

% === INSERT FLASHLIGHT DEBUGGING HERE ===
disp('--- procFluor CRITICAL FLASHLIGHT: Before Division ---');
fprintf('Baseline min/max: %.4f / %.4f\n', min(baseline(:)), max(baseline(:)));
fprintf('Number baseline==0: %d\n', sum(baseline(:)==0));
fprintf('Number fluorData==0: %d\n', sum(fluorData(:)==0));
disp('------------------------------------------------------');


% Compute the ratio before subtracting 1
ratioData = fluorData ./ repmat(baseline, [1 1 1 size(fluorData,4)]);
ratio_mean = mean(ratioData(:));
fprintf('DEBUG (procFluor): Ratio mean (should be near 1) = %g\n', ratio_mean);
% --- End DEBUG ---


fluorData = fluorData./repmat(baseline,[1 1 1 size(fluorData,4)]); % make the data ratiometric


fluorData = fluorData - 1; % make the data change from baseline (center at zero)
%fluorData = mouse.process.smoothImage(fluorData,5,1.2); % spatially smooth data
fluorData = single(fluorData);
end

