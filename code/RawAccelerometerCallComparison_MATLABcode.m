%%% Raw Accelerometer Signals for Focal Call Identification
%%% Dana Adcock
%%% 05/22/2025

%% Add a filter to get out noise - only for D4
% load in prh
load('path')

original_rate = 5000;  % Original sampling rate in Hz
target_rate = 250;     % Desired sampling rate in Hz
decimation_factor = original_rate / target_rate;  % Decimation factor

% Design a low-pass filter to avoid aliasing
nyquist_rate = target_rate / 2;
[b, a] = butter(4, nyquist_rate / (original_rate / 2), 'low');

% Apply the low-pass filter to each axis of the accelerometer data
data = double(A.data);
data(isnan(data)) = 0;
filtered_data = filtfilt(b, a, data);

% Downsample the filtered data
downsampled_data = downsample(filtered_data, decimation_factor);

% Plot the original and downsampled data for comparison
t_original = (0:length(data)-1) / original_rate;
t_downsampled = (0:length(downsampled_data)-1) / target_rate;

%% Convert from gravity units to m/s^2
Ag = A*9.81 %D3
Ag = downsampled_data*9.81 %D4

%% Perform mean subtraction - Goldbogen 2014

% Calculate the mean of each column, ignoring NaNs
meanAg = nanmean(Ag);

% Subtract the mean from each column, preserving NaNs
Ag_mean_subtracted = bsxfun(@minus, Ag, meanAg);


%% Remove linear trends - Goldbogen 2014

t = 1:size(Ag_mean_subtracted, 1);  

chunkSize = 10;  % Adjust as needed
numChunks = ceil(size(Ag_mean_subtracted, 1) / chunkSize);

Ag_detrended = zeros(size(Ag_mean_subtracted));  % Initialize

for chunk = 1:numChunks
    idxStart = (chunk - 1) * chunkSize + 1;
    idxEnd = min(chunk * chunkSize, size(Ag_mean_subtracted, 1));
    
    % Process chunk of data
    Ag_chunk = Ag_mean_subtracted(idxStart:idxEnd, :);
    t_chunk = t(idxStart:idxEnd)';
    
    % Remove linear trends
    for i = 1:3
        % Fit a linear trend (line) to each axis (column) of Ag_chunk
        o = polyfit(t_chunk, Ag_chunk(:, i), 1);  % Fit a first-degree polynomial (linear) to the data
        
        % Create the linear trend for this axis
        trend = polyval(o, t_chunk);
        
        % Remove the linear trend from this axis within the chunk
        Ag_detrended(idxStart:idxEnd, i) = Ag_chunk(:, i) - trend;
    end
end


%% plot each Acc axis separately
plott(-p, 50, Ag_detrended(:,1),2.5,Ag_detrended(:,2), 2.5,Ag_detrended(:,3), 2.5); %dtag
plott(Ag_detrended(:,1),400,Ag_detrended(:,2),400,Ag_detrended(:,3),400); %cats

