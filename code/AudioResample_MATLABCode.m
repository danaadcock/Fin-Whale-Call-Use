%% Downsampling to 2 kHz
clear;

% Define the folder containing the WAV files and the target sample rate
folderPath = 'path'; % Adjust this to your folder path
targetSampleRate = 2000; 

% Define the base name and the number of files to process
baseName = 'TagID'; % change to tag ID 
numFiles = 6; % Adjust this to the number of files you have

% Construct the first file name to check its sample rate
fileName = sprintf('%s%03d.wav', baseName, 1);
d4FilePath = fullfile(folderPath, fileName);

% Read the WAV file to get the sample rate
[~, originalSampleRate] = audioread(d4FilePath);
fprintf('Original sample rate of %s: %d Hz\n', fileName, originalSampleRate);

% Loop through each file
for i = 1:numFiles
    % Construct the file name with the sequential number (e.g., 001, 002)
    fileName = sprintf('%s00%d.wav', baseName, i);
    d4FilePath = fullfile(folderPath, fileName);
    
    % Read the WAV file
    [audioData, originalSampleRate] = audioread(d4FilePath);
    
    % Call postemph function
    try
        [y, fs, b, a] = postemph(audioData, originalSampleRate);
    catch ME
        fprintf('Error in postemph function: %s\n', ME.message);
        continue;
    end
    
    % Print durations to verify
    durOrig = length(audioData) / originalSampleRate;
    durPost = length(y) / fs;
    fprintf('Durations (s): Original = %.2f | Postemph = %.2f\n', durOrig, durPost);
    
    % Design a low-pass filter (cutoff slightly below Nyquist of target rate)
    [c, d] = butter(6, 900 / (fs / 2), 'low');
    audfil = filtfilt(c, d, y);
    
    % Compute resample ratio based on fs from postemph
    [P, Q] = rat(targetSampleRate / fs, 1e-6);
    
    % Resample to 2 kHz
    resampledAudio = resample(audfil, P, Q);
    
    % Define the output file name
    outputFileName = strrep(fileName, '.wav', '_2kHz.wav');
    outputFilePath = fullfile(folderPath, outputFileName);
    
    % Save the resampled and normalized audio
    audiowrite(outputFilePath, resampledAudio, targetSampleRate);
    
    fprintf('Resampled file saved as: %s\n', outputFilePath);
end
