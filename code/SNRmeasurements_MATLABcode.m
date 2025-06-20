% Program modified from Peter Madsen's batchmag.m script and Dana Cusano's
% sei whale RL script.
%Written to use Mark Johnson's script for postemph.m to correct or
% compensate for the highpass filter on the Dtag in measurements of the
% signals. All signals have been decimated to 2kHz sampling rate to
% eliminate high frequency noise sources (target calls are < 100 Hz generally). 
%Additional modification includes pasting of noise sample prior to signal clip.

%Compiled to make RL measurements from fin whale tags by Dana Adcock 2025.

%Define directory in which your sound clips are
myDirAU='path'
%Gets all files in myDirAU that are wav files and stores names in myAUFIles
s=dir(fullfile(myDirAU, '*.wav'));

% Enter the correct calibration value here for individual Dtags plus any
% gain. Otherwise use the 178 value from Holt et al. 2017 calibration of D3's. 
% 184 for D4s (Tonneson et al., 2020, 178 for D3, 177.8 for CATS
cal= 184; 
%create empty data frames
magsrms=zeros(length(s),1);
magSNR=zeros(length(s),1);
noiseRL=zeros(length(s),1);

for k = 1:length(s)
    baseFileNameAU = s(k).name;
    fullFileNameAU = fullfile(myDirAU, baseFileNameAU);
    fprintf(1, 'Working on sound file %s\n', fullFileNameAU); 

    % Load audio file. Get data and sampling frequency
    [y, fs] = audioread(fullFileNameAU);

    % Check if x has enough rows
    if size(y, 1) < 20
        fprintf('Warning: Audio data too short to prime filter states.\n');
        continue;
    end

    j = y(1:(0.25 * 16000)); % Measure the noise
    ce = cumsum(y.^2); % accumulate energy in window
    k1 = min(find(ce > 0.05 * max(ce), 1, 'first'));
    k2 = min(find(ce > 0.95 * max(ce), 1, 'first'));

    if isempty(k1) || isempty(k2) || k1 >= k2 || k2 > length(y)
        fprintf('Warning: Invalid indices k1 or k2 for energy window.\n');
        continue;
    end

    y = y(k1:k2); % window containing 90% energy

    % 2. Received rms level, dB re. 1uPa (rms)
    rms = round(10 * log10(mean(y.^2))) + cal;

    % SNR measurement, dB re 1uPa (rms)
    noise = round(10 * log10(mean(j.^2))) + cal;
    SNR = rms - noise;

    % RLV = 20 * log10(std(y)); % this is in ‘dB re Volts’ of the postemphasis modified signal after the noise measurement segment.
    % rms = RLV + cal; % Subtracting the negative value of the tag sensitivity is equivalent to adding the positive cal value above.
    % SNR = 20 * log10(std(y) / std(n));
    magsrms(k) = rms;
    magSNR(k) = SNR;
    noiseRL(k) = noise; % modification added in March 2018, calculated noiseRL on the tag.
end

writematrix(magSNR, 'bp24_204a_SNR.txt');

writematrix(magsrms, 'bp24_204a_RMS.txt');

clear;

%wavSL