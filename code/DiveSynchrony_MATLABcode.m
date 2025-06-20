%%% Fin Whale Pair Synchrony
% Dana Adcock
% Last updated 06.16.2025
%% Load data for both whales
W1 = load("path");
W2 = load("path");

%% Time and depth setup
% Whale 1
start_time_1 = datetime('2024/07/25 09:30:08', 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
fs1 = W1.fs;
depth_1 = -W1.p;
end_time_1 = start_time_1 + seconds((length(depth_1) - 1) / fs1);
time_1 = start_time_1:seconds(1/fs1):end_time_1;

% Whale 2
start_time_2 = datetime('2024/07/25 09:40:06', 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
fs2 = W2.fs;
depth_2 = -W2.p;
end_time_2 = start_time_2 + seconds((length(depth_2) - 1) / fs2);
time_2 = start_time_2:seconds(1/fs2):end_time_2;

%% Align and interpolate to common time base
% Find overlapping time window
common_start = max(time_1(1), time_2(1));
common_end   = min(time_1(end), time_2(end));
if common_start >= common_end
    error('No overlapping time period between deployments.');
end
% Or manually add
common_end = datetime('2024-07-25 10:12:47', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Use the slower sampling rate to avoid oversampling
fs_common = min(fs1, fs2);
common_time = common_start:seconds(1/fs_common):common_end;

% Interpolate both depth signals to common time vector
depth_1_interp = interp1(time_1, depth_1, common_time, 'linear');
depth_2_interp = interp1(time_2, depth_2, common_time, 'linear');

% Remove NaNs introduced by interpolation
valid_idx = ~isnan(depth_1_interp) & ~isnan(depth_2_interp);
depth_1_interp = depth_1_interp(valid_idx);
depth_2_interp = depth_2_interp(valid_idx);
common_time = common_time(valid_idx);

%% Normalize signals (zero mean, unit variance)
norm_1 = (depth_1_interp - mean(depth_1_interp)) / std(depth_1_interp);
norm_2 = (depth_2_interp - mean(depth_2_interp)) / std(depth_2_interp);

%% Plot aligned and normalized depth profiles
figure;
plot(common_time, norm_1, '-', 'Color', '#440154', 'DisplayName', 'Whale 1');
hold on;
plot(common_time, norm_2, '-', 'Color', '#22A884', 'DisplayName', 'Whale 2');
xlabel('Time');
ylabel('Normalized Depth');
title('Male Pair Dive Profile');
legend;

%% Correlation without time shift (original alignment)
fixed_corr = corr(norm_1(:), norm_2(:));
fprintf('Correlation: %.3f\n', fixed_corr);
% Correlation: 0.696

%% Second pair
% Load data for both whales
W3 = load("path");
W4 = load("path");

%% Time and depth setup
% Whale 1
start_time_1 = datetime('2024/07/22 17:35:18', 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
fs1 = W3.fs;
depth_1 = -W3.p;
end_time_1 = start_time_1 + seconds((length(depth_1) - 1) / fs1);
time_1 = start_time_1:seconds(1/fs1):end_time_1;

% Whale 2
start_time_2 = datetime('2024/07/22 18:46:33', 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
fs2 = W4.fs;
depth_2 = -W4.p;
end_time_2 = start_time_2 + seconds((length(depth_2) - 1) / fs2);
time_2 = start_time_2:seconds(1/fs2):end_time_2;

%% Align and interpolate to common time base
% Find overlapping time window
common_start = max(time_1(1), time_2(1));
common_end   = min(time_1(end), time_2(end));
if common_start >= common_end
    error('No overlapping time period between deployments.');
end
% Or manually add
common_start = datetime('2024-07-22 18:56:56', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
common_end = datetime('2024-07-22 19:07:07', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Use the slower sampling rate to avoid oversampling
fs_common = min(fs1, fs2);
common_time = common_start:seconds(1/fs_common):common_end;

% Interpolate both depth signals to common time vector
depth_1_interp = interp1(time_1, depth_1, common_time, 'linear');
depth_2_interp = interp1(time_2, depth_2, common_time, 'linear');

% Remove NaNs introduced by interpolation
valid_idx = ~isnan(depth_1_interp) & ~isnan(depth_2_interp);
depth_1_interp = depth_1_interp(valid_idx);
depth_2_interp = depth_2_interp(valid_idx);
common_time = common_time(valid_idx);

%% Normalize signals (zero mean, unit variance)
norm_1 = (depth_1_interp - mean(depth_1_interp)) / std(depth_1_interp);
norm_2 = (depth_2_interp - mean(depth_2_interp)) / std(depth_2_interp);

%% Plot aligned and normalized depth profiles
figure;
plot(common_time, norm_1, '-', 'Color', '#440154', 'DisplayName', 'Whale 3');  % Blue
hold on;
plot(common_time, norm_2, '-', 'Color', '#22A884', 'DisplayName', 'Whale 4');  % Red
xlabel('Time');
ylabel('Normalized Depth');
title('Male Pair Dive Profile');
legend;

%% Correlation without time shift (original alignment)
fixed_corr = corr(norm_1(:), norm_2(:));
fprintf('Correlation without lag: %.3f\n', fixed_corr);
% Correlation without lag: 0.527

%% Third pair
% Load data for both whales
W5 = load("path");
W6 = load("path");

%% Time and depth setup
% Whale 1
start_time_1 = datetime('2023/07/24 13:02:53', 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
fs1 = W5.fs;
depth_1 = -W5.p;
end_time_1 = start_time_1 + seconds((length(depth_1) - 1) / fs1);
time_1 = start_time_1:seconds(1/fs1):end_time_1;

% Whale 2
start_time_2 = datetime('2023/07/24 16:49:20', 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
fs2 = W6.fs;
depth_2 = -W6.p;
end_time_2 = start_time_2 + seconds((length(depth_2) - 1) / fs2);
time_2 = start_time_2:seconds(1/fs2):end_time_2;

%% Align and interpolate to common time base
% Find overlapping time window
common_start = max(time_1(1), time_2(1));
common_end   = min(time_1(end), time_2(end));
if common_start >= common_end
    error('No overlapping time period between deployments.');
end
% Or manually add
common_end = datetime('2023-07-24 18:12:33', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Use the slower sampling rate to avoid oversampling
fs_common = min(fs1, fs2);
common_time = common_start:seconds(1/fs_common):common_end;

% Interpolate both depth signals to common time vector
depth_1_interp = interp1(time_1, depth_1, common_time, 'linear');
depth_2_interp = interp1(time_2, depth_2, common_time, 'linear');

% Remove NaNs introduced by interpolation
valid_idx = ~isnan(depth_1_interp) & ~isnan(depth_2_interp);
depth_1_interp = depth_1_interp(valid_idx);
depth_2_interp = depth_2_interp(valid_idx);
common_time = common_time(valid_idx);

%% Normalize signals (zero mean, unit variance)
norm_1 = (depth_1_interp - mean(depth_1_interp)) / std(depth_1_interp);
norm_2 = (depth_2_interp - mean(depth_2_interp)) / std(depth_2_interp);

%% Plot aligned and normalized depth profiles
figure;
plot(common_time, norm_1, '-', 'Color', '#440154', 'DisplayName', 'Whale 5');
hold on;
plot(common_time, norm_2, '-', 'Color', '#22A884', 'DisplayName', 'Whale 6');
xlabel('Time');
ylabel('Normalized Depth');
title('Male Pair Dive Profile');
legend;

%% Correlation without time shift (original alignment)
fixed_corr = corr(norm_1(:), norm_2(:));
fprintf('Correlation without lag: %.3f\n', fixed_corr);
% Correlation without lag: 0.530

%% Load data for both whales
W7 = load("G:\AS-Filer\BIO\sparks\Students\General\DanaA\1_Finwhale_SENE\Data\Tag data\Dtag\2023\bp23\prh\bp23_212e_prh50_auto.mat");
W8 = load("G:\AS-Filer\BIO\sparks\Students\General\DanaA\1_Finwhale_SENE\Data\Tag data\Dtag\2023\bp23\prh\archive\bp23_212fprh25_2.mat");

%% Time and depth setup
% Whale 7
start_time_1 = datetime('2023/07/31 15:22:52', 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
fs1 = W7.fs;
pitch_1 = -W7.pitch;
end_time_1 = start_time_1 + seconds((length(pitch_1) - 1) / fs1);
time_1 = start_time_1:seconds(1/fs1):end_time_1;

% Whale 8
start_time_2 = datetime('2023/07/31 17:21:25', 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
fs2 = W8.fs;
pitch_2 = -W8.pitch;
end_time_2 = start_time_2 + seconds((length(pitch_2) - 1) / fs2);
time_2 = start_time_2:seconds(1/fs2):end_time_2;

%% Align and interpolate to common time base
% Find overlapping time window
common_start = max(time_1(1), time_2(1));
common_end   = min(time_1(end), time_2(end));
if common_start >= common_end
    error('No overlapping time period between deployments.');
end
% Or manually add
common_end = datetime('2023-08-01 05:23:17', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Use the slower sampling rate to avoid oversampling
fs_common = min(fs1, fs2);
common_time = common_start:seconds(1/fs_common):common_end;

% Interpolate both depth signals to common time vector
pitch_1_interp = interp1(time_1, pitch_1, common_time, 'linear');
pitch_2_interp = interp1(time_2, pitch_2, common_time, 'linear');

% Remove NaNs introduced by interpolation
valid_idx = ~isnan(pitch_1_interp) & ~isnan(pitch_2_interp);
pitch_1_interp = pitch_1_interp(valid_idx);
pitch_2_interp = pitch_2_interp(valid_idx);
common_time = common_time(valid_idx);

%% Normalize signals (zero mean, unit variance)
norm_1 = (pitch_1_interp - mean(pitch_1_interp)) / std(pitch_1_interp);
norm_2 = (pitch_2_interp - mean(pitch_2_interp)) / std(pitch_2_interp);

%% Plot aligned and normalized depth profiles
figure;
plot(common_time, norm_1, '-', 'Color', '#440154', 'DisplayName', 'Whale 7');
hold on;
plot(common_time, norm_2, '-', 'Color', '#22A884', 'DisplayName', 'Whale 8');
xlabel('Time');
ylabel('Normalized Pitch');
title('Male-Female Pair Pitch');
legend;

%% Correlation without time shift (original alignment)
fixed_corr = corr(norm_1(:), norm_2(:));
fprintf('Correlation: %.3f\n', fixed_corr);
% Correlation: 0.003
