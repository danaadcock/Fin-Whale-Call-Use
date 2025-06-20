%%%% Extract Dives
% Dana Adcock
% Last updated 05.22.2025
% Load prh
%load prh
clear;
file_path = 'path';
load(file_path);

%find dives
dives = find_dives(p2,fs,2);

%%
% Plot the dive profile
% Put depth in seconds
p_downsampled = downsample(p, fs);

% Convert sample rate to time in minutes
timeInMinutes = (0:length(p)-1) / fs / 3600; % Time in hours

figure;
plot(timeInMinutes, -p, 'k-', 'LineWidth', 1.5); % Plot -p
hold on;

% Mark the dive times on the graph
for i = 1:length(dives.start)
    % Convert call times from seconds to the corresponding index
    index = round(dives.start(i) * fs) + 1; % +1 for MATLAB's 1-based indexing
    if index <= length(p) % Ensure index is within bounds
        plot(timeInMinutes(index), -p(index), 'bo', 'MarkerSize', 8, 'LineWidth', 2); % Red circles
    end
end

for i = 1:length(dives.end)
    % Convert call times from seconds to the corresponding index
    index = round(dives.end(i) * fs) + 1; % +1 for MATLAB's 1-based indexing
    if index <= length(p) % Ensure index is within bounds
        plot(timeInMinutes(index), -p(index), 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Red circles
    end
end

% Add labels and title
xlabel('Time (hr)');
ylabel('Depth (m)');
title('Dives');
% grid on;
legend('Dive Profile', 'Dive Start', 'Dive End');

% Adjust axis limits if necessary
xlim([0 max(timeInMinutes)]);
ylim([min(-p) max(-p)]);
hold off;

%% remove short dives (surfacings)
% Remove observations with dive duration < 30
valid_indices = (dives.end - dives.start) >= 30;

% Apply indexing to all doubles in the struct
fields = fieldnames(dives);
for i = 1:numel(fields)
    if isnumeric(dives.(fields{i})) && ismatrix(dives.(fields{i}))
        dives.(fields{i}) = dives.(fields{i})(valid_indices);
    end
end

%% Add missed dives if needed
% Manually add additional first values to specified fields
dives.start = [dives.start; 377.52]; % Replace new_start_value with the desired value
dives.end = [dives.end; 495.96];       % Replace new_end_value with the desired value
dives.max = [dives.max; 58.514];       % Replace new_max_value with the desired value
dives.tmax = [dives.tmax; 391.72];    % Replace new_tmax_value with the desired value

%%
%save prh
save(file_path, 'dives','-append');

