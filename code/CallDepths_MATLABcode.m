%%%% Call Depths
% Dana Adcock
% Last updated 05.22.2025
% Load prh
clear;
load('path');

% Load call info
calls = readtable('path','Delimiter', '\t');

% Sample dive profile data (example)
beginTimes = calls.BeginTime_s_;

% Ensure indices are within bounds
callIndices = round(beginTimes * fs) + 1; % Convert seconds to indices

% Remove any out-of-bound indices
validMask = (callIndices <= length(p)) & (callIndices > 0);
validIndices = callIndices(validMask);

% Extract corresponding dive depths
callDepths = nan(height(calls), 1); % Preallocate with NaNs
callDepths(validMask) = -p(validIndices); % Assign valid depth values

% Add dive depth column to the calls table
calls.Depth_m = callDepths;

% Save the updated table to a CSV file
writetable(calls, 'path');

%%
%Merge csvs
% Define the folder containing the CSV files
folderPath = 'path';

% Get a list of all CSV files in the folder
csvFiles = dir(fullfile(folderPath, '*.csv'));

% Initialize an empty cell array to store tables
tables = cell(length(csvFiles), 1);
variableSets = cell(length(csvFiles), 1);

% Loop through each file and read the data
for i = 1:length(csvFiles)
    % Get full file path
    filePath = fullfile(folderPath, csvFiles(i).name);
    
    % Read the CSV file into a table
    tempTable = readtable(filePath);
    
    % Add a new column for the source file name
    tempTable.SourceFile = repmat({csvFiles(i).name}, height(tempTable), 1);
    
    % Store table and its variable names
    tables{i} = tempTable;
    variableSets{i} = tempTable.Properties.VariableNames;
end

% Start by assuming all columns from the first table are common
commonVarNames = variableSets{1};

% Iterate through the rest of the tables to find common columns
for i = 2:length(variableSets)
    commonVarNames = intersect(commonVarNames, variableSets{i});
end

% Filter each table to keep only the common columns, including the new SourceFile column
for i = 1:length(tables)
    tables{i} = tables{i}(:, [commonVarNames, 'SourceFile']); % Include 'SourceFile' column
end

% Concatenate all tables
mergedTable = vertcat(tables{:});

% Save the merged table to a new CSV file
outputFilePath = fullfile(folderPath, 'merged_calls_with_dive_depths.csv');
writetable(mergedTable, outputFilePath);

% Display the first few rows
disp(mergedTable);

%%
%Total time in depth bins
clear;

% Prompt user to select the source .mat file
[file_name, file_path] = uigetfile('*.mat', 'Select the Source File');
if isequal(file_name, 0)
    disp('No file selected. Exiting.');
    return;
end

% Load the .mat file
source_file = fullfile(file_path, file_name);
load(source_file);

% Depth Data - change based on tag type
depth_data = p;

% Define the depth bins (in meters)
depth_bins = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100];

%Plot depth
plot(-depth_data);

% Manually specify start and end time in samples
%LEAVE OUT FIRST TWO DIVES
start_sample = 13308;
end_sample = 367093; 

% Extract the relevant depth data
depth_segment = depth_data(start_sample:end_sample);

% Initialize an array to hold the time spent in each depth bin
time_in_bins = zeros(1, length(depth_bins)-1);

% Loop through each depth bin and calculate the time spent in each
for i = 1:length(depth_bins)-1
    % Find the samples within the current depth bin
    bin_samples = (depth_segment >= depth_bins(i)) & (depth_segment < depth_bins(i+1));
    
    % Calculate the time spent in this depth bin (in seconds)
    time_in_bins(i) = sum(bin_samples) / fs; % Convert samples to time
end

% Create a table to hold the results
result_table = table(depth_bins(1:end-1)', depth_bins(2:end)', time_in_bins', 'VariableNames', {'Depth_Lower', 'Depth_Upper', 'Time_Spent_sec'});

% Add the source file name to the table
result_table.Source_File = repmat({file_name}, size(result_table, 1), 1);

% Export the results to a CSV file
output_path = 'path'
output_file = fullfile(output_path, 'name'); % change each time
writetable(result_table, output_file);

disp(['Export complete! Results saved to: ', output_file]);
