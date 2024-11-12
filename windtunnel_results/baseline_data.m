clear;
close all;
clc;
GraphGood();



% Get list of all txt files in the GEARBOX DRAG folder and sort them
files = dir('GEARBOX DRAG/*.txt');
% Extract numbers from filenames and sort numerically (descending order)
numbers = zeros(1, length(files));
for i = 1:length(files)
    [~, filename, ~] = fileparts(files(i).name);
    numbers(i) = str2double(regexp(filename, '\d+', 'match'));
end
[~, order] = sort(numbers, 'descend');
files = files(order);

% Initialize arrays to store data
col2_data = [];
legend_names = {};

% First pass: determine minimum length of all data files
min_length = inf;
for i = 1:length(files)
    filepath = fullfile(files(i).folder, files(i).name);
    data = readmatrix(filepath);
    min_length = min(min_length, size(data, 1));
end

% Second pass: read and store truncated data
for i = 1:length(files)
    % Read the current file
    filepath = fullfile(files(i).folder, files(i).name);
    data = readmatrix(filepath);
    
    % Store only column 2, inverted and filtered
    filtered_data = filter(ones(1,100)/100, 1, -data(1:min_length,2));
    col2_data(:,i) = filtered_data;
    
    % Store filename for legend, removing any non-numeric characters
    [~, filename, ~] = fileparts(files(i).name);
    % Convert filename to number for proper ordering
    num_only = str2double(regexp(filename, '\d+', 'match'));
    legend_names{i} = sprintf('%d m/s', num_only);
end

% Create single plot with updated style
figure('Units', 'centimeters', 'Position', [10, 10, 12, 10])
plot(col2_data, 'LineWidth', 1.5)
% title('Gearbox Drag', 'Interpreter', 'latex')
xlabel('Data Points', 'Interpreter', 'latex')
ylabel('Drag (N)', 'Interpreter', 'latex')
legend(legend_names, 'Location', 'southeast', 'Interpreter', 'latex')
grid on
grid minor
box on
xlim([0,9000])

% Save as SVG
saveas(gcf, 'gearbox_drag.svg')

% Calculate mean drag for each velocity
mean_drags = mean(col2_data);

% After plotting, calculate and display means
fprintf('\nMean Drag Values:\n');
fprintf('----------------\n');
for i = 1:length(col2_data(1,:))
    mean_drag = mean(col2_data(:,i));
    fprintf('Velocity: %s \t Mean Drag: %.3f N\n', legend_names{i}, mean_drag);
end


% Read membrane benchmark data from folder
folder_path = '.\MEMBRANE RESULTS BENCHMARK';
membrane_files = dir(fullfile(folder_path, '*.txt'));

% Get filter object (using the same filter as in inertial_data_reader.m)
Hd = INERTIAL_FILTER();

% Read files and process data
if length(membrane_files) >= 2
    % Read both files
    data1 = readmatrix(fullfile(folder_path, membrane_files(1).name));
    data2 = readmatrix(fullfile(folder_path, membrane_files(2).name));
    
    % Find minimum length
    min_length = min(size(data1,1), size(data2,1));
    
    % Extract drag (column 1) and lift (column 2) from both files
    drag1 = data1(1:min_length,1);
    lift1 = data1(1:min_length,2);
    drag2 = data2(1:min_length,1);
    lift2 = data2(1:min_length,2);
    
    % Apply filtering (five passes to make it stronger)
    filtered_drag1 = filtfilt(Hd.Numerator, 1, drag1);
    filtered_drag1 = filtfilt(Hd.Numerator, 1, filtered_drag1); % second pass
    filtered_drag1 = filtfilt(Hd.Numerator, 1, filtered_drag1); % third pass
    filtered_drag1 = filtfilt(Hd.Numerator, 1, filtered_drag1); % fourth pass
    filtered_drag1 = filtfilt(Hd.Numerator, 1, filtered_drag1); % fifth pass
    
    filtered_lift1 = filtfilt(Hd.Numerator, 1, lift1);
    filtered_lift1 = filtfilt(Hd.Numerator, 1, filtered_lift1); % second pass
    filtered_lift1 = filtfilt(Hd.Numerator, 1, filtered_lift1); % third pass
    filtered_lift1 = filtfilt(Hd.Numerator, 1, filtered_lift1); % fourth pass
    filtered_lift1 = filtfilt(Hd.Numerator, 1, filtered_lift1); % fifth pass
    
    filtered_drag2 = filtfilt(Hd.Numerator, 1, drag2);
    filtered_drag2 = filtfilt(Hd.Numerator, 1, filtered_drag2); % second pass
    filtered_drag2 = filtfilt(Hd.Numerator, 1, filtered_drag2); % third pass
    filtered_drag2 = filtfilt(Hd.Numerator, 1, filtered_drag2); % fourth pass
    filtered_drag2 = filtfilt(Hd.Numerator, 1, filtered_drag2); % fifth pass
    
    filtered_lift2 = filtfilt(Hd.Numerator, 1, lift2);
    filtered_lift2 = filtfilt(Hd.Numerator, 1, filtered_lift2); % second pass
    filtered_lift2 = filtfilt(Hd.Numerator, 1, filtered_lift2); % third pass
    filtered_lift2 = filtfilt(Hd.Numerator, 1, filtered_lift2); % fourth pass
    filtered_lift2 = filtfilt(Hd.Numerator, 1, filtered_lift2); % fifth pass
    
    % Create time vector
    time = (0:min_length-1) * 0.001; % Assuming 1000 Hz sampling rate
    
    % Extract frequencies from filenames
    [~, filename1, ~] = fileparts(membrane_files(1).name);
    [~, filename2, ~] = fileparts(membrane_files(2).name);
    
    freq1 = regexp(filename1, '(\d+\.?\d*)Hz', 'tokens');
    freq2 = regexp(filename2, '(\d+\.?\d*)Hz', 'tokens');
    
    legend1 = sprintf('%.1f Hz', str2double(freq1{1}{1}));
    legend2 = sprintf('%.1f Hz', str2double(freq2{1}{1}));
    
    % Create figure for drag comparison
    figure('Units', 'centimeters', 'Position', [10, 10, 11, 7])
    plot(time, filtered_drag1, 'LineWidth', 1.5)
    hold on
    plot(time, filtered_drag2, 'LineWidth', 1.5)
    xlabel('Time (s)', 'Interpreter', 'latex')
    ylabel('Drag Force (N)', 'Interpreter', 'latex')
    legend({legend1, legend2}, 'Location', 'best', 'Interpreter', 'latex')
    grid on
    grid minor
    box on
    xlim([0,2])
    ylim([-8,4])
    saveas(gcf, 'membrane_drag.svg')  % Save drag plot
    
    % Create figure for lift comparison
    figure('Units', 'centimeters', 'Position', [10, 10, 11, 7])
    plot(time, filtered_lift1, 'LineWidth', 1.5)
    hold on
    plot(time, filtered_lift2, 'LineWidth', 1.5)
    xlabel('Time (s)', 'Interpreter', 'latex')
    ylabel('Lift Force (N)', 'Interpreter', 'latex')
    legend({legend1, legend2}, 'Location', 'best', 'Interpreter', 'latex')
    grid on
    grid minor
    box on
    xlim([0,2])
    ylim([-8,4])
    saveas(gcf, 'membrane_lift.svg')  % Save lift plot
    
else
    fprintf('Need at least two files in membrane results folder\n');
end

