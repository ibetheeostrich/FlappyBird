clear; clc; close all;
GraphGood();

folder_path = '.\INERTIAL (CLEAN)';

% Get list of CSV files
files = dir(fullfile(folder_path, '*.csv'));

% Sort files by number in filename
numbers = zeros(length(files), 1);
for i = 1:length(files)
    filename = files(i).name;
    matches = regexp(filename, '\d+\.?\d*', 'match');
    numbers(i) = str2double(matches{1});
end
[~, order] = sort(numbers);
files = files(order);

% Initialize arrays
drag_data = [];
lift_data = [];
legend_names = {};

% First pass to find minimum length
min_length = inf;
for i = 1:length(files)
    filepath = fullfile(files(i).folder, files(i).name);
    data = readmatrix(filepath);
    min_length = min(min_length, size(data, 1));
end

% Read and store data
for i = 1:length(files)
    % Read file
    filepath = fullfile(files(i).folder, files(i).name);
    data = readmatrix(filepath);
    
    % Store truncated drag and lift data
    drag_data(:,i) = data(1:min_length,1);
    lift_data(:,i) = data(1:min_length,2);
    
    % Create legend entry
    [~, filename, ~] = fileparts(files(i).name);
    matches = regexp(filename, '\d+\.?\d*', 'match');
    num_only = str2double(matches{1});
    legend_names{i} = sprintf('%.1f Hz', num_only);
end

% After loading all data but before plotting
% Find the index of 1.5 Hz data
target_freq = 1.5;
file_freqs = zeros(length(files), 1);
for i = 1:length(files)
    matches = regexp(files(i).name, '\d+\.?\d*', 'match');
    file_freqs(i) = str2double(matches{1});
end
freq_idx = find(abs(file_freqs - target_freq) < 0.01);
% Shift the drag data for 1.5 Hz down by 4
drag_data(:,freq_idx) = drag_data(:,freq_idx) - 4;

% After loading all data but before plotting
% Find the first peak for each frequency and align them
first_peaks = zeros(length(files), 1);
for i = 1:length(files)
    % Find first peak in drag data with more strict criteria
    [~, peak_idx] = findpeaks(drag_data(:,i), ...
        'NPeaks', 1, ...
        'SortStr', 'descend', ...
        'MinPeakHeight', mean(drag_data(:,i)), ...
        'MinPeakProminence', 0.5);  % Adjust this value as needed
    
    if ~isempty(peak_idx)
        first_peaks(i) = peak_idx(1);
    end
end

% Find the minimum peak location to use as reference
ref_peak = min(first_peaks);

% Calculate how many samples we need for 1 second of data
samples_per_second = 1000; % Since sampling rate is 1000 Hz
desired_length = samples_per_second;

% Create new arrays for the aligned data
aligned_drag = zeros(desired_length, length(files));
aligned_lift = zeros(desired_length, length(files));

% Trim and align all signals
for i = 1:length(files)
    shift_amount = first_peaks(i) - ref_peak;
    
    % Create indices for the section we want to keep
    start_idx = shift_amount + 1;
    end_idx = start_idx + desired_length - 1;
    
    % Make sure we have enough data after the peak
    if end_idx > size(drag_data, 1)
        % Shift the window back to get a full second of data
        end_idx = size(drag_data, 1);
        start_idx = end_idx - desired_length + 1;
    end
    
    % Extract the aligned section
    aligned_drag(:,i) = drag_data(start_idx:end_idx, i);
    aligned_lift(:,i) = lift_data(start_idx:end_idx, i);
end

% Update time vector for the trimmed data
time = (0:0.001:(desired_length-1)*0.001)';

% Find indices of frequencies we want to display (exclude 1 Hz)
display_idx = file_freqs > 1;

% Find common y-axis limits
y_min = min(min([aligned_drag(:,display_idx), aligned_lift(:,display_idx)]));
y_max = max(max([aligned_drag(:,display_idx), aligned_lift(:,display_idx)]));

% Plot drag data
figure('Units', 'centimeters', 'Position', [10, 10, 12, 7])
plot(time, aligned_drag(:,display_idx), 'LineWidth', 1.5)
xlabel('Time (s)', 'Interpreter', 'latex')
ylabel('Drag Axis Inertial Force (N)', 'Interpreter', 'latex')
legend(legend_names(display_idx), 'Location', 'best', 'Interpreter', 'latex', ...
    'NumColumns', 4)  % Changed to 4 columns
xlim([0 1])
ylim([-6 10])
grid on
grid minor
box on
saveas(gcf, 'inertial_drag.svg')

% Plot lift data
figure('Units', 'centimeters', 'Position', [10, 10, 12, 7])
plot(time, aligned_lift(:,display_idx), 'LineWidth', 1.5)
xlabel('Time (s)', 'Interpreter', 'latex')
ylabel('Lift Axis Inertial Force (N)', 'Interpreter', 'latex')
legend(legend_names(display_idx), 'Location', 'best', 'Interpreter', 'latex', ...
    'NumColumns', 4)  % Changed to 4 columns
xlim([0 1])
ylim([-6 10])
grid on
grid minor
box on
saveas(gcf, 'inertial_lift.svg')

% %% Section 2: Period Analysis and Cycle Averaging
% % Find indices of frequencies we want to display (exclude 1 Hz)
% display_idx = file_freqs > 1;

% % Initialize arrays to store averaged cycles
% samples_per_cycle = 1000;  % We'll interpolate to this many points per cycle
% avg_drag_cycles = zeros(samples_per_cycle, sum(display_idx));
% avg_lift_cycles = zeros(samples_per_cycle, sum(display_idx));

% % Process each frequency
% for i = 1:length(files)
%     if ~display_idx(i)
%         continue;  % Skip 1 Hz data
%     end
    
%     % Find peaks for period calculation
%     [~, peak_locs] = findpeaks(aligned_drag(:,i), 'MinPeakDistance', 100);
    
%     % Calculate periods between peaks
%     periods = diff(peak_locs);
%     mean_period = mean(periods);
    
%     % Extract and average complete cycles
%     cycles_drag = [];
%     cycles_lift = [];
    
%     for j = 1:(length(peak_locs)-1)
%         % Extract one complete cycle
%         cycle_data_drag = aligned_drag(peak_locs(j):peak_locs(j+1)-1, i);
%         cycle_data_lift = aligned_lift(peak_locs(j):peak_locs(j+1)-1, i);
        
%         % Interpolate to common length
%         t_original = linspace(0, 1, length(cycle_data_drag));
%         t_interp = linspace(0, 1, samples_per_cycle);
        
%         cycles_drag(:,j) = interp1(t_original, cycle_data_drag, t_interp);
%         cycles_lift(:,j) = interp1(t_original, cycle_data_lift, t_interp);
%     end
    
%     % Average the cycles
%     avg_drag_cycles(:,i) = mean(cycles_drag, 2);
%     avg_lift_cycles(:,i) = mean(cycles_lift, 2);
% end

% % Create time vector for one cycle
% t_cycle = linspace(0, 1, samples_per_cycle)';

% % Plot averaged cycles
% figure('Units', 'centimeters', 'Position', [10, 10, 12, 10])
% plot(t_cycle, avg_drag_cycles(:,display_idx), 'LineWidth', 1.5)
% xlabel('Normalized Time (Period)', 'Interpreter', 'latex')
% ylabel('Drag Axis Inertial Force (N)', 'Interpreter', 'latex')
% legend(legend_names(display_idx), 'Location', 'best', 'Interpreter', 'latex', ...
%     'NumColumns', 4)  % Changed to 4 columns
% grid on
% grid minor
% box on
% ylim([-6 8])
% saveas(gcf, 'inertial_drag_averaged.svg')

% % Plot lift data
% figure('Units', 'centimeters', 'Position', [10, 10, 12, 10])
% plot(t_cycle, avg_lift_cycles(:,display_idx), 'LineWidth', 1.5)
% xlabel('Normalized Time (Period)', 'Interpreter', 'latex')
% ylabel('Lift Axis Inertial Force (N)', 'Interpreter', 'latex')
% legend(legend_names(display_idx), 'Location', 'best', 'Interpreter', 'latex', ...
%     'NumColumns', 4)  % Changed to 4 columns
% grid on
% grid minor
% box on
% ylim([-6 8])
% saveas(gcf, 'inertial_lift_averaged.svg')


