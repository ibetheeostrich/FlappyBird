clear;
close all;
clc;
GraphGood();

extrac1p = false;

% File Directories
folder_path = '.\DYNAMIC (SELECTED)';
clean_folder_path = '.\DYNAMIC (SELECTED)';
clean_folder_path_1p = '.\DYNAMIC ALL (CLEAN1P)';
victor_dynamic_dodgy = '.\DYNAMIC DODGY';
figs_folder_path = '.\DYNAMIC FIGS (DODGY)';
% Get list of all victor's files
cd(victor_dynamic_dodgy)
victor_files = dir('*.csv');
cd ..

% % Loop through all victor's files and plot
% for i = 1:length(victor_files)
%     file_path = fullfile(victor_dynamic_dodgy, victor_files(i).name);
%     data = readtable(file_path);
%     time = data.t;
%     lift = data.lift;
%     % drag = data.drag;
% 
%     figure('Units', 'centimeters', 'Position', [10, 10, 12, 10]);  
%     plot(time, lift, 'LineWidth', 1.7);
%     hold on;
%     title(['Victor Dynamic Data: ', char(victor_files(i).name)]);
%     xlabel('Time (s)');
%     ylabel('Force (N)');
%     legend('Lift', 'Location', 'best');
%     grid minor;
% end

% Get a list of all .txt files in the selected folder
cd(folder_path)
file_list = dir('*.csv');
cd ..

% Get filter object
Hd = FILTER();

% Define x-limits for each index
xlims = [
    1720-100, 2172-100;
    1734-50, 2134-50;
    1805-50, 2305-50;
    1886-100, 2286-100;
    2067-50, 2567-50;
    1245, 2245
];

% Loop through each file in the folder
for i = 1:length(file_list)
    % Construct full file path
    file_path = fullfile(folder_path, file_list(i).name);

    % Extract conditions from filename
    [~, filename, ~] = fileparts(file_list(i).name);

    % Use regular expressions to extract information
    pattern = '(-?\d+(?:\.\d+)?)deg_?(\d+(?:\.\d+)?)?ms_?(\d+(?:\.\d+)?)?Hz?';
    tokens = regexp(filename, pattern, 'tokens');

    if ~isempty(tokens)
        aoa = str2double(tokens{1}{1});
        velocity = str2double(tokens{1}{2});
        frequency = str2double(tokens{1}{3});
    else
        aoa = NaN;
        velocity = NaN;
        frequency = NaN;
    end

    % Read the current file
    data = readtable(file_path);

    drag = data.Drag_N_';
    lift = data.Lift_N_';

    aero_dat = [drag; lift];

    % Rotation matrix of lift and drag
    R = [cosd(aoa), -sind(aoa); sind(aoa), cosd(aoa)];
    aero_rot = R \ aero_dat;
    drag_rot = aero_rot(1,:)';
    lift_rot = aero_rot(2,:)';

    time = [0:0.001:(length(lift)-1)*0.001]';

    new_lift = filtfilt(Hd.Numerator, 1, lift');
    new_drag = filtfilt(Hd.Numerator, 1, drag');

    % Create a new table with filtered lift and drag data
    filtered_data = table(new_drag, new_lift, time, 'VariableNames', {'Drag (N)', 'Lift (N)', 'Time (s)'});

    % Construct the output filename
    output_filename = fullfile(clean_folder_path, [filename '.csv']);

    % Write the filtered data to a new CSV file
    writetable(filtered_data, output_filename);

    fprintf('Saved filtered data to: %s\n', output_filename);

    % Plot the filtered data
    figure('Units', 'centimeters', 'Position', [10, 10, 12, 10]);  
    plot(time(xlims(i, 1):xlims(i, 2))-time(xlims(i,1)), new_lift(xlims(i, 1):xlims(i, 2)), 'LineWidth', 1.7);
    hold on;
    title(['Dynamic: ', ' | Velocity: ', num2str(velocity), ' m/s | AOA: ', num2str(aoa), ' deg| Frequency: ', num2str(frequency), ' Hz']);
    xlabel('Time (s)');
    ylabel('Force (N)');
    % xlim(xlims(i, :)); % Set specific x-limits based on index
    grid minor;
    % legend('Lift', 'Location', 'best');

    % Check if the same file exists in victor_files and plot if it does
    victor_file_path = fullfile(victor_dynamic_dodgy, [filename '.csv']);
    if isfile(victor_file_path)
        victor_data = readtable(victor_file_path);
        victor_time = victor_data.t;
        victor_lift = victor_data.lift;

        % Plot Victor's data
        plot(victor_time(5:end), victor_lift(5:end), 'LineWidth', 1.7);
        xlim([0, victor_time(end)]);
        legend('Wind Tunnel', 'Analytical', 'Location', 'best');
    end

    % Display progress
    fprintf('Processed file %d of %d: %s\n', i, length(file_list), file_list(i).name);

    saveas(gcf, fullfile(figs_folder_path, [filename '.png']));
end

% Add any final operations or data compilation here if needed
