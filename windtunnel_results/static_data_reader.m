clear; clc; close all;
GraphGood();

folder_path = '.\STATIC TESTS (RAW)';
new_folder_path = '.\STATIC TESTS (CLEANED)';

% Get a list of all .txt files in the selected folder
cd(folder_path)

file_list = dir('*.txt');

cd ..


% Open or create a CSV file for writing
csv_filename = 'static_data_summary.csv';
if ~isfile(csv_filename)
    % Write header if the file does not exist
    writecell({'Lift Mean', 'Drag Mean', 'Phase Angle', 'Velocity', 'AOA'}, csv_filename);
end

% Initialize progress bar with a custom window name
h = waitbar(0, 'Processing files...', 'Name', 'your ');

% Loop through each file in the folder
for i = 1:length(file_list)
    % Update progress bar
    waitbar(i / length(file_list), h, sprintf('Processing file %d of %d...', i, length(file_list)));
    
    % Construct full file path
    file_path = fullfile(folder_path, file_list(i).name);

        % Extract conditions from filename
    [~, filename, ~] = fileparts(file_list(i).name);
    
    % Use regular expressions to extract information
    pattern = '(-?\d+(?:\.\d+)?)deg_?(\d+(?:\.\d+)?)?ms_?(\d+(?:\.\d+)?)?deg?';
    tokens = regexp(filename, pattern, 'tokens');
    
    if ~isempty(tokens)
        aoa = str2double(tokens{1}{1});
        velocity = str2double(tokens{1}{2});
        phase_angle = str2double(tokens{1}{3});
    else
        aoa = NaN;
        velocity = NaN;
        phase_angle = NaN;
    end
    
    % fprintf('Filename: %s\n', filename);
    fprintf('Angle of Attack: %.1f deg\n', aoa);
    fprintf('Wind Tunnel Velocity: ');
    if isnan(velocity)
        fprintf('NaN m/s\n');
    else
        fprintf('%.1f m/s\n', velocity);
    end
    fprintf('Phase Angle: ');
    if isnan(phase_angle)
        fprintf('NaN deg\n');
    else
        fprintf('%.1f deg\n', phase_angle);
    end
    
    % Read the current file
    data = readtable(file_path);

    lift = -data.Var2;
    drag = data.Var1;
    lift_mean = mean(lift(1:3000));
    drag_mean = mean(drag(1:3000));

    aero_dat = [drag_mean; lift_mean];
    R = [cosd(aoa), -sind(aoa); sind(aoa), cosd(aoa)];
    aero_rot = R \ aero_dat;
    drag_rot = aero_rot(1,:);
    lift_rot = aero_rot(2,:);

    % Prepare data for CSV
    data_row = [lift_rot, drag_rot, phase_angle, velocity, aoa];
    
    % Append data to CSV
    writematrix(data_row, csv_filename, 'WriteMode', 'append');

    % Print results with conditions

    fprintf('Lift: %.4f N, Drag: %.4f N\n', lift_rot, drag_rot);
    fprintf('------------------------\n');
end

% Close progress bar
close(h);

% Display completion message
disp('Complete!');






