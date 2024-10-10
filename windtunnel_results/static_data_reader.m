clear; clc; close all;
GraphGood();

folder_path = '.\STATIC TESTS (RAW)';
new_folder_path = '.\STATIC TESTS (CLEANED)';

% Get a list of all .txt files in the selected folder
cd(folder_path)

file_list = dir('*.txt');

cd ..


% Loop through each file in the folder
for i = 1:1%length(file_list)
    % Construct full file path
    file_path = fullfile(folder_path, file_list(i).name);
    
    % Read the current file
    data = readtable(file_path);

    lift = data.Var2;
    drag = data.Var1;

    lift_mean = mean(lift);
    drag_mean = mean(drag);

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

    % create a structure, using the same name as the file
    data_struct = struct('aoa', aoa, 'velocity', velocity, 'phase_angle', phase_angle, 'lift', lift_mean, 'drag', drag_mean);
    
    % Print results with conditions
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
    fprintf('Lift: %.4f N, Drag: %.4f N\n', lift_mean, drag_mean);
    fprintf('------------------------\n');

    % store all 5 data points in a new structure in matlab, with the same 
    % name as the file

    % plot the data
    % Initialize figures outside the loop if they don't exist yet
    if ~exist('lift_fig', 'var')
        lift_fig = figure('Units', 'centimeters', 'Position', [10, 10, 20, 7.5]);
        title('Lift vs Phase Angle for Various AOA and Velocities', 'Interpreter', 'latex');
        xlabel('Phase Angle (deg)', 'Interpreter', 'latex');
        ylabel('Lift (N)', 'Interpreter', 'latex');
        hold on;
        box on; grid on; grid minor;
    end
    
    if ~exist('drag_fig', 'var')
        drag_fig = figure('Units', 'centimeters', 'Position', [40, 10, 20, 7.5]);
        title('Drag vs Phase Angle for Various AOA and Velocities', 'Interpreter', 'latex');
        xlabel('Phase Angle (deg)', 'Interpreter', 'latex');
        ylabel('Drag (N)', 'Interpreter', 'latex');
        hold on;
        box on; grid on; grid minor;
    end
    
    % Plot lift data
    figure(lift_fig);
    field_name = sprintf('AOA%.1f_V%.1f', abs(aoa), velocity);
    field_name = strrep(field_name, '.', '_');  % Replace dots with underscores
    if ~isfield(lift_fig, 'UserData') || ~isfield(lift_fig.UserData, field_name)
        color = rand(1,3);
        lift_fig.UserData.(field_name) = color;
    else
        color = lift_fig.UserData.(field_name);
    end
    plot(phase_angle, lift_mean, 'o', 'LineWidth', 1.7, 'Color', color);
    
    % Plot drag data
    figure(drag_fig);
    if ~isfield(drag_fig, 'UserData') || ~isfield(drag_fig.UserData, sprintf('AOA%.1f_V%.1f', aoa, velocity))
        color = rand(1,3);
        drag_fig.UserData.(sprintf('AOA%.1f_V%.1f', aoa, velocity)) = color;
    else
        color = drag_fig.UserData.(sprintf('AOA%.1f_V%.1f', aoa, velocity));
    end
    plot(phase_angle, drag_mean, 'o', 'LineWidth', 1.7, 'Color', color);

end

% phase angle always X axis
% plot parametric data e.g. lift vs phase angle for each velocity and overlay


% plot paramatric data e.g. drag vs phase angle for each velocity and overlay
% plot parametric data e.g. lift vs phase angle for each AOA and overlay    
% plot parametric data e.g. drag vs phase angle for each AOA and overlay

