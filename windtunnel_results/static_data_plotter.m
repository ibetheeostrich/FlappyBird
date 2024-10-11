close all; clear; clc;
% Read the CSV data

GraphGood();


data = readmatrix('static_data_summary.csv', 'NumHeaderLines', 1);

% Sort data by Velocity, AOA, and then Phase Angle
sorted_data = sortrows(data, [4, 5, 3]); % Sort by Velocity (4th column), AOA (5th column), Phase Angle (3rd column)
 
% Extract sorted columns
lift_mean = sorted_data(:, 1);
drag_mean = sorted_data(:, 2);
phase_angle = sorted_data(:, 3);
velocity = sorted_data(:, 4);
aoa = sorted_data(:, 5);

% Extract unique AOAs
unique_aoas = unique(aoa);

%% LIFT + DRAG VS PHASE ANGLE VS VELOCITY
% close all; clc;
for a = 1:length(unique_aoas)
    current_aoa = unique_aoas(a);
    
    % Filter data for the current AOA
    current_aoa_data = sorted_data(aoa == current_aoa, :);
    
    % Extract relevant columns
    lift_mean_current_aoa = current_aoa_data(:, 1);
    drag_mean_current_aoa = current_aoa_data(:, 2);
    phase_angle_current_aoa = current_aoa_data(:, 3);
    velocity_current_aoa = current_aoa_data(:, 4);
    
    % Plot lift vs phase angle for each velocity
    figure('Units', 'centimeters', 'Position', [10, 10, 12, 10]);
    hold on;
    unique_velocities = unique(velocity_current_aoa);
    
    for v = 1:length(unique_velocities)
        % Filter data for the current velocity
        current_velocity_data = current_aoa_data(velocity_current_aoa == unique_velocities(v), :);
        
        % Extract phase angle and lift for the current velocity
        phase_angle_current = current_velocity_data(:, 3);
        lift_current = current_velocity_data(:, 1);
        
        % Plot lift
        plot(phase_angle_current, lift_current, 'o-', 'LineWidth', 1.5, 'DisplayName', sprintf('%.1f m/s', unique_velocities(v)));
    end
    
    title(sprintf('Lift vs Phase Angle for AOA = %.1f deg', current_aoa), 'Interpreter', 'latex');
    xlabel('Phase Angle (deg)', 'Interpreter', 'latex');
    ylabel('Lift Mean (N)', 'Interpreter', 'latex');
    legend('show', 'Location', 'best');
    grid on; box on; grid minor;
    xlim([0, 360]);
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:40:360;
    ax.XTick = 0:40:360; % Set major axis to have 40 deg increments
    hold off;
    
    % Plot drag vs phase angle for each velocity
    figure('Units', 'centimeters', 'Position', [10, 10, 12, 10]);
    hold on;
    
    for v = 1:length(unique_velocities)
        % Filter data for the current velocity
        current_velocity_data = current_aoa_data(velocity_current_aoa == unique_velocities(v), :);
        
        % Extract phase angle and drag for the current velocity
        phase_angle_current = current_velocity_data(:, 3);
        drag_current = current_velocity_data(:, 2);
        
        % Plot drag
        plot(phase_angle_current, drag_current, 'x-', 'LineWidth', 1.5, 'DisplayName', sprintf('%.1f m/s', unique_velocities(v)));
    end
    
    title(sprintf('Drag vs Phase Angle for AOA = %.1f deg', current_aoa), 'Interpreter', 'latex');
    xlabel('Phase Angle (deg)', 'Interpreter', 'latex');
    ylabel('Drag Mean (N)', 'Interpreter', 'latex');
    legend('show', 'Location', 'best');
    grid on; box on; grid minor;
    xlim([0, 360]);
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:40:360;
    ax.XTick = 0:40:360; % Set major axis to have 40 deg increments
    hold off;
end

%% LIFT + DRAG VS PHASE ANGLE VS VELOCITY (Carpet Plot)
% close all; clc;
for a = 1:length(unique_aoas)
    current_aoa = unique_aoas(a);
    
    % Filter data for the current AOA
    current_aoa_data = sorted_data(aoa == current_aoa, :);
    
    % Extract relevant columns
    phase_angle_current_aoa = current_aoa_data(:, 3);
    velocity_current_aoa = current_aoa_data(:, 4);
    lift_mean_current_aoa = current_aoa_data(:, 1);
    drag_mean_current_aoa = current_aoa_data(:, 2);
    
    % Example preprocessing to remove duplicates
    [unique_points, ~, idx] = unique([phase_angle_current_aoa, velocity_current_aoa], 'rows');
    lift_mean_current_aoa = accumarray(idx, lift_mean_current_aoa, [], @mean);
    drag_mean_current_aoa = accumarray(idx, drag_mean_current_aoa, [], @mean);
    
    % Now use unique_points for griddata
    phase_angle_current_aoa = unique_points(:, 1);
    velocity_current_aoa = unique_points(:, 2);
    
    % Create meshgrid for phase angle and velocity
    [phase_grid, velocity_grid] = meshgrid(unique(phase_angle_current_aoa), unique(velocity_current_aoa));
    
    % Interpolate lift and drag data
    lift_grid = griddata(phase_angle_current_aoa, velocity_current_aoa, lift_mean_current_aoa, phase_grid, velocity_grid);
    drag_grid = griddata(phase_angle_current_aoa, velocity_current_aoa, drag_mean_current_aoa, phase_grid, velocity_grid);
    
    % Plot lift carpet plot
    figure('Units', 'centimeters', 'Position', [10, 10, 12, 10]);
    surf(phase_grid, velocity_grid, lift_grid);
    hold on;
    % % Add dashed lines from the XY plane to the surface
    % for i = 1:numel(phase_grid)
    %     plot3([phase_grid(i), phase_grid(i)], [velocity_grid(i), velocity_grid(i)], [0, lift_grid(i)], 'k--');
    % end
    hold off;
    title(sprintf('Lift Carpet Plot for AOA = %.1f deg', current_aoa), 'Interpreter', 'latex');
    xlabel('Phase Angle (deg)', 'Interpreter', 'latex');
    ylabel('Velocity (m/s)', 'Interpreter', 'latex');
    zlabel('Lift Mean (N)', 'Interpreter', 'latex');
    grid on; box on;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:40:360;
    ax.XTick = 0:40:360; % Set major axis to have 40 deg increments

    % Plot drag carpet plot
    figure('Units', 'centimeters', 'Position', [10, 10, 12, 10]);
    surf(phase_grid, velocity_grid, drag_grid);
    hold on;
    % % Add dashed lines from the XY plane to the surface
    % for i = 1:numel(phase_grid)
    %     plot3([phase_grid(i), phase_grid(i)], [velocity_grid(i), velocity_grid(i)], [0, drag_grid(i)], 'k--');
    % end
    hold off;
    title(sprintf('Drag Carpet Plot for AOA = %.1f deg', current_aoa), 'Interpreter', 'latex');
    xlabel('Phase Angle (deg)', 'Interpreter', 'latex');
    ylabel('Velocity (m/s)', 'Interpreter', 'latex');
    zlabel('Drag Mean (N)', 'Interpreter', 'latex');
    grid on; box on;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:40:360;
    ax.XTick = 0:40:360; % Set major axis to have 40 deg increments
end

% ... existing code ...

%% LIFT + DRAG VS PHASE ANGLE VS AOA

% close all; clc;

for v = 1:length(unique_velocities)
    current_velocity = unique_velocities(v);
    
    % Filter data for the current velocity
    current_velocity_data = sorted_data(velocity == current_velocity, :);
    
    % Extract relevant columns
    lift_mean_current_velocity = current_velocity_data(:, 1);
    drag_mean_current_velocity = current_velocity_data(:, 2);
    phase_angle_current_velocity = current_velocity_data(:, 3);
    aoa_current_velocity = current_velocity_data(:, 5);
    
    % Plot lift vs phase angle for each AOA
    figure('Units', 'centimeters', 'Position', [10, 10, 12, 10]);
    hold on;
    unique_aoas = unique(aoa_current_velocity);
    
    for a = 1:length(unique_aoas)
        % Filter data for the current AOA
        current_aoa_data = current_velocity_data(aoa_current_velocity == unique_aoas(a), :);
        
        % Extract phase angle and lift for the current AOA
        phase_angle_current = current_aoa_data(:, 3);
        lift_current = current_aoa_data(:, 1);
        
        % Plot lift
        plot(phase_angle_current, lift_current, 'o-', 'LineWidth', 1.5, 'DisplayName', sprintf('AOA = %.1f deg', unique_aoas(a)));
    end
    
    title(sprintf('Lift vs Phase Angle for Velocity = %.1f m/s', current_velocity), 'Interpreter', 'latex');
    xlabel('Phase Angle (deg)', 'Interpreter', 'latex');
    ylabel('Lift Mean (N)', 'Interpreter', 'latex');
    legend('show', 'Location', 'best');
    grid on; box on; grid minor;
    xlim([0, 360]);
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:40:360;
    ax.XTick = 0:40:360; % Set major axis to have 40 deg increments
    hold off;
    
    % Plot drag vs phase angle for each AOA
    figure('Units', 'centimeters', 'Position', [10, 10, 12, 10]);
    hold on;
    
    for a = 1:length(unique_aoas)
        % Filter data for the current AOA
        current_aoa_data = current_velocity_data(aoa_current_velocity == unique_aoas(a), :);
        
        % Extract phase angle and drag for the current AOA
        phase_angle_current = current_aoa_data(:, 3);
        drag_current = current_aoa_data(:, 2);
        
        % Plot drag
        plot(phase_angle_current, drag_current, 'x-', 'LineWidth', 1.5, 'DisplayName', sprintf('AOA = %.1f deg', unique_aoas(a)));
    end
    
    title(sprintf('Drag vs Phase Angle for Velocity = %.1f m/s', current_velocity), 'Interpreter', 'latex');
    xlabel('Phase Angle (deg)', 'Interpreter', 'latex');
    ylabel('Drag Mean (N)', 'Interpreter', 'latex');
    legend('show', 'Location', 'best');
    grid on; box on; grid minor;
    xlim([0, 360]);
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:40:360;
    ax.XTick = 0:40:360; % Set major axis to have 40 deg increments
    hold off;
end

% ... existing code ...

%% LIFT + DRAG VS PHASE ANGLE VS AOA (Carpet Plot)
% close all; clc;
for v = 1:length(unique_velocities)
    current_velocity = unique_velocities(v);
    
    % Filter data for the current velocity
    current_velocity_data = sorted_data(velocity == current_velocity, :);
    
    % Extract relevant columns
    phase_angle_current_velocity = current_velocity_data(:, 3);
    aoa_current_velocity = current_velocity_data(:, 5);
    lift_mean_current_velocity = current_velocity_data(:, 1);
    drag_mean_current_velocity = current_velocity_data(:, 2);
    
    % Example preprocessing to remove duplicates
    [unique_points, ~, idx] = unique([phase_angle_current_velocity, aoa_current_velocity], 'rows');
    lift_mean_current_velocity = accumarray(idx, lift_mean_current_velocity, [], @mean);
    drag_mean_current_velocity = accumarray(idx, drag_mean_current_velocity, [], @mean);
    
    % Now use unique_points for griddata
    phase_angle_current_velocity = unique_points(:, 1);
    aoa_current_velocity = unique_points(:, 2);

    % Create meshgrid for phase angle and AOA
    [phase_grid, aoa_grid] = meshgrid(unique(phase_angle_current_velocity), unique(aoa_current_velocity));
    
    % Interpolate lift and drag data
    lift_grid = griddata(phase_angle_current_velocity, aoa_current_velocity, lift_mean_current_velocity, phase_grid, aoa_grid);
    drag_grid = griddata(phase_angle_current_velocity, aoa_current_velocity, drag_mean_current_velocity, phase_grid, aoa_grid);
    
    % Plot lift carpet plot
    figure;
    surf(phase_grid, aoa_grid, lift_grid);
    hold on;
    % % Add dashed lines from the XY plane to the surface
    % for i = 1:numel(phase_grid)
    %     plot3([phase_grid(i), phase_grid(i)], [aoa_grid(i), aoa_grid(i)], [0, lift_grid(i)], 'k--');
    % end
    hold off;
    title(sprintf('Lift Carpet Plot for Velocity = %.1f m/s', current_velocity), 'Interpreter', 'latex');
    xlabel('Phase Angle (deg)', 'Interpreter', 'latex');
    ylabel('AOA (deg)', 'Interpreter', 'latex');
    zlabel('Lift Mean (N)', 'Interpreter', 'latex');
    grid on; box on;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:40:360;
    ax.XTick = 0:40:360; % Set major axis to have 40 deg increments
    ax.YTick = -5:2.5:5;

    % Plot drag carpet plot
    figure;
    surf(phase_grid, aoa_grid, drag_grid);
    hold on;
    % % Add dashed lines from the XY plane to the surface
    % for i = 1:numel(phase_grid)
    %     plot3([phase_grid(i), phase_grid(i)], [aoa_grid(i), aoa_grid(i)], [0, drag_grid(i)], 'k--');
    % end
    hold off;
    title(sprintf('Drag Carpet Plot for Velocity = %.1f m/s', current_velocity), 'Interpreter', 'latex');
    xlabel('Phase Angle (deg)', 'Interpreter', 'latex');
    ylabel('AOA (deg)', 'Interpreter', 'latex');
    zlabel('Drag Mean (N)', 'Interpreter', 'latex');
    grid on; box on;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:40:360;
    ax.XTick = 0:40:360; % Set major axis to have 40 deg increments
    ax.YTick = -5:2.5:5;
end

% ... existing code ...