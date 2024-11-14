clear;
close all;
clc;
GraphGood();



extrac1p = false;



% File Directories
folder_path = '.\DYNAMIC ALL (RAW)';
clean_folder_path = '.\DYNAMIC ALL (CLEAN)';
clean_folder_path_1p = '.\DYNAMIC ALL (CLEAN1P)';

figs_folder_path = '.\DYNAMIC FIGS';


% Get a list of all .txt files in the selected folder
cd(folder_path)

file_list = dir('*.txt');

cd ..

% Get filter object
Hd = FILTER();

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

    drag = data.Var1';
    lift = data.Var2';

    aero_dat = [drag; lift];

    % extract angle of attack from filename (first 3 characters)
    % aoa = str2double(file_list(i).name(1:3));

    % rotation matrix of lift and drag
    R = [cosd(aoa), -sind(aoa); sind(aoa), cosd(aoa)];
    aero_rot = R \ aero_dat;
    drag_rot = aero_rot(1,:)';
    lift_rot = aero_rot(2,:)';

    time = [0:0.001:(length(lift)-1)*0.001]';

    new_lift = filtfilt(Hd.Numerator, 1, lift');%_rot);
    new_drag = filtfilt(Hd.Numerator, 1, drag');%_rot);

    % Extract the original filename without extension
    [~, filename, ~] = fileparts(file_list(i).name);

    % Create a new table with filtered lift and drag data
    filtered_data = table(new_drag, new_lift, time, 'VariableNames', {'Drag (N)', 'Lift (N)' , 'Time (s)'});

    % Construct the output filename
    output_filename = fullfile(clean_folder_path, [filename '.csv']);

    % Write the filtered data to a new CSV file
    writetable(filtered_data, output_filename);

    fprintf('Saved filtered data to: %s\n', output_filename);

    if extrac1p
        testparams = strsplit(filename,'_');

        % Extract frequency of experiment
        freq = testparams(3);
        freq = char(freq);
        freq = str2double(freq(1:end-2));

        % Find intercepts for one wavelength
        ind = [];
        for j = 1:length(new_lift)-1
            if new_lift(j) > 0 && new_lift(j+1) < 0
                ind(max(size(ind))+1) = j;
            end
        end

        % Check that intercepts exist
        if size(ind) ==0

            % Check the domain is actually approximately one wavelength
        elseif max(size(ind)) > 0 && (abs((time(ind(3)) - time(ind(2))) - 1/freq) < 0.25/freq)

            % Extract lift and drag for one period
            new_lift1p = new_lift(ind(2):ind(3));
            new_drag1p = new_drag(ind(2):ind(3));
            new_time1p = time(ind(2):ind(3));

            % Frequency scaling of time domain
            time_scaled = linspace(0,1/freq,length(new_time1p))';

            output_filename = fullfile(clean_folder_path_1p, [filename '.csv']);
            filtered_data_1p = table(new_drag1p, new_lift1p, time_scaled, 'VariableNames', {'Drag (N)', 'Lift (N)' , 'Time (s)'});
            writetable(filtered_data_1p, output_filename);
            figure(i)
            plot(time_scaled,new_lift1p,'LineWidth',1.7)
            hold on
            % plot(time_scaled,lift(ind(2):ind(3)), 'Color', [0,0,1,0.1])
            plot(time_scaled,new_drag1p,'LineWidth',1.7)
            % plot(time_scaled,drag(ind(2):ind(3)), 'Color', [1,0,0,0.1])
            title(['Dynamic: ', ' | Velocity: ', num2str(velocity), ' m/s | AOA: ', num2str(aoa), ' deg | Frequency: ', num2str(frequency), ' Hz'])
            xlabel('Time (s)')
            ylabel('Force (N)')
            xlim([0,5]);
            grid minor
            legend('Filtered Lift', 'Filtered Drag', 'Location', 'best', 'Orientation', 'horizontal', 'NumColumns', 2)

        end
    else
        figure('Units', 'centimeters', 'Position', [10, 10, 12, 8]);  
        plot(time,new_lift,'LineWidth',1.7)
        hold on
        % plot(time,lift, 'Color', [0,0,1,0.1])
        plot(time,new_drag,'LineWidth',1.7)
        % plot(time,drag, 'Color', [1,0,0,0.1])
%         title(['Vel: ', num2str(velocity), ' m/s | AOA: ', num2str(aoa), ' deg| Freq: ', num2str(frequency), ' Hz'])
        xlabel('Time (s)')
        ylabel('Force (N)')
        xlim([0,2]);
        grid minor
        legend('Lift', 'Drag', 'Location', 'best', 'Orientation', 'horizontal', 'NumColumns', 2)
        
        saveas(gcf, fullfile(figs_folder_path, ['dynamic' filename '.svg']));
    end

    % Display progress
    fprintf('Processed file %d of %d: %s\n', i, length(file_list), file_list(i).name);
end

% Add any final operations or data compilation here if needed



