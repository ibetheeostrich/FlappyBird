clear;
close all;
clc;
GraphGood();

folder_path = 'C:\Users\Reagan\OneDrive - The University of Sydney (Students)\2024 SEM 1\FlappyBird\FlappyBird\windtunnel_results\DYNAMIC ALL (RAW)';
new_folder_path = 'C:\Users\Reagan\OneDrive - The University of Sydney (Students)\2024 SEM 1\FlappyBird\FlappyBird\windtunnel_results\DYNAMIC ALL (CLEAN)';

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
    
    % Read the current file
    data = readtable(file_path);

    lift = data.Var2;
    drag = data.Var1;

    time = [0:0.001:(length(lift)-1)*0.001]';

    new_lift = filtfilt(Hd.Numerator, 1, lift);
    new_drag = filtfilt(Hd.Numerator, 1, drag);

    % Extract the original filename without extension
    [~, filename, ~] = fileparts(file_list(i).name);
    
    % Create a new table with filtered lift and drag data
    filtered_data = table(new_drag, new_lift, time, 'VariableNames', {'Drag (N)', 'Lift (N)' , 'Time (s)'});
    
    % Construct the output filename
    output_filename = fullfile(new_folder_path, [filename '.csv']);
    
    % Write the filtered data to a new CSV file
    writetable(filtered_data, output_filename);
    
    fprintf('Saved filtered data to: %s\n', output_filename);

    figure(i)
    plot(time,new_lift,'LineWidth',1.7)
    hold on

    plot(time,lift, 'Color', [0,0,1,0.1])
    plot(time,new_drag,'LineWidth',1.7)
    plot(time,drag, 'Color', [1,0,0,0.1])
    xlim([0,2])
    
    % Process the data (you can add your processing steps here)
    % For example:
    % processed_data = process_data(a, Hd);
    
    % Optionally, you can save the processed data or perform further operations
    % save_processed_data(processed_data, file_list(i).name);
    
    % Display progress
    fprintf('Processed file %d of %d: %s\n', i, length(file_list), file_list(i).name);
end

% Add any final operations or data compilation here if needed
