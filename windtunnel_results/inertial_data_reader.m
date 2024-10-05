clear; clc; close all;
GraphGood();

folder_path = '.\INERTIAL MEASUREMENTS';
new_folder_path = '.\INERTIAL (CLEAN)';

% Get a list of all .txt files in the selected folder
cd(folder_path)

file_list = dir('*.txt');

cd ..

% Get filter object 
Hdd = INERTIAL_FILTER();

% Loop through each file in the folder
for i = 1:length(file_list)
    % Construct full file path
    file_path = fullfile(folder_path, file_list(i).name);
    
    % Read the current file
    data = readtable(file_path);

    lift_inertia = data.Var2;
    drag_inertia = data.Var1;
    time = [0:0.001:(length(lift_inertia)-1)*0.001]';

    new_lift = filtfilt(Hdd.Numerator, 1, lift_inertia);
    new_drag = filtfilt(Hdd.Numerator, 1, drag_inertia);

    % output to csv
    % Combine new_lift and new_drag into a single matrix
    combined_data = [new_drag, new_lift];
    
    % Write the combined data to a single CSV file
    writematrix(combined_data, fullfile(new_folder_path, [file_list(i).name(1:end-4) '.csv']));

    plot(time, new_lift);
    hold on;
    plot(time, new_drag);
    hold off;

end



