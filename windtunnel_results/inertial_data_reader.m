clear; clc; close all;
GraphGood();

folder_path = '.\INERTIAL MEASUREMENTS';
new_folder_path = '.\INERTIAL (CLEAN)';

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

    lift_inertia = data.Var2;
    drag_inertia = data.Var1;
    time = [0:0.001:(length(lift_inertia)-1)*0.001]';

end



