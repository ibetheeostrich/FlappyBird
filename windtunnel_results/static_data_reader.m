clear; clc; close all;
GraphGood();

folder_path = 'C:\Users\Reagan\OneDrive - The University of Sydney (Students)\2024 SEM 1\FlappyBird\FlappyBird\windtunnel_results\STATIC TESTS (RAW)';

% Get a list of all .txt files in the selected folder
cd(folder_path)

file_list = dir('*.txt');

cd ..


% Loop through each file in the folder
for i = 1:length(file_list)
    % Construct full file path
    file_path = fullfile(folder_path, file_list(i).name);
    
    % Read the current file
    data = readtable(file_path);



end