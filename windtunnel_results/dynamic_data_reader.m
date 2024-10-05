clear;
close all;
clc;
GraphGood();

folder_path = '.\DYNAMIC ALL (RAW)';
new_folder_path = '.\DYNAMIC ALL (CLEAN)';

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
        time_scaled = linspace(0,1/freq,length(new_time1p));

        figure(i)
        plot(time_scaled,new_lift1p,'LineWidth',1.7)
        hold on
    
        plot(time_scaled,lift(ind(2):ind(3)), 'Color', [0,0,1,0.1])
        plot(time_scaled,new_drag1p,'LineWidth',1.7)
        plot(time_scaled,drag(ind(2):ind(3)), 'Color', [1,0,0,0.1])

    end
    
    % Process the data (you can add your processing steps here)
    % For example:
    % processed_data = process_data(a, Hd);
    
    % Optionally, you can save the processed data or perform further operations
    % save_processed_data(processed_data, file_list(i).name);
    
    % Display progress
    fprintf('Processed file %d of %d: %s\n', i, length(file_list), file_list(i).name);
end

% Add any final operations or data compilation here if needed
