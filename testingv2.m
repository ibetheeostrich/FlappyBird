clear;
close all;
clc;

folder_path = '.\windtunnel_results\DYNAMIC ALL (CLEAN1P)';
filename = '0.0deg_12.0ms_2.5Hz.csv'; 

file_path = fullfile(folder_path, filename);

data_wt = readtable(file_path);

folder_path = '.\';
filename = '0.0010000_0.0deg_12.0ms_2.5Hz_0.24LESP.csv'; 

file_path = fullfile(folder_path, filename);

data_ldvm = readtable(file_path);

% Get filter object
Hd = FILTER();

new_lift_ldvm = filtfilt(Hd.Numerator, 1, data_ldvm.Var2');%_rot);
new_drag_ldvm = filtfilt(Hd.Numerator, 1, data_ldvm.Var1');%_rot);

figure(1)
plot(data_ldvm.Var3,new_lift_ldvm)
hold on
plot(data_wt.Time_s_, data_wt.Lift_N_)