clear; clc; close all;
GraphGood();

filename = '3hz_14ms.xlsx'; % replace with your actual file name
data = readtable(filename);

inertia = 'inertial_3hz.xlsx';
inertial_data = readtable(inertia);

% To access specific columns like 'X', 'Y', 'Z', 'M_X', 'M_Y', 'M_Z', 'Time', you can do:
X = data.Y;
Time = data.Time;

X_inertial = inertial_data.X;
Time_inertial = inertial_data.Time;

% Apply the filter to both lift data and inertial data
Hd = FILTER(); % Get the filter object
Hdd = INERTIAL_FILTER();
% X_filtered = filtfilt(Hd.Numerator, 1, X);
% X_inertial_filtered = filtfilt(Hdd.Numerator, 1, X_inertial);

% filter both samples
inertial_lift = filtfilt(Hdd.Numerator, 1, X_inertial);
lift = filtfilt(Hd.Numerator, 1, X);

% extract 1000 samples
inertial_lift_filtered = inertial_lift(300+116:1300+116);
lift_filtered = lift(88+250:1088+250);



% plot and compare the two
figure('Units', 'centimeters', 'Position', [10, 10, 20, 7.5]);  
plot(Time(1:1001), lift_filtered, 'LineWidth', 1.7);
hold on;
plot(Time_inertial(1:1001), inertial_lift_filtered, 'LineWidth', 1.7);
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Lift (N)', 'Interpreter', 'latex');
title('Lift vs Time (14 m/s @ 3 Hz)', 'Interpreter', 'latex');
legend('Lift', 'Inertial Lift', 'Interpreter', 'latex');

box on; grid on; grid minor;

% plot filtered acuatl lift
vz = readtable('0.0006250_14.0_3.0test.csv');

% filter victors shit
filter_vz = filtfilt(Hd.Numerator, 1, vz.Var1);


actual_lift = lift_filtered - inertial_lift_filtered;

figure('Units', 'centimeters', 'Position', [10, 10, 12, 10]);  
plot(Time(1:1001), actual_lift, 'LineWidth', 1.7);
hold on
plot(vz.Var3,filter_vz, 'LineWidth', 1.7)
hold off
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Lift (N)', 'Interpreter', 'latex');
title('Lift vs Time (14 m/s @ 3 Hz)', 'Interpreter', 'latex');
xlim([0,1/3]);
box on; grid on; grid minor;
legend('Wind Tunnel Lift', 'Panel Method Lift', 'Interpreter', 'latex', 'Location', 'best');
% xlim([0.287, 0.287+0.33]);
% saveas(gcf, 'filtered_actual_lift_vs_time.svg');

% % actual lift
% actual_lift = lift - inertial_lift;
% 
% % filter actual lift
% actual_lift_filtered = filtfilt(Hd.Numerator, 1, actual_lift);
% 
% % plot filtered acuatl lift
% figure('Units', 'centimeters', 'Position', [10, 10, 20, 7.5]);  
% plot(Time(1:1001), actual_lift_filtered, 'LineWidth', 1.7);
% xlabel('Time (s)', 'Interpreter', 'latex');
% ylabel('Lift (N)', 'Interpreter', 'latex');
% title('Filtered Actual Lift vs Time (14 m/s @ 3 Hz)', 'Interpreter', 'latex');
% box on; grid on; grid minor;
% % xlim([0.287, 0.287+0.33]);
% % saveas(gcf, 'filtered_actual_lift_vs_time.svg');
% 
% % Plot original Lift (X) against Time
% figure('Units', 'centimeters', 'Position', [10, 10, 20, 7.5]);  
% plot(Time, X, 'LineWidth', 1.7);
% xlabel('Time (s)', 'Interpreter', 'latex');
% ylabel('Lift (N)', 'Interpreter', 'latex');
% title('Original Lift vs Time (14 m/s @ 3 Hz)', 'Interpreter', 'latex');
% box on; grid on; grid minor;
% xlim([0, 4]);
% % saveas(gcf, 'original_lift_vs_time.svg');
% 
% % Plot filtered Lift (X) against Time
% figure('Units', 'centimeters', 'Position', [10, 10, 20, 7.5]);  
% plot(Time, X_filtered, 'LineWidth', 1.7);
% xlabel('Time (s)', 'Interpreter', 'latex');
% ylabel('Lift (N)', 'Interpreter', 'latex');
% title('Filtered Lift vs Time (14 m/s @ 3 Hz)', 'Interpreter', 'latex');
% box on; grid on; grid minor;
% xlim([0, 4]);
% % saveas(gcf, 'filtered_lift_vs_time.svg');
% 
% % Plot original Inertial Lift against Time
% figure('Units', 'centimeters', 'Position', [10, 10, 20, 7.5]);  
% plot(Time_inertial, X_inertial, 'LineWidth', 1.7);
% xlabel('Time (s)', 'Interpreter', 'latex');
% ylabel('Inertial Lift (N)', 'Interpreter', 'latex');
% title('Original Inertial Lift vs Time (14 m/s @ 3 Hz)', 'Interpreter', 'latex');
% box on; grid on; grid minor;
% xlim([0, 4]);
% % saveas(gcf, 'original_inertial_lift_vs_time.svg');
% 
% % Plot filtered Inertial Lift against Time
% figure('Units', 'centimeters', 'Position', [10, 10, 20, 7.5]);  
% plot(Time_inertial, X_inertial_filtered, 'LineWidth', 1.7);
% xlabel('Time (s)', 'Interpreter', 'latex');
% ylabel('Inertial Lift (N)', 'Interpreter', 'latex');
% title('Filtered Inertial Lift vs Time (14 m/s @ 3 Hz)', 'Interpreter', 'latex');
% box on; grid on; grid minor;
% xlim([0, 4]);
% % saveas(gcf, 'filtered_inertial_lift_vs_time.svg');
% 
% % Perform FFT on original and filtered Lift data
% Fs = 1000;  % Known sampling frequency of 1000 Hz
% L = length(X);              % Length of signal
% 
% % FFT of original data
% Y = fft(X);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% % FFT of filtered data
% Y_filtered = fft(X_filtered);
% P2_filtered = abs(Y_filtered/L);
% P1_filtered = P2_filtered(1:L/2+1);
% P1_filtered(2:end-1) = 2*P1_filtered(2:end-1);
% 
% f = Fs*(0:(L/2))/L;
% 
% % Plot FFT of original and filtered Lift
% figure('Units', 'centimeters', 'Position', [10, 10, 20, 7.5]);
% plot(f, P1, 'LineWidth', 1.7);
% hold on;
% plot(f, P1_filtered, 'LineWidth', 1.7);
% xlabel('Frequency (Hz)', 'Interpreter', 'latex');
% ylabel('Magnitude', 'Interpreter', 'latex');
% title('FFT of Lift (14 m/s @ 3 Hz)', 'Interpreter', 'latex');
% legend('Original', 'Filtered', 'Interpreter', 'latex');
% box on; grid on; grid minor;
% xlim([0, 4]);
% % saveas(gcf, 'lift_fft.svg');
% 
% % Perform FFT on original and filtered Inertial Lift data
% L_inertial = length(X_inertial);
% 
% % FFT of original inertial data
% Y_inertial = fft(X_inertial);
% P2_inertial = abs(Y_inertial/L_inertial);
% P1_inertial = P2_inertial(1:L_inertial/2+1);
% P1_inertial(2:end-1) = 2*P1_inertial(2:end-1);
% 
% % FFT of filtered inertial data
% Y_inertial_filtered = fft(X_inertial_filtered);
% P2_inertial_filtered = abs(Y_inertial_filtered/L_inertial);
% P1_inertial_filtered = P2_inertial_filtered(1:L_inertial/2+1);
% P1_inertial_filtered(2:end-1) = 2*P1_inertial_filtered(2:end-1);
% 
% f_inertial = Fs*(0:(L_inertial/2))/L_inertial;
% 
% % Plot FFT of original and filtered Inertial Lift
% figure('Units', 'centimeters', 'Position', [10, 10, 20, 7.5]);
% plot(f_inertial, P1_inertial, 'LineWidth', 1.7);
% hold on;
% plot(f_inertial, P1_inertial_filtered, 'LineWidth', 1.7);
% xlabel('Frequency (Hz)', 'Interpreter', 'latex');
% ylabel('Magnitude', 'Interpreter', 'latex');
% title('FFT of Inertial Lift (14 m/s @ 3 Hz)', 'Interpreter', 'latex');
% legend('Original', 'Filtered', 'Interpreter', 'latex');
% box on; grid on; grid minor;
% xlim([0, 4]);
% % saveas(gcf, 'inertial_lift_fft.svg');