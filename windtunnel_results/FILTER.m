function Hd = FILTER
%FILTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.14 and Signal Processing Toolbox 9.2.
% Generated on: 04-Oct-2024 14:39:19

% FIR Window Lowpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 1000;  % Sampling Frequency

N    = 100;       % Order (further reduced)
Fc   = 12;       % Cutoff Frequency (significantly increased)
flag = 'scale';  % Sampling Flag
Beta = 0.2;      % Window Parameter (further reduced)

% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'low', win, flag);
Hd = dfilt.dffir(b);

% [EOF]
