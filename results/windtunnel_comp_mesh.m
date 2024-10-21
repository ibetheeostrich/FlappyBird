clear;
close all;
clc;


myDir = './ldvm_results';%gets directory
myFiles = dir(fullfile(myDir,'*.csv')); %gets all wav files in struct

for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  [wavData, Fs] = readtable(fullFileName);
  % all of your actions for filtering and plotting go here
end