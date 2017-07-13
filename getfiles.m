%% little wrapper to get all csv files (containing matrices) into a single matlab matrix
function [matrix files] = getfiles(direc) 

files = dir(fullfile(direc,'*.csv'));

for i = 1:length(files)
    
    temp = csvread(files(i).name);
    matrix(:,:,i) =  temp;
    
end