%% Hyper_Analysis
clear all
close all
clc

%% initial parameters

fold = '240920/';       % insert search path here

backper = 0.1;                 %fraction of spectra to use in background average
lowercut=50;                   %Pixels to cut from the blue end of spectra
uppercut=50;                   %Pixels to cut from the red side of spectra
lower = 0;                    %lower bound for particle identification (controls how dim of an object is a particle)
upper = 0.5;                      %upper bound for particle identification (controls how bright of an object is selected)
nhood = 1;                      %odd number: nhood by nhood pixels binned to make spectra
rsquarelim = 0.5;
binfac = 'one';                 %binning factor 'one' 'two' or 'four'
bf = 2;                         %1,2,4filename2


%% Opening files, analyzing and saving data in mat  structures
clc;
addpath(fold)               
dataloc = fold;
addpath(dataloc)                            % adds directory to path
fid = fopen([dataloc,'/mydata.txt']);       % opens mydata and creates string array
names = textscan(fid, '%s');    
fclose(fid);
nsamp = numel(names{1,1})-2;                % counts number of files to analyze

fcheck=exist([dataloc '/standark.mat'],'file'); % is there a standard already here?
if fcheck==2
    load([dataloc '/standark.mat'])         % if yes, load it
    stanbig = stan;

else
    stan = standardreadindark(strcat(names{1,1}(1),'.tdms'),strcat(names{1,1}(2),'.tdms'),dataloc); % if no, read in the files and make one
    
end
stan = stan - min(min(min(stan))) + 0.1;

for c3 = 1:nsamp        
    anfunc_lorentz_fit(dataloc,names,c3,stan,backper,lowercut,uppercut,lower,upper,nhood,rsquarelim,binfac)  % analyze all of the files listed in mydata using another function for memory purposes
end


