%% ------------------------------------------------------------------------
%% zoomout_closestmatch_calculation_2025.m
%  Modified from Steven's old code by Bobae
%  Purpose of code:
%  Given an experiment, this code loads all of the zoomed-in saved 
%  boundaries for that experiment.
%  Then the radius of curvature is computed for a given peakbox size.
%  The information is then plotted.
%  Useful to run a peakbox stability analysis to choose the peakbox size.
%  It seems that the peakbox stability region may have changed.
%  Not the most urgent though... start looking at zoomed-out POV.
%
%  To be run after zoomin_boundary_collection_2025.m has saved boundaries.
%% ------------------------------------------------------------------------

clear; close all;

RR = 0.03; % cm
HH = 10; % cm
conv_index = 2;

filepath = 'blahfornow';

% color
lightBLUE = [0.356862745098039,0.811764705882353,0.956862745098039];
darkBLUE = [0.0196078431372549,0.0745098039215686,0.670588235294118];
blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));

%% Experiment information
addpath('functions');
basePath = '../../../experiments/300micron/';
expName = '2025-01-20-bubblyice/';
subfolder = 'zoomout_boundaries/';
pathToBoundaries = [basePath,expName,subfolder];

fc = func_curve();

%% Identify boundary file names in pathToBoundaries
files = dir(pathToBoundaries);
fileNames = {files(~[files.isdir] & ~strcmp({files.name}, '.DS_Store')).name}; % cells {'0.mat'}, etc.
% Delete 'data.mat' from fileNames:
idx = strcmp(fileNames, 'data.mat'); fileNames(idx) = [];
ts_arr = zeros(size(fileNames));
for j = 1:length(fileNames)
    ts_arr(j) = str2double(fileNames{j}(1:end-4));
end
clear files; clear fileNames;

for j = 1:length(ts_arr) % iterate over frames
    ts = ts_arr(j);
    % Load data:
    load([pathToBoundaries,num2str(ts),'.mat']);

    [optparamfit, x_interp, z_interp, x_opttheory, z_opttheory] = zo_intp(x_cm,z_cm,filepath,ts_arr(1));
    disp(optparamfit)
    
    filename = [pathToBoundaries,'bestfitdata/',num2str(ts),'_changetotalheight.mat'];
    save(filename,'optparamfit','x_interp','z_interp','x_opttheory','z_opttheory','x_cm','z_cm');
    disp(['Saved: ',num2str(ts),' seconds']);
end