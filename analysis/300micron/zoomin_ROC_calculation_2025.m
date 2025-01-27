%% ------------------------------------------------------------------------
%% zoomin_ROC_calculation_2025.m
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

%% Experiment information
addpath('functions');
basePath = '../../../experiments/300micron/';
expName = '2025-01-20-bubblyice/';
subfolder = 'zoomin_boundaries/';
pathToBoundaries = [basePath,expName,subfolder];

fc = func_curve();

%% Identify boundary file names in pathToBoundaries
files = dir(pathToBoundaries);
fileNames = {files(~[files.isdir] & ~strcmp({files.name}, '.DS_Store')).name}; % cells {'0.mat'}, etc.
ts_arr = zeros(size(fileNames));
for j = 1:length(fileNames)
    ts_arr(j) = str2double(fileNames{j}(1:end-4));
end
clear files; clear fileNames;

%% Iterate over all times and compute radius of curvature
peakbox = 400;%601.6/2; % in pixels, for one side
% Seems to be stable for 330-1000 pixels.
for j = 1:length(ts_arr)
    ts = ts_arr(j); % specific time stamp
    load([pathToBoundaries,num2str(ts),'.mat']); % load data from that time stamp
    
    % Find points around the peak according to peakbox size.
    [peak_x,peak_z,xpeakloc,zpeakloc] = isolatePeakPoints(x_pxl,z_pxl,peakbox);
    % % Symmetrize if off-balance:
    % epsilon = 1; % about 4 microns
    % [peak_x,peak_z,peakbox] = symmetrize(peak_x,peak_z,epsilon);
    % Convert to centimeter
    peak_x = peak_x * convratio; peak_z = peak_z * convratio;
    
    % Polynomial fit
    rr1 = linspace(min(peak_x),max(peak_x),10000);
    [p,opt_deg] = degree_check(peak_x,peak_z,rr1,16);
    zz = polyval(p,rr1);

    % compute tip curvature (from S.W.)
    rk = computeTipCurvature(p,rr1,opt_deg);
    
    disp(['box size: ', num2str(peakbox), ' pixels'...
        newline 'radius of curvature: ',num2str(rk*1e4),' microns'])    

    %% Plot peakbox with rectangle
    hold on;
    xleft = (xpeakloc-peakbox)*convratio;
    xright = (xpeakloc+peakbox)*convratio;
    xline(xleft); xline(xright);
    plot(peak_x,peak_z,'o');
    axis equal
    plot(rr1,zz,'Linewidth',1.2)
    xlim([1.2*xleft,1.2*xright])
    title(['Peakbox with interpolant degree ',num2str(opt_deg),'. Peakbox size = ',num2str(peakbox*convratio),' cm. Local R = ',num2str(round(rk*1e4)),' microns'])
    xlabel('cm'); ylabel('cm')
end