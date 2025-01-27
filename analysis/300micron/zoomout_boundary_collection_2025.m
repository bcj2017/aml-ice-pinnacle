%% --------------------------------------------------------------------------
%% zoomout_boundary_collection_2025.m
%  Modified from Steven's old code by Bobae
%  Purpose of code:
%  Given an experiment, this code will iterate across photos (line 31
%  determines interval between photos) and allow computer+manual boundary
%  collection for each zoomed-out boundary.
%
%  For each boundary, deletes excess ice near base (takes 10 centimeter
%  from t = 0 tip and passes to all other frames).
%
%  Might need to change time spacing 'dt' depending on interval timer
%  shooting...
%
%  Boundaries are saved as 'ts.mat' ('ts' is time in seconds) so they can 
%  be reused for all calculations later.
%
%  To be run first following experiments to save boundaries!
%% --------------------------------------------------------------------------

close all
clear

addpath('functions');

RR = 0.03; % cm
HH = 10; % cm

basepath = '../../../experiments/300micron/';
filepath = '2025-01-20-bubblyice'; % insert experiment name
basepath = [basepath,filepath];
subfolder = '/zoomout/nicest/'; % NOTE: remove 'nicest' for regular experiments

convratio = calibra_zoomout(basepath,subfolder,1); % pixel -> cm ratio

dt = 10; % s % NOTE: change this if we change interval timer shooting interval.
path = [basepath,subfolder];
S = dir(fullfile(path,'**','*.jpeg'));
names = {S.name};
names = names(1:end-1); 

upptime = 1200; % ending time (s) % NOTE: change this later.
intv = 10; % NOTE: difference between this & dt...?
callen = round((upptime/dt - 5)/intv); % 
boundcoll = cell(1,callen); % save tracking result
datacoll = []; % data registration

% functionality index
overlay = 1; 
inorout = 1;
intplot = 1; 

% main loop
lln = length(names);
cnt = 1; 
for i = 1:lln % change back later!
    ts = (i-1) * dt; 
    spimage = char(names(i));
    if i < upptime/dt
        pic = imread([basepath,subfolder,spimage]);
        disp(['==',spimage,'=='])
        % dimension
        xdim = size(pic,1); 
        ydim = size(pic,2);
        % cutoff dimension, relates to size of pinnacle
        leftcut = 50; 
        rightcut = ydim-leftcut;
        % [HYPER]: ignore part of bottom which is likely to be unclear
        basecut = 500; 
        
        pic = rgb2gray(pic);
        pic = imadjust(pic);
        pic = imsharpen(pic);
        pic = imrotate(pic,90);
        
        [BW,~,Gh,Gv] = edge(pic,'sobel');
        G = Gh.^2 + Gv.^2;
        
        % find contour line
        ll = contourc(double(G)); 
        sub = ll(1,:) <= 0 | ll(2,:) <= 0 | ll(1,:) >= xdim | ll(2,:) >= ydim ...
            | ll(1,:) <= leftcut | ll(1,:) >= rightcut | ll(2,:) >= ydim-basecut; % ignore base
        
        xx = ll(1,~sub);
        yy = ll(2,~sub);
        
        [xx,idx] = sort(xx); yy = yy(idx);

        % clear out noise
        [ex,ey] = ext_xmost(xx,yy);
        
        fr = func_aux();
        [newbdx,newbdy] = func_manual(pic,ex,ey,fr);

        %% Rescale boundary 
        % Sort boundary to be increasing in x:
        [newbdx,sort_idx] = sort(newbdx); newbdy = newbdy(sort_idx); 
        newbdy = -newbdy;
        
        % Horizontally move peak to about x = 0
        [~,peak_idx] = max(newbdy); % identify peak index from y values
        newbdx = newbdx - newbdx(peak_idx); % horizontal shift
     
        newbdx_cm = newbdx * convratio; % convert to cm
        newbdy_cm = newbdy * convratio; % convert to cm, flip rightside up
        
        %% Rename variables to be nice
        x_cm = newbdx_cm; z_cm = newbdy_cm;
        x_pxl = newbdx; % z_pxl in line 123

        %% Delete bit near base:
        if i == 1 % first frame: identify y coordinate of "bottom"
            zmin = max(z_cm) - HH;
            save([basepath,'/zoomout_boundaries/data.mat'],'zmin','convratio');
        else
            load([basepath,'/zoomout_boundaries/data.mat']);
        end
        % All frames: delete all elements outside of this range
        indices = find(z_cm > zmin);
        x_cm = x_cm(indices); z_cm = z_cm(indices);
        x_pxl = x_pxl(indices); 
        % Move bottom of pinnacle to be at z = 0:
        z_cm = z_cm - zmin; z_pxl = z_cm / convratio;
        
        %% Save boundary file
        filename = [basepath,'/zoomout_boundaries/',num2str(ts),'.mat'];
        save(filename,'x_cm','z_cm','x_pxl','z_pxl','convratio','indices');
        disp(['Saved Data to ', filename,' in pixels and centimeters']); 
    end
end