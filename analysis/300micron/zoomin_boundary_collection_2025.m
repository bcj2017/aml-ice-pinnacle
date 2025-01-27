%% ------------------------------------------------------------------------
%% zoomin_boundary_collection.m
%  Modified from Steven's old code by Bobae
%  Purpose of code:
%  Given an experiment, this code will iterate across photos (line 31
%  determines interval between photos) and allow computer+manual boundary
%  collection for each zoomed-in boundary.
%
%  Line 32 sets the time spacing between each photo (interval timer
%  shooting).
%
%  Boundaries are saved as 'ts.mat' ('ts' is time in seconds) so they can 
%  be reused for all calculations later.
%
%  To be run first following experiments to save boundaries!
%% ------------------------------------------------------------------------

close all
clear
warning('off')
addpath('functions');

basepath = '../../../experiments/300micron/';
filepath = '2025-01-20-bubblyice'; % insert experiment name
basepath = [basepath,filepath];
subfolder = '/zoomin/nicest/'; % NOTE: remove 'nicest' for regular experiments
conv_index = 1; 

fc = func_curve();

convratio = calibra(basepath,subfolder,1); % run calibration function, convratio is cm/pixel

path = [basepath,subfolder];
S = dir(fullfile(path,'**','*.JPG'));
names = {S.name};
disp(names);

time_spacing = 10; % 10 second interval timer shooting
ts = 0;
photo_interval = 1
for jjj = 1:length(names)
    disp(['Current time stamp: ',num2str(ts),' seconds']);
    if mod(jjj,photo_interval) == 0 % used for plotting (every photo_interval photos)
        disp(['Boundary will be computed for current time stamp: ',num2str(ts),' seconds.'])
        % === User Manual === 
        % 0: set as 0; 2: run hist; 3: normal; 
        % 4: pass; other value: change sparse point threshold
        % === End === 
        disp("*************************")
        disp("=== Start Calculation ===")
        disp("*************************")
        execind = 3; 
        
            spimage = char(names(jjj));
            disp(['Image name: ',spimage]) % prints i.e. DSC_3874.JPG

            pic = imread([basepath,subfolder,spimage]);

            ydim = size(pic,1); % dimension
            xdim = size(pic,2);
            % subject to change for every trail, but consistent within the trail
            leftcut = 1700; % cutoff dimension
            rightcut = 5000;
            deg = 4; 

            % figure(1);
            pic = rgb2gray(pic);

            if execind == 2 || execind > 4
                pic = histeq(pic);
            end
            % imshow(pic)

            [BW,~,Gh,Gv] = edge(pic,'sobel');
            G = Gh.^2 + Gv.^2;
            % G(Gh < 0) = 0;

            % figure(2)
            ll = contour(G); % find contour line
            % ll = rmoutliers(ll);
            sub = ll(1,:) <= 0 | ll(2,:) <= 0 | ll(1,:) >= xdim | ll(2,:) >= ydim ...
                | ll(1,:) <= leftcut | ll(1,:) >= rightcut; 
            xx = ll(1,~sub); % x val
            yy = ll(2,~sub); % y val

            [xx,idx] = sort(xx); yy = yy(idx);

            % remove spa
            % rse points
            if execind == 3 || execind == 2
                [ex,ey] = rm_outliers(xx,yy,xdim,ydim,50,50);
                sparthold = 500;
            elseif execind > 4
                sparthold = execind;
                [ex,ey] = rm_outliers(xx,yy,xdim,ydim,50,sparthold);
            end

            % extract boundary strategies
            ext_stra = 2;
            if ext_stra == 1
                [ex,ey] = ext_whole(xx,yy,xdim,ydim,sparthold);
            elseif ext_stra == 2
                [ex,ey] = ext_xmost(ex,ey);
            end

            fr = func_aux();
            [newbdx,newbdy] = func_manual(pic,ex,ey,fr);
    
    %% Rescale boundary 
    % Sort boundary to be increasing in x:
    [newbdx,sort_idx] = sort(newbdx); newbdy = newbdy(sort_idx); 
    newbdy = -newbdy;

    newbdy = newbdy - max(newbdy); % Vertically move peak to y = 0
    % Horizontally move peak to about x = 0
    [~,peak_idx] = max(newbdy); % identify peak index from y values
    newbdx = newbdx - newbdx(peak_idx); % horizontal shift

    newbdx_cm = newbdx * convratio; % convert to cm
    newbdy_cm = newbdy * convratio; % convert to cm, flip rightside up
    
    %% Rename variables to be nice
    x_cm = newbdx_cm; z_cm = newbdy_cm;
    x_pxl = newbdx; z_pxl = newbdy;

    %% Save boundary file
    filename = [basepath,'/zoomin_boundaries/',num2str(ts),'.mat'];
    % saveBool = input('Press 1 to save: ');
    % if saveBool == 1
    save(filename,'x_cm','z_cm','x_pxl','z_pxl','convratio');
    disp(['Saved Data to ', filename,' in pixels and centimeters']); 
    
    end
    
    %% Increase time stamp after
    ts = ts + time_spacing; 
end




















