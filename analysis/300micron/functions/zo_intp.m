%--------------------------------------------------------------------------
% Optimization function to the family of curve (Huang & Moore, 2022) in zoom-out view
%
% Edited by Bobae 2025
%--------------------------------------------------------------------------

function [param,xxx,yyy,opttheoxxx,opttheoyyy] = zo_intp(bdx,bdy,filepath,cctime)
    t1 = datetime('now');
    % theory of curve generated rightside up
    updown = 1; 
    bdx_cm = bdx; bdy_cm = bdy; 
    % num of interpolation pt
    numpt = 1000; 
    [~,ia] = unique(bdx_cm);
    bdx_cm = bdx_cm(ia); bdy_cm = bdy_cm(ia); 
    % linearly interpolation result
    xxx = linspace(min(bdx_cm),max(bdx_cm),numpt); 
    yyy = interp1(bdx_cm,bdy_cm,xxx);
    % height of current pinnacle
    hei = max(yyy)-min(yyy); 
    maxval = [max(yyy)-0.001,max(yyy)+0.001];
    % might have multiple values if the tip is unclear/blunted
    ind = (yyy > maxval(1) & yyy < maxval(2));
    shiftx = xxx(ind);
    if length(shiftx) > 1
        shiftx = mean(shiftx); % peak mean value
    end
    xxx = xxx-shiftx; % x-axis centered at 0
    bdx_cm = bdx_cm-shiftx; 

    index = 1; 
    if index == 1
        % brute force, param setup
        % r accuracy: 0.0005cm = 5 microns
        % radius range
        rrrange = linspace(0.015,0.045,61); % Previous: rrrange = linspace(0.0005,0.03,60)
        % x-y accuracy: 0.05cm 
        % vertical shift range
        vhrange = -1:.01:2; % 101 values, yspacing of 0.05cm. Previous: (4,20,321);
        % horizontal shift
        xhrange = linspace(-.1,.1,41); % spacing of 0.05cm. Previous: (-.1,.1,41), changedto (-1,1,41)...
        % xhrange = zeros(1,1);
        [m,n,x] = ndgrid(rrrange,vhrange,xhrange);
        % [m,n] = ndgrid(rrrange,vhrange);
        Z = [m(:),n(:),x(:)];
        mindistset = zeros(length(Z),1);
        
        % parallel computing
        parfor  (i=1:length(Z),20)
            param = Z(i,:);
            [fxx,fyy,dist] = attractor3d(param,hei,numpt,updown,xxx,yyy);
            mindistset(i) = dist;
        end
    
        % choose the optimal param set
        optind = (mindistset == min(mindistset));
        assert(sum(optind) == 1);
        param = Z(optind,:);
        [opttheoxxx,opttheoyyy,~] = attractor3d(param,hei,numpt,updown,xxx,yyy);

    elseif index == 2
        x0 = [0.01,5,0.1];
        f = @(x)attractor3d(x0,hei,numpt,updown,xxx,yyy);
        [x,fval,exitflag,output,grad] = fminunc(f,x0);

    end
end




