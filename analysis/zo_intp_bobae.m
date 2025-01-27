%--------------------------------------------------------------------------
% Optimization function to the family of curve (Huang & Moore, 2022) in zoom-out view
%
% Steven Zhang, Courant Institute
% Updated Jan 2023
%--------------------------------------------------------------------------
%exp_x, exp_y {j}: spaced linearly interpolated simulation results
%param_cells{j}: [r verticalshift]
%x_cells, y_cells {j}: theoretical best fit
function [exp_x,exp_y,param_cells,x_cells,y_cells] = zo_intp_bobae(Xsave,spacing,h0,epsilon)
    [x,y] = cutoff_tails(Xsave,spacing,h0,epsilon); %cutoff extra tails and reduce number of frames
    
%     t1 = datetime('now');
    updown = 2; % theory of curve generated upside down
    
    
    numpt = 1000;
    
    %For each frame:
    for j = 1:length(x) 
        %linearly interpolate
        xj = x{j}; yj = y{j}; 
        yj = yj - min(yj); %recenter to be sitting on x-axis.
        xx = linspace(min(xj),max(xj),numpt); % linearly interpolation result
        yy = interp1(xj,yj,xx);
        hei = max(yy)-min(yy); % height of current pinnacle 
        exp_x{j} = xj; exp_y{j} = yj;
        
        %NEW STUFF
        rrrange = linspace(0.005,0.025,100); % radius range
        % x-y accuracy: 0.05cm 
        vhrange = linspace(-.25,.25,50); % vertical shift range
        %Deleting the below for now because in terms of simulation, we know
        %we are always aligned horizontally.
%         
%         xhrange = linspace(-0.1,0.1,41); % horizontal shift
%         [m,n,x] = ndgrid(rrrange,vhrange,xhrange); %length(x)
%         Z = [m(:),n(:),x(:)]; length(Z)

        Z = rrrange;
        [m,n] = ndgrid(rrrange,vhrange);
        Z = [m(:),n(:)];size(Z);
        mindistset = zeros(length(Z),1);
        
%         % parallel computing
        parfor (i=1:length(Z),2000)
            param = Z(i,:);
            [fxx,fyy,dist] = attractor3d(param,hei,numpt,updown,xx,yy);
            mindistset(i) = dist;
        end
%         'hi'
%         % choose the optimal param set
        optind = (mindistset == min(mindistset));
%         assert(sum(optind) == 1);
        param = Z(optind,:);
%         param = Z(optind,:);
        [opttheoxxx,opttheoyyy,~] = attractor3d(param,hei,numpt,updown,xx,yy);
%         param
%         hei
        j
        param_cells{j} = param;
        x_cells{j} = opttheoxxx;
        y_cells{j} = opttheoyyy;
    end
end




