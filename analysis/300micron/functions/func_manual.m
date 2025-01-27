%--------------------------------------------------------------------------
% Manual operation on boundary data points
% Based on ginput2.m 
% Delete unwanted points -> Add Points -> Check
% Use mouse to operate: 
% Single left click: zoomin; 
% Double click: back to original view
% Right click: select datapoint and mark
% [The frame will be automatically centered based on the direction on the last two added
% points for better user experience]
% Left click (outside of frame) + Enter: stop and move on
% (Be careful otherwise may add unwanted point :()
% [Unfortunately no rollback option is provided internally in this function, so please be
% as possible in operations]
% Each finished frame will be required to be saved (or not) outside this function, so if
% error encountered, the user could change the starting frame and rerun the main function
% 
% ==== 20240702 Update ====
% New methodology for deleting points by using 4 points to construct boundaries and
% deleting all points inside of it (accept for multiple iterations). 
%
% Zihan Zhang, Courant Institute
% Updated July 2024
%--------------------------------------------------------------------------

function [thisx,thisy] = func_manual(pic,datax,datay,fr)
% INPUT
% pic - background image
% datax, datay - extracted datapoints using alphashape-algo
% OUTPUT
% thisx, thisy - edited datapoints

    close all
    figure()
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    %% Delete Points
    imshow(pic); hold on; title('Deleted Points');
    scatter(datax,datay,'green','filled');

    % Delete points by region
    while true
        rectx = [];
        recty = [];
        
        for i = 1:4
            [x, y] = ginput2(1,'PlotOpt');
            rectx = [rectx; x];
            recty = [recty; y];
            scatter(rectx, recty, 'red', 'filled');
            if length(rectx) > 1
                plot(rectx, recty, 'r-');
            end
        end
        
        plot([rectx; rectx(1)], [recty; recty(1)], 'r-');
        
        min_x = min(rectx);
        max_x = max(rectx);
        min_y = min(recty);
        max_y = max(recty);
        
        in = datax >= min_x & datax <= max_x & datay >= min_y & datay <= max_y;
        datax(in) = [];
        datay(in) = [];
        
        hold on; title('Delete Points');
        disp(['Number of points in region: ', num2str(sum(in))]);
        scatter(datax, datay, 'green', 'filled');
        
        cont = input('Continue deleting? ');
        if cont == 0
            break;
        end
    end

    close


    %% Add Points
    figure()
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    imshow(pic); hold on; title('Added Points');
    scatter(datax,datay,'green','filled');
    [addx,addy] = ginput2('PlotOpt'); 
    datax = datax'; datay = datay'; % Bobae: added because of a concatenation error
    disp(size(addx));
    disp(size(datax));
    thisx = [datax;addx];
    thisy = [datay;addy];
    disp(['Added ', num2str(length(addx)), ' Points!']);
    close all

    %% CheckPoint
    figure()
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    imshow(pic); hold on; title('Final Result')
    scatter(thisx,thisy,'blue','filled');
    pause
    close all

    [thisx,thisy] = fr.angle_sort(thisx,thisy);
end





