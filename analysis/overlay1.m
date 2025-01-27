%% Overlay Plots

% Every spacing = 5 frames, overlays pinnacle boundaries aligned at the tip
%% Compute local radius of curvature, converting everything to dimensions first.
% close all;
% Let n be the number of points used in interpolation.
% ^Renormalize to get R = 0.01 in first data set
close all; 
N = length(Xsave);
spacing = 100;
figure(); hold on

%initial pinnacle
h0 = 8; %cm
L = h0*1e-2 / .8;
x0 = Xsave{1}(1,:); y0 = Xsave{1}(2,:);
x0 = x0 * L * 1e2; y0 = y0 * L * 1e2; %convert to cm
epsilon = 3e-2; %some tolerance level
for j = 1:spacing:N
    x1 = Xsave{j}(1,:); y1 = Xsave{j}(2,:); %get data from Xsave
    x1 = x1 * L * 1e2; y1 = y1 * L * 1e2; %convert everything to cm
    
%     peakboxrange = .2202; %in cm, approximately 600 pixels (Steven's peakboxrange)
    
    %Deleting bottom bits
    if j == spacing+1 %second in iteration
        %Try deleting manually for this one
        logicalIndices = (y1 > .75); %Delete bottom 2 centimeters %HAVE TO DECIDE MANUALLY!
        x1 = x1(logicalIndices);
        y1 = y1(logicalIndices);
    end
    
    %Move pinnacle to be at same initial height.
    ymax = max(y1); y1 = y1 + h0 - ymax;
    
    if j > spacing+1 %after second iteration
        %Delete extraneous bits based on previous pinnacle
        %align centers of pinnacles:
        len0 = length(y0) / 2; %half-length of previous pinnacle
        len1 = length(y1) / 2; %half-length of current pinnacle
        %Shorten current pinnacle (old pinnacle is shorter because we deleted some already)
        if len0 > len1
            x0 = x0(len0-len1+1: len0+len1);
            y0 = y0(len0-len1+1: len0+len1);
        else
            x1 = x1(len1-len0+1: len0+len1);
            y1 = y1(len1-len0+1: len0+len1); 
        end
        %length(y1)-length(y0) %sanity check that we did this right
        logicalIndices = (abs(y1-y0) < epsilon);
        x1 = x1(logicalIndices); y1 = y1(logicalIndices);
        
    end
    
    frame_legend = strcat('Frame',{' '},num2str(j));
    plot(x1,y1,'.','DisplayName',frame_legend{1})
    
    
    xlabel('cm'); ylabel('cm')
    title(strcat('Pinnacles Overlayed at Tip Height of', {' '}, num2str(h0), {' '},'cm, Deleted Portions of Tolerance',{' '},num2str(epsilon),{' '},'cm'))
    ylim([0,h0+1]);
    axis equal
    
    x0 = x1; y0 = y1; %have previous pinnacle stored
    
    lgd_entry{j} = frame_legend{1};
end
legend()