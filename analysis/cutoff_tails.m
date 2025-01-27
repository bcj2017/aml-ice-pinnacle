%% This function takes in Xsave, a frame interval, and some tolerance epsilon,
% then cuts off the tails produced by the phase-field simulations
% not allowing melting at the base.

function [xxx,yyy] = cutoff_tails(Xsave,spacing,h0,epsilon)

N = length(Xsave);
L = h0*1e-2 / .8;
x0 = Xsave{1}(1,:); y0 = Xsave{1}(2,:);
x0 = x0 * L * 1e2; y0 = y0 * L * 1e2; %convert to cm
% epsilon = 3e-2; %some tolerance level

xxx{1} = x0; yyy{1} = y0;
for j = 1:spacing:N
    x1 = Xsave{j}(1,:); y1 = Xsave{j}(2,:); %get data from Xsave
    x1 = x1 * L * 1e2; y1 = y1 * L * 1e2; %convert everything to cm
    
%     peakboxrange = .2202; %in cm, approximately 600 pixels (Steven's peakboxrange)
    
    %Deleting bottom bits
    if j == spacing+1 %second in iteration
        %Try deleting manually for this one
        logicalIndices = (y1 > .25); %Delete bottom 2 centimeters %HAVE TO DECIDE MANUALLY!
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
    
    x0 = x1; y0 = y1; %have previous pinnacle stored
    
    xxx{(j-1)/spacing+2} = x1; yyy{(j-1)/spacing + 2} = y1;
end
end