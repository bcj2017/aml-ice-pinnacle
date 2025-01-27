% Find peak
function [peak_x,peak_z,xpeakloc,zpeakloc] = isolatePeakPoints(x_pxl,z_pxl,peakbox)
    [~,peak_idx] = max(z_pxl);
    xpeakloc = x_pxl(peak_idx); zpeakloc = z_pxl(peak_idx);
    
    peak_x = []; peak_z = [];
    for i = 1:length(x_pxl)
        if x_pxl(i) > xpeakloc-peakbox && x_pxl(i) < xpeakloc+peakbox && ...
            z_pxl(i) > zpeakloc-peakbox && z_pxl(i) < zpeakloc+peakbox
            peak_x(end+1) = x_pxl(i);
            peak_z(end+1) = z_pxl(i);
        end
    end
end