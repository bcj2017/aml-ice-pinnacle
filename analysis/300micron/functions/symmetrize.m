function [peak_x,peak_z,peakbox] = symmetrize(peak_x,peak_z,epsilon)
    while abs(peak_z(1)-peak_z(end)) > epsilon
        if peak_z(1) > peak_z(end)
            peak_x = peak_x(1:end-1);
            peak_z = peak_z(1:end-1);
        else
            peak_x = peak_x(2:end);
            peak_z = peak_z(2:end);
        end
    end
    peakbox = (max(peak_x)-min(peak_x))/2;
    disp(['Corrected peakbox size: ',num2str(peakbox)])
end
