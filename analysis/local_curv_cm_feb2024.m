%% Compute local radius of curvature, converting everything to dimensions first.
% close all;
% Let n be the number of points used in interpolation.
% ^Renormalize to get R = 0.01 in first data set
close all; clear rk_arr;
L = .08 / .8;
N = length(Xsave); rk_arr = zeros(1,N); %spacing = floor(145/37);
for j = 1:N
    x1 = Xsave{j}(1,:); y1 = Xsave{j}(2,:); %get data from Xsave
    x1 = x1 * L * 1e2; y1 = y1 * L * 1e2; %convert everything to cm
    
    peakboxrange = .2202; %in cm, approximately 600 pixels (Steven's peakboxrange)
    peakboxrange = .04; %ugh trying smaller
    if j == 1
        n = 30; %just for a better fit
        xred = x1(length(x1)/2-n/2+1:length(x1)/2+n/2); yred = y1(length(y1)/2-n/2+1:length(y1)/2+n/2);
    else
        logicalIndices = (x1 > -peakboxrange) & (x1 < peakboxrange); %take values in peakboxrange
        xred = x1(logicalIndices);
        yred = y1(logicalIndices);
    end
    
    
    %Fourth-order polynomial interpolation, from Scott / Steven:
    p = polyfit(xred,yred,4);
    x = linspace(min(xred),max(xred),100); 

    rr1 = linspace(min(xred),max(xred),10000);
    rr2 = linspace(min(yred),max(yred),10000);
    opt_deg = 4;
    px = p(opt_deg); pxx = 0;
    for m=1:opt_deg-1
        px = px+(opt_deg-m+1)*p(m)*rr1.^(opt_deg-m);
        pxx = pxx+(opt_deg-m+1)*(opt_deg-m)*p(m)*rr1.^(opt_deg-m-1);
    end
    k0 = max(abs(pxx)./(1+px.^2).^(3/2)); % curvature
    rk = 1/k0; % radius of curvature in cm
%     disp(['box size: ', num2str(smallran), ...
%         newline 'radius of curvature: ',num2str(rk)])
    % rk_lst(end+1) = rk;
    rk_arr(j) = rk;
    if mod(j,50) == 1
        figure();
        plot(x1,y1,'.','MarkerSize',10); hold on
        xxx = linspace(min(xred),max(xred),500);
        plot(xxx,polyval(p,xxx),'Linewidth',1.2)
        
        xmin = 2*min(xred); xmax = 2*max(xred);
        xline(-peakboxrange); xline(peakboxrange)
        axis equal
        
        title(strcat('Frame',{' '},string(j),', Radius of Curvature: ',{' '},string(rk), {' '},'cm' )); xlabel('cm'); ylabel('cm') 
        xlim([-peakboxrange*1.2,peakboxrange*1.2]);
        legend({'Simulation Boundary','Fourth-Order Polynomial Fit','Peakboxrange'},'Location','South')
    end
end
figure();
plot(linspace(0,1,N),rk_arr(1:N),'*','MarkerSize',10)
hold on
title('Computed radius of curvature in cm')
xlabel('time'); ylabel('Local Radius of Curvature');
yline(0.01,'b','Linewidth',1.2); yline(0.06,'r','Linewidth',1.2);
xline(3/4,'Linewidth',1.2)
legend({'Measured in Simulation','Expected Theory','Expected Experimental Attractor','3/4 Time'});
total_avg = mean(rk_arr)
middle_avg = mean(rk_arr(floor(N/4):ceil(3*N/4)))