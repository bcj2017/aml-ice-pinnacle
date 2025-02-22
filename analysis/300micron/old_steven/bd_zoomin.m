%--------------------------------------------------------------------------
% Main function of data analysis of pinnacle experiment (zoomed in view)
%
% Steven Zhang, Courant Institute
% Updated Jul 2023
% Modified 2025 by Bobae
%--------------------------------------------------------------------------

close all
clear
warning('off')

addpath('functions');
% addpath('functions');

% declare analyze object 
basepath = '../../../experiments/300micron/';
filepath = '2025-01-20-bubblyice'; %experiment name
basepath = [basepath,filepath];
subfolder = '/zoomin/nicest/'; % NOTE: remove 'nicest' for regular experiments
conv_index = 1; 

fc = func_curve();

convratio = calibra(basepath,subfolder,1); % run calibration function

path = [basepath,subfolder];
S = dir(fullfile(path,'**','*.JPG'));
names = {S.name};

% functionality index
overlay = 0;
inorout = 0;
testdemo = 0;
testrange = 1;

% define constants, iterative object, and other data structures
boxnumcount = 50;
boxstabthreshold = 0.02;
% [HYPER]
% peak box half size [~601.6]
% peakboxrange = linspace(100,1500,boxnumcount); % Bobae: comment out, uncomment 43 to be 601.6/2
peakboxrange = 400;[601.6/2];
dt = 10; % s
time_register = []; % time register
curv_register = []; % curvature register
arclength_register = []; 

localcurv_register = []; % local curvature register
callen = floor((length(names)+5)/10);
boundCell = cell(1,callen); % save tracking result
boxidxCell = cell(1,callen);
boxidxMatchValue = [];

% main loop
cnt = 1;
for jjj = 1:length(names)
%     if jjj < 20 || mod(jjj,5) == 0 % used for calculation
    if mod(jjj,5) == 1 % used for plotting
        % === User Manual === 
        % 0: set as 0; 2: run hist; 3: normal; 
        % 4: pass; other value: change sparse point threshold
        % === End === 
        disp("*************************")
        disp("=== Start Calculation ===")
        disp("*************************")
        execind = 3; 
        % while execind ~= 4 && execind ~= 0
            % try
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

            % figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
            % subplot(1,6,1)
            [xx,idx] = sort(xx); yy = yy(idx);
            % plot(xx,yy,'o','MarkerSize',2); axis equal
            % title('Original')

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

            % Bobae added manual read stuff:
            % spimage = char(names(jjj));
            % disp(['Image name: ',spimage]) % prints i.e. DSC_3874.JPG
            % pic = imread([basepath,subfolder,spimage]);
            fr = func_aux();
            [newbdx,newbdy] = func_manual(pic,ex,ey,fr);
            % Bobae end
                
            subplot(1,5,1)
            [newbdx,idx] = sort(newbdx);
            newbdy = newbdy(idx);
            plot(newbdx,newbdy,'o','MarkerSize',2); axis equal
            title('Raw boundary points')

            % find X and Y max gap; box size larger than this value is 
            maxDiffX = max(newbdx) - min(newbdx);
            maxDiffY = max(newbdy) - min(newbdy); 
            takeVal = min(maxDiffX, maxDiffY);
            if jjj == 1
                peakboxrange = linspace(100,takeVal,boxnumcount);
            end

            % find peak
            prr = [0.9*mean(newbdx),1.1*mean(newbdx)];
            indr = (newbdx >= prr(1) & newbdx <= prr(2));
            kxx = newbdx(indr); kyy = newbdy(indr);
            [~,ind] = min(kyy);
            r_peak_x = kxx(ind); r_peak_y = kyy(ind);

            % different box size
            rk_lst = []; 
            for iii = 1:length(peakboxrange)
                smallran = peakboxrange(iii);

                peak_x = []; peak_y = [];
                for i = 1:length(newbdx)
                    if newbdx(i) > r_peak_x-smallran && newbdx(i) < r_peak_x+smallran && ...
                        newbdy(i) > r_peak_y-smallran && newbdy(i) < r_peak_y+smallran
                        peak_x(end+1) = newbdx(i);
                        peak_y(end+1) = newbdy(i);
                    end
                end

                
                
                subplot(1,5,2)
                plot(peak_x,peak_y,'o','MarkerSize',2); axis equal
                hold on
                plot(r_peak_x,r_peak_y,'*','MarkerSize',2)
                % draw indicator rectangle
                if (iii == length(peakboxrange))
                    rectangle('Position',[r_peak_x-smallran,r_peak_y-smallran, 2*smallran, ...
                        2*smallran]);
                end
                title('Raw boundary points in peakbox')

                % change unit
                peak_x = peak_x * convratio; 
                peak_y = peak_y * convratio; 
                
                rr1 = linspace(min(peak_x),max(peak_x),10000);
                rr2 = linspace(min(peak_y),max(peak_y),10000);
                % polynomial degree check
                [p,opt_deg] = degree_check(peak_x,peak_y,rr1,16);
                zz = polyval(p,rr1);
                
                plot(rr1/convratio,zz/convratio)
                title(num2str(opt_deg))
                hold off
                
                % compute tip curvature (from S.W.)
                px = p(opt_deg); pxx = 0;
                for m=1:opt_deg-1
                    px = px+(opt_deg-m+1)*p(m)*rr1.^(opt_deg-m);
                    pxx = pxx+(opt_deg-m+1)*(opt_deg-m)*p(m)*rr1.^(opt_deg-m-1);
                end
                % curvature
                k0 = max(abs(pxx)./(1+px.^2).^(3/2)); 
                % radius of curvature
                rk = 1/k0; 
                disp(['box size: ', num2str(smallran), ...
                    newline 'radius of curvature: ',num2str(rk)])
                rk_lst(end+1) = rk;

                % compute local curvature (sensitive to tracking quality)
                % find arclength
                N = length(peak_x);
                ds = sqrt((peak_x(2:end)-peak_x(1:end-1)).^2+...
                    (peak_y(2:end)-peak_y(1:end-1)).^2);
                % arclength accumulation
                s = zeros(N,1); 
                for n = 1:N-1
                    s(n+1) = s(n) + ds(n);
                end
                finalL = s(end);
                % equidistributed arclength
                xxss = linspace(0,finalL,4096);
                % polynomial interpolation
                xys_fitdeg = 8;
                [p1,optdeg1] = degree_check(s,peak_x,xxss,xys_fitdeg);
                [p2,optdeg2] = degree_check(s,peak_y,xxss,xys_fitdeg);
                % % Bobae: seems redundant?
            end

            
            subplot(1,5,3)
            plot(peakboxrange,rk_lst)
            xlabel('Box Length / 2 (pixel)')
            ylabel('Radius of Curvature (cm)')
            title('BoxLength Sensitivity')
            [minidx,maxidx] = find_longest_interval([peakboxrange;rk_lst]',boxstabthreshold);
            % visually 
            xline(peakboxrange(minidx),'LineStyle','--');
            xline(peakboxrange(maxidx),'LineStyle','--');

            ax5 = subplot(1,5,4);
            xxx = polyval(p1,xxss); % x = f(s)
            yyy = polyval(p2,xxss); % y = g(s)
            plot(xxx,yyy,'*');
            axis equal; hold on;
            plot(peak_x,peak_y,'o')
            title(['X/Y as f(s) ---',num2str(optdeg1)]);% ,';',num2str(optdeg2),'.'])

            % https://en.wikipedia.org/wiki/Curvature
            p1_fd = polyder(p1); p1_sd = polyder(p1_fd);
            p2_fd = polyder(p2); p2_sd = polyder(p2_fd);
            res1 = abs(polyval(p1_fd,xxss).*polyval(p2_sd,xxss)...
                -polyval(p2_fd,xxss).*polyval(p1_sd,xxss));
            res2 = (polyval(p1_fd,xxss).^2+polyval(p2_fd,xxss).^2).^(3/2);
            cuv = res1./res2; % curvature

            subplot(1,5,5)
            norxxss = normalize(xxss);
            plot(norxxss,cuv)
            bd1 = -1; bd2 = 1;
            xlim([bd1 bd2])
            secind = (norxxss > bd1) & (norxxss < bd2);
            leftpos = sum((norxxss <= bd1)); 
            cuvsub = cuv(secind); xxssxub = norxxss(secind);
            [~,maxcuvind] = max(cuvsub); 
            maxcuvindkk = maxcuvind + leftpos;
            hold on
            plot(norxxss(maxcuvindkk),cuv(maxcuvindkk),'*','MarkerSize',5)
            xlabel('Normalized Cumulative Arclength')
            ylabel('Curvature')

            plot(ax5,xxx(maxcuvindkk),yyy(maxcuvindkk),'o','MarkerSize',5,...
                'MarkerFaceColor','green')
            legend(ax5,'Intp Result','Rough Result','Curviest Point')

            % another approach
            [L,R,K] = fc.curvature([polyval(p1,xxss);polyval(p2,xxss)]');
           
            execind = input('Manual Check Result = ');
            close all
            % catch exception
            %     msgText = getReport(exception);
            %     disp('PROBLEM ENCOUNTERED!')
            %     disp(msgText)

            % execind = input('Manual Check Result = ');
            % close all
            % end
        % end

        % edge case
        % if execind == 0
        %     rk_lst = zeros(1,length(peakboxrange)); 
        % end

        % data registeration
        boundCell{cnt} = [peak_x;peak_y];
        time_register(end+1) = jjj*dt;
        curv_register(end+1,:) = rk_lst; 

        arclength_register(end+1,:) = norxxss;
        localcurv_register(end+1,:) = cuv;

        % box range iteration
        stableRange = peakboxrange(minidx:maxidx);
        boxidxCell{cnt} = stableRange;
        matchedValue = mean(rk_lst(minidx:maxidx));
        boxidxMatchValue(end+1) = matchedValue;

        cnt = cnt + 1;

    end
end


% plot the radius of curvature vs. t [when not overlaying]
checkind = (curv_register < 0.02);
curv_register(checkind) = [];
time_register(checkind) = [];
% IGNORE BECAUSE overlay = 0
if overlay ~= 1 && testdemo == 1
    figure()
    plot(time_register,curv_register,'-o');
    xlabel('Time (s)','FontSize',14)
    ylabel('Radius of Curvature (cm)','FontSize',14)
    title(['Length of Peak Box = ', num2str(2*smallran)],'FontSize',14)
    
    if ext_stra == 1
        ext_method = 'whole';
    elseif ext_stra == 2
        ext_method = 'xmost';
    end
    ttt = time_register;
    curv_rad = curv_register;
    save(['curve-in-data/',filepath,'-',ext_method],'ttt','curv_rad')
end

% IGNORE BECAUSE OVERLAY = 0
if overlay == 1 && testdemo == 1
    usedboundcoll = cell(1,length(checkind)-sum(checkind));
    for i = 1:length(checkind)
        if checkind(i) == 0
            usedboundcoll{length(checkind(1:i))-sum(checkind(1:i))} = boundCell{i};
        end
    end
    overlay_func(usedboundcoll,overlay,time_register,filepath(1:end),inorout)
end

localcurv_ind = 1;
if localcurv_ind == 1 && testdemo == 0
    datasize = size(localcurv_register);
    Ntot = datasize(1);
    figure()
    hold on
    xlim([bd1/2 bd2/2])
    for i = 1:Ntot
        lcsel = localcurv_register(i,:);
        arcsel = arclength_register(i,:);
        ssrr = arcsel > bd1 & arcsel < bd2;
        [maxval,maxind] = max(lcsel(ssrr));
        partarcsel = arcsel(ssrr);
        sftarc = partarcsel(maxind);
        plot(arcsel-sftarc,lcsel,'LineWidth',3,'Color',rb_color(i,Ntot))
    end
    timemarkerCell = cell(1,Ntot);
    for iN = 1:Ntot
        timemarkerCell{iN} = num2str(time_register(iN),'N=%-d s');
    end
    legend(timemarkerCell,'FontSize',15);
    title(['Max Order of Poly Fitting=',num2str(xys_fitdeg),';Box Size=',num2str(peakboxrange(1))], ...
        'FontSize',14)
    xlabel('Shifted Normalized Cumulative Arclength (z-score)','FontSize',14)
    ylabel('Curvature (cm^{-1})','FontSize',14)
end

% test stability under various box length
assert(length(time_register) == length(boxidxCell));

if testrange == 1
    figure(1)
    hold on
    for iii = 1:12
        boxrr = boxidxCell{iii};
        tt = time_register(iii);
        plot([tt tt],[boxrr(1) boxrr(end)],'LineStyle','-.','LineWidth',2);
    end
    for iii = 2:12
        oldtt = time_register(iii-1);
        curtt = time_register(iii);
        oldrr = boxidxCell{iii-1};
        currr = boxidxCell{iii};
        oldtop = [oldtt,oldrr(end)];
        oldbot = [oldtt,oldrr(1)]; 
        curtop = [curtt,currr(end)];
        curbot = [curtt,currr(1)];
        plot([oldtop(1) curtop(1)],[oldtop(2) curtop(2)],'LineStyle','--','LineWidth',2,'Color','red')
        plot([oldbot(1) curbot(1)],[oldbot(2) curbot(2)],'LineStyle','--','LineWidth',2,'Color','red')
    end
    title(['Box Half Length Stability over Time, threshold=',num2str(boxstabthreshold)]);
    xlabel('Time')
    ylabel('Half-Length of Box')
    hold off

    
    figure(2)
    plot(time_register,boxidxMatchValue,'LineStyle','-.')
    title("Average Radius of Curvature calculated in the stability range")
    xlabel('Time')
    ylabel('Average R')

end




















