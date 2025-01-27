close all; clear;

R = 0.03; %radius of curvature, (cm)
H = 10; %total height, (cm)
nnpt = 1000; % number of points

% t0 = asin(sqrt((R-sqrt(R^2 - 4*(R/4-H)*(3*R/4)))/(2*(R/4-H)))); %solve for final tangent angle
% theta = linspace(pi/2,t0,nnpt/2); %setup grid

t0 = asin(sqrt((R-sqrt(R^2 - 4*(R/4-H)*(3*R/4)))/(2*(R/4-H)))); %solve for final tangent angle
theta = linspace(pi/2,t0,nnpt); %setup grid

r = @(theta) R*(cos(theta)./sin(theta).^3); %r coordinate
z = @(theta) H-R*(-1./sin(theta).^2 + (3/4)*1./sin(theta).^4 + (1/4)); %z coordinate

xx = [-fliplr(r(theta)) r(theta)];
yy = [fliplr(z(theta)) z(theta)];

% flip = 1;
% if flip == 2
%     yy = -yy+H;
% end

[~,ind] = unique(xx);
xx = xx(ind); yy = yy(ind);

x_pinnacle = xx; y_pinnacle = yy; % equispaced in theta

% == redistribute points to be equispaced on x axis == %
% x_pinnacle = linspace(min(xx),max(xx),nnpt);
% y_pinnacle = interp1(xx,yy,x_pinnacle);

figure();
p = plot(x_pinnacle,y_pinnacle);%,'.');
p.LineWidth = 1;
ylim([0,H+1]);
axis equal;
xlim([5,5]);