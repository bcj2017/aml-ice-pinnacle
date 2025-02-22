close all

path = 'curve-in-data/';
S = dir(fullfile(path,'**','*.mat'));
names = {S.name};
figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on
for i = 1:length(names)
    name = names{i};
    data = load([path,name]);
    time = data.ttt;
    curv_rad = data.curv_rad;
    plot(time,curv_rad,'-o','color', [.5 .5 .5],'LineWidth',2)
end
% yline(0.05,'--','LineWidth',2)
yline(0.01,'--','LineWidth',2)
yline(0.06,'--','LineWidth',2)
xline(1275,'--','LineWidth',2)
xlabel('Time (s)','FontSize',14)
ylabel('Radius of Curvature (cm)','FontSize',14)
% legend([names],'Location','Northeast') %Bobae comment out
title('Zoom-In Radius of Curvature','FontSize',14)
%ylim([0 0.15])
% xlim([0,800]); ylim([0,.2])%Bobae
% xlim([0 1200])
% set(gca, 'YScale', 'log')


%% Bobae
% legend({'2023-02-28-a','2023-03-01-a','2023-03-16-b','2023-03-16-c','2023-03-16-d','2023-03-16-e','2023-04-05-a','2023-04-05-b','2023-04-05-c','2023-04-05-d'},'Location','Northeast');




