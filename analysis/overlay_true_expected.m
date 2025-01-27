%% Plot overlaid true and expected shapes.
% Given param_cells, x_cells, y_cells
%param_cells{j} = [r verticalshift]
%x_cells, y_cells are theoretical best fit

close all; figure();
hold on

for j=1:length(x_cells)
    plot(exp_x{j},exp_y{j},'Linewidth',1.2,'Color','r')
    plot(x_cells{j},y_cells{j},'Linewidth',1.2,'Color','b')
    axis equal
end
legend('Simulation Boundary','Theoretical Best Fit')