function [] = nBodyPlot(dx, nb)
% Plots nb iBodies and their supports.
% nBodyPlot(dx, nb, d)
% Physical spacing dx allows distances to be converted to lattice units.
% d is the number of dimensions

close all

% Plotting flags
cols = {'r', 'k', 'b', 'g', 'm','c'};
style = {'rx','-ko','b.','gs','md','c*'};

for n = 1:nb
    
    % Read body coordinates
    Body = dlmread(['../Output/IBbody_' num2str(n-1) '.out'],'\t',1,0);
    
    if (n == 1)
        figure
        subplot(1,2,1);
    end
    
    % Plot just the markers
    for i = 1:size(Body,1)
        plot3(Body(i,1)/dx, Body(i,2)/dx, Body(i,3)/dx, [cell2mat(cols(mod(i,length(cols))+1)) '^'],'MarkerSize',10)
        hold on
    end
    view(d)
    axis tight
    axis equal
    grid on
    
end


for n = 1:nb

    % Read body coordinates
    Body = dlmread(['./Output/IBbody_' num2str(n-1) '.out'],'\t',1,0);

    % Loop over markers and read in support
    for i = 0:size(Body,1)-1
        eval(['Supp_' num2str(i) ' = dlmread(''../Output/Supp_' num2str(n-1) ...
            '_' num2str(i) '.out'',''\t'',1,0);'])
    end

    % Plot markers and support
    if (n == 1)
        subplot(1,2,2);
    end
    for i = 1:size(Body,1)
        plot3(Body(i,1)/dx, Body(i,2)/dx, Body(i,3)/dx, [cell2mat(cols(mod(i,length(cols))+1)) '^'],'MarkerSize',10)
        grid on
        hold on
        eval(['plot3(Supp_' num2str(i-1) '(:,1)/dx, Supp_' num2str(i-1) '(:,2)/dx, Supp_' num2str(i-1) '(:,3)/dx,'''...
            cell2mat(style(mod(i,length(style))+1)) ''',''MarkerSize'',10)'])    
    end
    view(d)
    axis tight
    axis equal
end
