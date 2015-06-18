clear
close all

% Read body coordinates
Body = dlmread('IBbody.out','\t',1,0);

% Loop over markers and read in support
for i = 0:size(Body,1)-1
    eval(['Supp_' num2str(i) ' = dlmread(''Supp_' num2str(i) '.out'',''\t'',1,0);'])
end

% Plot markers and support
figure;
cols = {'r', 'k', 'b', 'g', 'm','c'};
style = {'rx','-ko','b.','gs','md','c*'};
for i = 1:size(Body,1)
    plot3(Body(i,1), Body(i,2), Body(i,3), [cell2mat(cols(mod(i,length(cols))+1)) '^'],'MarkerSize',10)
    if all(Body(1,3) == Body(:,3)) % If 2D then look top down
        view(2)
    end
    grid on
    hold on
    eval(['plot3(Supp_' num2str(i-1) '(:,1), Supp_' num2str(i-1) '(:,2), Supp_' num2str(i-1) '(:,3),'''...
        cell2mat(style(mod(i,length(style))+1)) ''',''MarkerSize'',10)'])    
    % disp('Press a key for next marker')
    % pause
end
axis equal

figure
for i = 1:size(Body,1)
    plot3(Body(i,1), Body(i,2), Body(i,3), [cell2mat(cols(mod(i,length(cols))+1)) '^'],'MarkerSize',10)
    hold on
end
axis equal
grid on
