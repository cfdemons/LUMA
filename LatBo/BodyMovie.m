function [] = BodyMovie(bod, n, out)
% Generate animation of IBbody marker movement.
% bod is the number of bodies. n is the number of timesteps to be
% included in the movie. out is how timesteps were performed before writing
% out (see defintions.h).

% Create movie object
Vid = VideoWriter('./Output/BodyMovie.avi');
% Vid.FrameRate = 25;
open(Vid);
exit_flag = false;
xy0 = 0;
close all

% Cycle through number of text files and read in positions
for c = 1:n+1
    
    if mod(c,out) == 0 % Only check for files which we know exist
        
        for b = 0:bod-1
            try % If simulation crashes still allows movie to be built up to that point
            eval(['xy' num2str(b) ' = csvread(''./Output/Body_' num2str(b)... 
                '_position_' num2str(c) '.out'',1,0);']) % Read in values
            catch
                exit_flag = true;
                break                
            end
        end
        if exit_flag == true
            break
        end
        
        % Plot positions
        figure;
        hold on;
        for b = 0:bod-1
            eval(['plot3(xy' num2str(b) '(:,1), xy' num2str(b) '(:,2), xy' num2str(b) '(:,3),''^'');'])
        end
        axis([1.4 2 .1 .9 .1 .9]);
        view(3)
        grid
        
        % Add frame to movie
        set(gcf,'Renderer','zbuffer');
        writeVideo(Vid,getframe);
        close;
    end
end