function [] = JacMovie(bod, n, out)
% Generate animation of filament marker movement.
% bod is the number of filaments. n is the number of timesteps to be
% included in the movie. out is how timesteps were performed before writing
% out (see defintions.h).

% Create movie object
Vid = VideoWriter('./Output/Movie.avi');
% Vid.FrameRate = 25;
open(Vid);
exit_flag = false;

% Cycle through number of text files and read in positions
for c = 1:n+1
    
    if mod(c,out) == 0 % Only check for files which we know exist
        
        for b = 0:bod-1
            try % If simulation crashes still allows movie to be built up to that point
            eval(['xy' num2str(c) ' = csvread(''./Output/Jac_' num2str(b)... 
                '_position_' num2str(c) '.out'',1,0);'])
            catch
                exit_flag = true;
                break                
            end
        end
        if exit_flag == true
            break
        end
        
        % Plot positions
        figure('Visible','off');
        eval(['plot(xy' num2str(c) '(:,1), xy' num2str(c) '(:,2),''^'');'])
        axis([1.4 2 .1 .9]);
        
        % Add frame to movie
        set(gcf,'Renderer','zbuffer');
        writeVideo(Vid,getframe);
        close;
    end
end