function [] = VtkBodyMaker(bod, n, out, loopedbody)
% Generate vtk file of IBbody geometry for paraview.
% VtkBodyMaker(bod, n, out, loopedbody)
% bod is the number of bodies. n is the number of timesteps to be
% included in the movie. out is how timesteps were performed before writing
% out (see defintions.h).
% loopedbody is an array of flags to indicate whether the body is a cloased
% geometry or not.

exit_flag = false;

% Cycle through number of text files and read in positions
for c = 1:n+1
    
    if mod(c,out) == 0 % Only check for files which we know exist
        fprintf('Time %d / %d \r',c,n);
        
        for b = 0:bod-1
            try % If simulation crashes still allows movie to be built up to that point
            eval(['xyz' num2str(b) ' = csvread(''../Output/Body_' num2str(b)... 
                '_position_' num2str(c) '.out'',1,0);']) % Read in values ignoring first line
            catch
                exit_flag = true;
                break                
            end
        end
        if exit_flag == true
            break
        end
        
        
        for b = 0:bod-1
            % Create a vtk file with required header
            [fID, status] = fopen(['../Output/Body_' num2str(b) ...
                '_vtk_' num2str(c) '.vtk'],'w');
            fprintf(fID,'# vtk DataFile Version 3.0\n');
            fprintf(fID,'IBbody position\n');
            fprintf(fID,'ASCII\n');
            fprintf(fID,'DATASET POLYDATA\n');

            % Write the positions
            eval(['len = length(xyz' num2str(b) ');'])
            fprintf(fID,['POINTS ' num2str(len) ' float\n']);
            for i = 1:len
               eval(['fprintf(fID,''%8.5f %8.5f %8.5f\n'',xyz' num2str(b) ...
                   '(i,1),xyz' num2str(b) '(i,2),xyz' num2str(b) '(i,3));'])
            end

            % Write the lines
            if (loopedbody(b+1))
                lenl = len;
            else
                lenl = len-1;
            end
            cpt = 0;
            fprintf(fID,['LINES ' num2str(lenl) ' ' num2str(3*lenl) '\n']);
            for i = 1:lenl
                fprintf(fID,['2 ' num2str(cpt) ' ' num2str(cpt+1) '\n']);
                cpt = cpt + 1;
            end
        status = fclose(fID);
        end
        
    end
end