% STL to Point Cloud converter
% Faces and Vertices are 3 column matrices of vertex numbers
% and vertex positions respectively. Resolution is how many points
% per length the surface should be discretised into using the object
% bounding box dimension in the x-direction as the length.
% Requires stlread from The Mathworks to work.
%
% Copyright Adrian Harwood, The University of Manchester, UK
function PC = StlToPc(infile, outfile, resolution)

% Use STL read first
[faces, vertices] = stlread(infile);

% Start counters
count = 0;
count_discrete = 0;
count_notdiscrete = 0;
point_added = false;

% Length
model_length = max(vertices(:,1)) - min(vertices(:,1));
dh = model_length / resolution;

% Display progress
fprintf('Discretising %d faces...\n', size(faces,1));
charcount = 1;

% Loop over faces
for i = 1:size(faces,1)
    
    % Progress
    fprintf(1, repmat('\b',1,charcount));   % Delete line before
    charcount = fprintf('\rDiscretising face %d / %d. %f %% complete.', i, size(faces,1), i / size(faces,1) * 100);

    % Create vectors for each vertex of triangle
    v1 = [vertices(faces(i,1),1); vertices(faces(i,1),2); vertices(faces(i,1),3)];
    v2 = [vertices(faces(i,2),1); vertices(faces(i,2),2); vertices(faces(i,2),3)];
    v3 = [vertices(faces(i,3),1); vertices(faces(i,3),2); vertices(faces(i,3),3)];

    % Translate to origin
    v10 = v1 - v1;
    v20 = v2 - v1;
    v30 = v3 - v1;

    % Rotate to reference triangle
    u = v20;            % First side
    w = cross(u, v30);  % Get normal to both sides
    
    if (norm(u) < eps || norm(w) < eps)
      % Just add a point at the centre for now
      untrans_c = [...
        (v1(1) + v2(1) + v3(1)) / 3; ...
        (v1(2) + v2(2) + v3(2)) / 3; ...
        (v1(3) + v2(3) + v3(3)) / 3 ...
        ];
        count = count + 1;
        count_notdiscrete = count_notdiscrete + 1;
        PC(count,:) = untrans_c';
        continue;
    end
    
    % Normalise both to get unit vectors
    U = u / norm(u);
    W = w / norm(w);

    % Final unit normal to previous normal and first side
    V = cross(U, W);

    % Rotation matrix assembly
    R = [U'; V'; W'];
    v1f = R * v10;
    v2f = R * v20;
    v3f = R * v30;

    % Find centre of triangle
    c = [...
        (v1f(1) + v2f(1) + v3f(1)) / 3; ...
        (v1f(2) + v2f(2) + v3f(2)) / 3; ...
        (v1f(3) + v2f(3) + v3f(3)) / 3 ...
        ];
    
    % Area
    AB = v2f - v1f;
    AC = v3f - v1f;
    Area = norm(cross(AB,AC)) / 2;

    % If almost no area then just add centre point
    if (Area < eps)
      bkmap = Ri * c;   % Reverse rotation
      bkmap = bkmap + v1;       % Reverse shift
      count = count + 1
      count_notdiscrete = count_notdiscrete + 1;
      PC(count,:) = bkmap';
      continue;      
    end
    
    % Inverse for later
    Ri = inv(R);

    % Seed points
    xmin = min([v1f(1) v2f(1) v3f(1)]);
    xmax = max([v1f(1) v2f(1) v3f(1)]);
    ymin = min([v1f(2) v2f(2) v3f(2)]);
    ymax = max([v1f(2) v2f(2) v3f(2)]);

    % Loop over block of points
    for x = xmin:dh:xmax
        for y = ymin:dh:ymax
        
          % Reference coords of point
          vp = [x; y; 0];
          AP = vp - v1f;
                  
          % Compute barycentric coordinates of point
          d00 = dot(AB, AB);
          d01 = dot(AB, AC);
          d11 = dot(AC, AC);
          d20 = dot(AP, AB);
          d21 = dot(AP, AC);
          denom = d00 * d11 - d01 * d01;
          if (abs(denom) < eps)
            continue;
          end
          alpha = (d11 * d20 - d01 * d21) / denom;
          beta = (d00 * d21 - d01 * d20) / denom;
          gamma = 1 - alpha - beta;
                
          % Test for within the triangular plane and not on edges
          if (...
              alpha > 0 && alpha < 1 &&...
              beta > 0 && beta < 1 &&...
              gamma > 0 && gamma < 1 ...
              )
             
            % Reverse mapping and store
            bkmap = Ri * vp;
            bkmap = bkmap + v1;
            count = count + 1;
            
            if (~point_added)
              count_discrete = count_discrete + 1;
              point_added = true;
            end
            PC(count,:) = bkmap';
            
          end
            
        end
    end

    % If too small for resolution just add centre
    if (~point_added)
      bkmap = Ri * c;
      bkmap = bkmap + v1;
      count = count + 1;
      count_notdiscrete = count_notdiscrete + 1;
      PC(count,:) = bkmap';
    end

    point_added = false;   
    
end

PC = [PC(:,1) PC(:,3) PC(:,2)];

fprintf('\n');
disp(['dh = ' num2str(dh)]);
disp(['Number of faces discretised = ' num2str(count_discrete)]);
disp(['Number of faces not discretised = ' num2str(count_notdiscrete)]);
disp(['Total number of points = ' num2str(count)]);

% Swap columns
tmp = PC(:,2);
PC(:,2) = PC(:,3);
PC(:,3) = tmp;

figure
plot3(PC(:,1),PC(:,2),PC(:,3),'.','MarkerSize',1);
axis equal
axis tight
grid on

dlmwrite(outfile,PC,'Delimiter','\t');

end