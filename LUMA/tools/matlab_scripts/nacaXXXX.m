% Computes points for NACA MPXX aerofoil
clear
code = [0 0 1 2];
M = code(1) / 100;
P = code(2) / 10;
T = (10 * code(3) + code(4)) / 100;

% Resolution
res = 800;

% Uneven distribution to ensure even surface points
x = 1 - cos(linspace(0, pi / 2, res / 2));

% Thickness
yt = (T / 0.2) * ...
    (.2969 * sqrt(x) - .1260 * x - .3516 * x.^2 + .2843 * x.^3 -.1015 * x.^4);

% Camber line
for i = 1 : length(x)
    if x(i) < P
        yc(i) = (M / P^2) * (2 * P * x(i) - x(i)^2);
        dydx(i) = (2 * M / P^2) * (P - x(i));
    else
        yc(i) = (M / (1 - P)^2) * (1 - 2 * P + 2 * P * x(i) - x(i)^2);
        dydx(i) = (2 * M / (1 - P)^2) * (P - x(i)); 
    end
end

% Points
theta = atan(dydx);
xupper = x - yt .* sin(theta);
xlower = x + yt .* sin(theta);
yupper = yc + yt .* cos(theta);
ylower = yc - yt .* cos(theta);

% Create the arrays
out(:,1) = vertcat(xupper', xlower');
out(:,2) = vertcat(yupper', ylower');
out(:,3) = zeros(size(out,1),1);
out = unique(out,'rows');

% Write file
dlmwrite(['naca_' ...
    num2str(1000 * code(1) + 100 * code(2) + 10 * code(3) + code(4)) ...
    '.pointcloud'], out, '\t');

    
