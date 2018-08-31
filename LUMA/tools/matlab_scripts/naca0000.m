clear
close all
res = 800;
AR = 2;
t = 0.12;
x = 1 - cos(linspace(0,pi / 2,res));
x_store = x;
y_upper = 5*t*(.2969*sqrt(x) - .1260*x - .3516*x.^2 + .2843*x.^3 -.1015*x.^4);
y_lower = -y_upper;
y_lower = fliplr(y_lower);
y_lower(1) = [];
y_lower(end) = [];
y = [y_upper y_lower];
x = x';
x_lower = flipud(x);
x_lower(1) = [];
x_lower(end) = [];
x = [x; x_lower];
y = y';
plot(x,y,'-x'), axis equal, axis([0 1 -t/2 t/2]);
z = linspace(0, AR, AR * res);
for zz = 1 : length(z);
    for p = 1 : length(x)
        z_ex(p + (zz - 1) * length(x)) = z(zz);
    end
end
z_ex = z_ex';
x_ex = repmat(x,length(z),1);
y_ex = repmat(y,length(z),1);

% Create a rectangle of points with same spacing as x and y
num_y = round(t * res);
if (num_y < 3)
    num_y = 3;
end
[x_ends, y_ends] = meshgrid(x_store, linspace(min(y), max(y), num_y ));
y_upper_mat = repmat(y_upper, size(y_ends, 1), 1);
y_lower_mat = repmat(-y_upper, size(y_ends, 1), 1);
logical_end = logical(y_ends > y_lower_mat) .* logical(y_ends < y_upper_mat);
re_log = logical(reshape(logical_end, size(logical_end,1) * size(logical_end,2), 1));
x_end_allow = reshape(x_ends, length(re_log), 1);
x_end_allow = x_end_allow(re_log);
y_end_allow = reshape(y_ends, length(re_log), 1);
y_end_allow = y_end_allow(re_log);
z_end_allow = [min(z) * ones(length(x_end_allow), 1); max(z) * ones(length(x_end_allow), 1)];
x_end_allow = [x_end_allow; x_end_allow];
y_end_allow = [y_end_allow; y_end_allow];


% Tag onto vectors
x_ex = [x_ex; x_end_allow];
y_ex = [y_ex; y_end_allow];
z_ex = [z_ex; z_end_allow];

figure, plot3(x_ex, y_ex, z_ex, '.', 'MarkerSize',2), axis equal
dlmwrite(['naca00' num2str(t * 100) '_2D.pointcloud'], [x y ones(size(x))], 'Delimiter', '\t');
dlmwrite(['naca00' num2str(t * 100) '_3D.pointcloud'], [x_ex y_ex z_ex], 'Delimiter', '\t');