function pts = FilledCuboidCreator(filename, resolution, sizeXYZ)
numpts = round(sizeXYZ * resolution);
pts = zeros(numpts(1) * numpts(2) * numpts(3), 3);
c = 0;
for i = linspace(0, sizeXYZ(1), numpts(1))
    for j = linspace(0, sizeXYZ(2), numpts(2))
        for k = linspace(0, sizeXYZ(3), numpts(3))
            c = c + 1;
            pts(c, 1) = i;
            pts(c, 2) = j;
            pts(c, 3) = k;
        end
    end
end
dlmwrite(filename, pts, 'delimiter', '\t');
plot3(pts(:,1), pts(:,2),pts(:,3),'.','MarkerSize',1)
axis equal
grid on