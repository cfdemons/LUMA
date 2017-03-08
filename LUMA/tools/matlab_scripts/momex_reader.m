clear
close all

out_every = 60000;
last_time = 60000;
last_rank = 7;
scale = 0.2;

% Read the CSV file set and stitch together
for t = out_every : out_every : last_time

  % New figure
  figure
    
  for rank = 0 : last_rank
    data = csvread(['momex_debug_' num2str(t) '_Rnk' num2str(rank) '.csv'],1,0);
    
    % If no points on this rank then continue
    if (isempty(data))
      continue
    else
      disp(['Time = ' num2str(t) ', Rank' num2str(rank) ', Data count = ' num2str(size(data,1))]);
    end
    
    % Sort data
    pos = data(:,1:3);
    Fx = data(:,4:3:size(data,2));
    Fy = data(:,5:3:size(data,2));
    Fz = data(:,6:3:size(data,2));
    ForceX = sum(Fx,2);
    ForceY = sum(Fy,2);
    ForceZ = sum(Fz,2);

    % Plot contributing Forces from each site
    h1 = quiver3(pos(:,1),pos(:,2),pos(:,3),ForceX,ForceY,ForceZ,scale,'Autoscale','off');
    hU = get(h1,'UData') ;
    hV = get(h1,'VData') ;
    set(h1,'UData',scale*hU,'VData',scale*hV)
    grid on
    hold on

    % Plot contributing markers
    plot3(pos(:,1),pos(:,2),pos(:,3),'rs','MarkerSize',14);

    view(2)

    % Flip arrowheads
    L = get(h1,'children');
    D = get(L(2),{'xdata','ydata'}); % Arrowhead data 
    C = get(L(3),{'xdata','ydata'}); % End point data

    % Get X points excluding the NaNs
    NPx = D{1}(~isnan(D{1}));
    MPx = C{1}(~isnan(C{1}));

    % Get Y points excluding the NaNs
    NPy = D{2}(~isnan(D{2}));
    MPy = C{2}(~isnan(C{2}));

    % Debugging
    %plot3(NPx,NPy,repmat(pos(1,3),length(NPx),1),'ko','MarkerSize',10)
    %plot3(MPx(1:2:length(MPx)),MPy(1:2:length(MPx)),repmat(pos(1,3),length(MPx)/2,1),'kx','MarkerSize',10)
    %plot3(MPx(2:2:length(MPx)),MPy(2:2:length(MPx)),repmat(pos(1,3),length(MPx)/2,1),'k*','MarkerSize',10)

    % Get length of arrows based on "M" data
    LenX = MPx(2:2:length(MPx)) - MPx(1:2:length(MPx));
    LenY = MPy(2:2:length(MPx)) - MPy(1:2:length(MPx));

    % Shift coordinates of arrows
    Op = zeros(size(MPx));
    Op(1:2:length(MPx)) = LenX;
    Op(2:2:length(MPx)) = LenX;
    MPx = MPx - Op;

    Op = zeros(size(MPy));
    Op(1:2:length(MPy)) = LenY;
    Op(2:2:length(MPy)) = LenY;
    MPy = MPy - Op;

    % Shift arrowheads
    Op = zeros(size(NPx));
    Op(1:3:length(NPx)) = LenX;
    Op(2:3:length(NPx)) = LenX;
    Op(3:3:length(NPx)) = LenX;
    NPx = NPx - Op;

    Op = zeros(size(NPy));
    Op(1:3:length(NPy)) = LenY;
    Op(2:3:length(NPy)) = LenY;
    Op(3:3:length(NPy)) = LenY;
    NPy = NPy - Op;

    %plot3(NPx,NPy,repmat(pos(1,3),length(NPx),1),'go','MarkerSize',10)

    % Put NaNs back in
    MPx = reshape(MPx,2,length(MPx)/2);
    MPx(3,:) = NaN;
    MPx = reshape(MPx,1,size(MPx,2)*3);

    MPy = reshape(MPy,2,length(MPy)/2);
    MPy(3,:) = NaN;
    MPy = reshape(MPy,1,size(MPy,2)*3);

    NPx = reshape(NPx,3,length(NPx)/3);
    NPx(4,:) = NaN;
    NPx = reshape(NPx,1,size(NPx,2)*4);

    NPy = reshape(NPy,3,length(NPy)/3);
    NPy(4,:) = NaN;
    NPy = reshape(NPy,1,size(NPy,2)*4);

    set(L(2),{'xdata','ydata'},{NPx,NPy}) % Update arrowheads
    set(L(3),{'xdata','ydata'},{MPx,MPy}) % Update arrows

    axis equal    
    axis tight

    end
end