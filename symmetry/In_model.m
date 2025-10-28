%% Boundary of the valve
clear
clc
close all

dmin1 = 3e-3;
addpath('boundary\');
boundary_all = readmatrix('boundary-100mm.txt');
boundary_xy = boundary_all(find((abs(boundary_all(:,4))<=1e-3) |...
    (abs(boundary_all(:,4))<=2e-3) & (abs(boundary_all(:,2))>=0.2125) |...
    (abs(boundary_all(:,4))<=2e-3) & (abs(boundary_all(:,2))<=0.13) |...
    (abs(boundary_all(:,4))<=2e-3) & (boundary_all(:,3)>=0.135) |...
    (abs(boundary_all(:,4))<=2e-3) & (boundary_all(:,3)<=-0.03)),:);
boundary_xy = boundary_xy(find(abs(boundary_xy(:,2))<=1),:);
% bvalve = [boundary_xy(:,2),boundary_xy(:,3)];
scatter(boundary_xy(:,2),boundary_xy(:,3),10,'filled');
hold on

hole_x = [-0.22,-0.172,-0.172,-0.165,0.165,0.172,0.172,0.22,0.22,0.172,0.172,-0.172,-0.172,-0.22,-0.22];
hole_y = [0.01,0.01,0.005,-0.002,-0.002,0.005,0.01,0.01,0.1,0.1,0.085,0.085,0.1,0.1,0.01];
%plot(hole_x,hole_y);

in_hole = inpolygon(boundary_xy(:,2),boundary_xy(:,3),hole_x,hole_y);

boundary_hole = boundary_xy(in_hole,:);
boundary_wall = boundary_xy(~in_hole,:);
% scatter(boundary_wall(:,2),boundary_wall(:,3))
hold on

bwall = [boundary_wall(:,2),boundary_wall(:,3)];

up = find(bwall(:,2)>0.05 | ((bwall(:,1)>0) & (bwall(:,2)>-0.1)));
wall_up = bwall(up,:);
wall_down = setdiff(bwall,wall_up,"rows");
wall_up = sortrows(wall_up,1);
wall_down = sortrows(wall_down,1);

valve_up_sorted = wall_up(1,:);
remain_point = wall_up(2:end,:);
potential_idx = [];
while ~isempty(remain_point)
    distance = sqrt((remain_point(:,1)-valve_up_sorted(end,1)).^2+(remain_point(:,2)-valve_up_sorted(end,2)).^2);
    potential_idx = find(distance<=dmin1);
    if length(potential_idx)>=2
        temp_idx = potential_idx;
        while ~isempty(potential_idx)
            temp_d = distance(potential_idx);
            [~,idx] = min(temp_d);
            valve_up_sorted=[valve_up_sorted;remain_point(potential_idx(idx),:)];
            potential_idx(idx) = [];
        end
        remain_point(temp_idx,:) = [];
    else
        [~,idx] = min(distance);
        valve_up_sorted = [valve_up_sorted;remain_point(idx,:)];
        remain_point(idx,:) = [];
    end
end

% plot(valve_up_sorted(:,1),valve_up_sorted(:,2));


valve_down_sorted = wall_down(1,:);
remain_point_down = wall_down(2:end,:);
while ~isempty(remain_point_down)
    distance = sqrt((remain_point_down(:,1)-valve_down_sorted(end,1)).^2+...
        (remain_point_down(:,2)-valve_down_sorted(end,2)).^2);
    potential_idx = find(distance<=dmin1);
    if length(potential_idx)>=2
        temp_idx = potential_idx;
        while ~isempty(potential_idx)
            temp_d = distance(potential_idx);
            [~,idx] = min(temp_d);
            valve_down_sorted=[valve_down_sorted;remain_point_down(potential_idx(idx),:)];
            potential_idx(idx) = [];
        end
        remain_point_down(temp_idx,:) = [];
    else
        [~,idx] = min(distance);
        valve_down_sorted = [valve_down_sorted;remain_point_down(idx,:)];
        remain_point_down(idx,:) = [];
    end
end
valve_sorted = [valve_up_sorted;flipud(valve_down_sorted);valve_up_sorted(1,:)];
plot(valve_sorted(:,1),valve_sorted(:,2));
%% hole
hole_sorted = [];
x_ini = -0.2125;
y_ini = 0.014;
for i = 1:9
    hole_sorted((i-1)*6+1,:) = [x_ini,y_ini];
    hole_sorted((i-1)*6+2,:) = [x_ini+0.0375,y_ini];
    hole_sorted((i-1)*6+3,:) = [x_ini+0.0375,y_ini+0.004];
    hole_sorted((i-1)*6+4,:) = [x_ini,y_ini+0.004];
    hole_sorted((i-1)*6+5,:) = [x_ini,y_ini];
    hole_sorted((i-1)*6+6,:) = [nan,nan];
    y_ini = y_ini+0.01;
end

x_ini = 0.175;
y_ini = 0.014;
for i = 1:9
    hole_sorted = [hole_sorted;x_ini,y_ini];
    hole_sorted = [hole_sorted;x_ini+0.0375,y_ini];
    hole_sorted = [hole_sorted;x_ini+0.0375,y_ini+0.004];
    hole_sorted = [hole_sorted;x_ini,y_ini+0.004];
    hole_sorted = [hole_sorted;x_ini,y_ini];
    hole_sorted = [hole_sorted;nan,nan];
    y_ini = y_ini+0.01;
end

x_ini = -0.165;
y_ini = -0.0005;
for i = 1:9
    hole_sorted = [hole_sorted;x_ini,y_ini];
    hole_sorted = [hole_sorted;x_ini+0.0325,y_ini];
    hole_sorted = [hole_sorted;x_ini+0.0325,y_ini+0.003];
    hole_sorted = [hole_sorted;x_ini,y_ini+0.003];
    hole_sorted = [hole_sorted;x_ini,y_ini];
    hole_sorted = [hole_sorted;nan,nan];
    y_ini = y_ini+0.01;
end

x_ini = 0.1325;
y_ini = -0.0005;
for i = 1:9
    hole_sorted = [hole_sorted;x_ini,y_ini];
    hole_sorted = [hole_sorted;x_ini+0.0325,y_ini];
    hole_sorted = [hole_sorted;x_ini+0.0325,y_ini+0.003];
    hole_sorted = [hole_sorted;x_ini,y_ini+0.003];
    hole_sorted = [hole_sorted;x_ini,y_ini];
    hole_sorted = [hole_sorted;nan,nan];
    y_ini = y_ini+0.01;
end

valve_sorted = [valve_sorted;nan,nan;hole_sorted];

%%
tic
in_mode = inpolygon(x_grid,y_grid,valve_sorted(:,1),valve_sorted(:,2));
toc
%% hole
in_temp = in_mode;
x_ini = -0.2125;
y_ini = 0.014;
for i = 1:9
    for ix = 1:size(in_temp,1)
        for iy = 1:size(in_temp,2)
            if((x_grid(ix,iy)>=x_ini & x_grid(ix,iy)<=x_ini+0.0375 & y_grid(ix,iy)>=y_ini & y_grid(ix,iy)<=y_ini+0.004)|...
                    (x_grid(ix,iy)>=-x_ini-0.0375 & x_grid(ix,iy)<=-x_ini & y_grid(ix,iy)>=y_ini & y_grid(ix,iy)<=y_ini+0.004))
                in_temp(ix,iy) = false;
            end
        end
    end
    y_ini = y_ini+0.01;
end

x_ini = -0.165;
y_ini = -0.0005;
for i = 1:9
    for ix = 1:size(in_temp,1)
        for iy = 1:size(in_temp,2)
            if((x_grid(ix,iy)>=x_ini & x_grid(ix,iy)<=x_ini+0.0325 & y_grid(ix,iy)>=y_ini & y_grid(ix,iy)<=y_ini+0.003)|...
                    (x_grid(ix,iy)>=-x_ini-0.0325 & x_grid(ix,iy)<=-x_ini & y_grid(ix,iy)>=y_ini & y_grid(ix,iy)<=y_ini+0.003))
                in_temp(ix,iy) = false;
            end
        end
    end
    y_ini = y_ini+0.01;
end
save('model.mat','in_temp')