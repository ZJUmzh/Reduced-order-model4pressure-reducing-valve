clear
clc
close all

nh = 5;
np = 2;

for ih = 1:nh
    for ip = 1:np
        tic
        origin(np*ih+ip-np).data = readmatrix([num2str(20*ih),'mm-2-',num2str(0.9+0.4*ip),'.txt']);
        [num2str(20*ih),'mm-2-',num2str(0.9+0.4*ip),'.txt']
        origin(np*ih+ip-np).data(:,1) = [];
        origin(np*ih+ip-np).data = origin(2*ih+ip-2).data(find(abs(origin(2*ih+ip-2).data(:,1))<=0.5),:);
        pause(0.000001);
        toc
    end
end

%% Downsampling
n = nh*np;
dx = 0.002;
dy = 0.001;
dz = 0.005;

o_x = origin(n).data(:,1);
o_y = origin(n).data(:,2);
o_z = origin(n).data(:,3);

nx = round((max(o_x)-min(o_x))/dx);
ny = round((max(o_y)-min(o_y))/dy);
nz = round((max(o_z)-min(o_z))/dz);

xq = linspace(min(o_x),max(o_x),nx);
yq = linspace(min(o_y),max(o_y),ny);
zq = linspace(min(o_z),max(o_z),nz);

[x_grid,y_grid,z_grid] = ndgrid(xq,yq,zq);

for i =1:n
    %vel = sqrt(origin_data(i).data(:,6).^2+origin_data(i).data(:,7).^2+origin_data(i).data(:,8).^2);
    rhou = origin(i).data(:,4).*origin(i).data(:,5);
    data_pod = zeros(nx,ny,nz,n);
    tic
    data_pod(:,:,:,i) = griddata(origin(i).data(:,1),origin(i).data(:,2),origin(i).data(:,3),...
        rhou,x_grid,y_grid,z_grid,'linear');
    data_pod(isnan(data_pod))=0;
    pause(0.00001);
    toc
end
%%
clear origin;
%% 3D-POD
clc
tic
[An,Ds,phi,U0x] = POD_3d(data_pod);
toc
%%
energy.ds = Ds;
semilogy(1:n,energy.ds/sum(energy.ds),'o-','LineWidth',1.5);
save('energy.mat','energy');

%% Energy distribution
figure(1)
semilogy(1:n,energy.ds/sum(energy.ds),'o-','LineWidth',1.5);
hold on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
yticks = get(gca, 'ytick'); % 获取当前y轴刻度值
newLabels = arrayfun(@(x) sprintf('%.2f%%', x*100), yticks, 'UniformOutput', false); % 将刻度值转换为百分比形式
set(gca, 'yticklabel', newLabels); % 设置新的y轴刻度标签
ylabel('Energy');
xlabel('k');
set(gcf, 'Units', 'centimeters', 'Position', [15,10,14,10]);
set(gcf,'Color',[1 1 1]);
%grid on
legend('Velocity','Pressure','Temperature');

figure(2)
semilogy(1:n,cumsum(energy.ds)/sum(energy.ds),'o-','LineWidth',1.5);
hold on
set(gca, 'FontName', 'Times New Roman', 'FontSize',18);
ylabel('Total Energy');
xlabel('k');
yticks = get(gca, 'ytick'); % 获取当前y轴刻度值
newLabels = arrayfun(@(x) sprintf('%.2f%%', x*100), yticks, 'UniformOutput', false); % 将刻度值转换为百分比形式
set(gca, 'yticklabel', newLabels); % 设置新的y轴刻度标签
set(gcf, 'Units', 'centimeters', 'Position', [30,10,14,10]);
set(gcf,'Color',[1 1 1]);
ylabel('Total Energy');
xlabel('k');
%grid on
legend('Velocity','Pressure','Temperature');



%% Reconstruction
nmodes = 2;%The first 4 modes are choosed
U_new = zeros(size(phiU));
tic
for i =1:nmodes
    % tic
    V{i}=An(:,i).*phiU(:,i)';
    U_new = U_new + V{i}';
    pause(0.0000001);
    % toc
end
toc
%%
nmodes = 4;
clc
tic
Rec = recon_zh(An,phiU,nmodes,U0x,nx,ny,nz);
toc
%% Validation of reconstruction
clc
tic
for i = 1:n
    u_temp = squeeze(Rec(:,:,:,i));
    pcolor(squeeze(u_temp(6,:,:)));
    shading interp
    colorbar;
    colormap;
    pause(0.1);
    q_out_pod(i,1) = trapz(zq,trapz(yq,squeeze(u_temp(6,:,:))));
end
toc