clear
clc
close all

addpath('symdata\')
nh = 6;
np = 4;
n = nh*np;
%%
for ih = 1:nh
    for ip = 1:np
        tic
        sym(ih,ip).data = readmatrix(['sym-',num2str(20*ih),'mm-2-',num2str(0.9+0.2*ip),'.txt']);
        sym(ih,ip).data(:,1) = [];%序号不要
        sym(ih,ip).data(:,3) = [];
        sym(ih,ip).data = sym(ih,ip).data(find(abs(sym(ih,ip).data(:,1))<=0.5),:);
        pause(0.000001);
        toc
    end
end
save('sym_data.mat','sym');
%% Sampling
load('sym_data.mat');
dx = 0.001;
dy = 0.001;

o_x = sym(nh,np).data(:,1);
o_y = sym(nh,np).data(:,2);

nx = round((max(o_x)-min(o_x))/dx);
ny = round((max(o_y)-min(o_y))/dy);

xq = linspace(min(o_x),max(o_x),nx);
yq = linspace(min(o_y),max(o_y),ny);

[x_grid,y_grid] = ndgrid(xq,yq);
%%
samp_sym.entropy = zeros(nx,ny,nh,np);
for ih =1:nh
    for ip = 1:np
        %rhou = origin(i).data(:,4).*origin(i).data(:,5);
        tic
        % number = ih*np+ip-np
        samp_sym.entropy(:,:,ih,ip) = griddata(sym(ih,ip).data(:,1),sym(ih,ip).data(:,2),...
            sym(ih,ip).data(:,10),x_grid,y_grid,'linear');
        pause(0.00001);
        toc
    end
end
samp_sym.entropy(isnan(samp_sym.entropy)) = 0;

save('samp_sym.mat','samp_sym')
%% Read testdata
addpath('testdata\');
load('samp_test.mat');
for i = 1:4
    test_data = readmatrix([num2str(1+0.2*i),'MPa.txt']);
    test_data(:,1) = [];%序号不要
    test_data(:,3) = [];
    test_data = test_data(find(abs(test_data(:,1))<=0.5),:);
    samp_test(i).entropy= griddata(test_data(:,1),test_data(:,2),...
                test_data(:,10),x_grid,y_grid,'linear');
    samp_test(i).entropy(isnan(samp_test(i).entropy)) = 0;
end
save('samp_test.mat',"samp_test");
%% fit
load('samp_sym.mat');
load('samp_test.mat');
%% physics error
mu = 3.17e-05;
u_phy = squeeze(samp_sym.u(:,:,1,1));
v_phy = squeeze(samp_sym.v(:,:,1,1));
p_phy = squeeze(samp_sym.pressure(:,:,1,1));
rho_phy = squeeze(samp_sym.density(:,:,1,1));
[ux,uy] = gradient(u_phy);
[uxx,uxy] = gradient(ux);
[uxy,uyy] = gradient(uy);
[vx,vy] = gradient(v_phy);
[vxx,vxy] = gradient(vx);
[vxy,vyy] = gradient(vy);
[px,py] = gradient(p_phy);
% error_phy = u_phy.*ux + v_phy.*uy + px./rho_phy - mu./rho_phy.*(uxx+vyy)+1./3*(mu./rho_phy).*(uxx+vxy);
% error_phy = abs(error_phy) ./ mean(abs(u_phy.*ux + v_phy.*uy),'all')*100;
[rhoux,rhouy] = gradient(rho_phy.*u_phy);
[rhovx,rhovy] = gradient(rho_phy.*v_phy);
error_phy = (rhoux+rhovy)/ mean(abs(rhoux + rhovy),'all')*100;
error_phy(~in_temp) = nan;
%surf(error_phy(200:end-199,200:end-199));
range1 = 50:150;
range2 = 200:300;
subplot(2,1,1)
pcolor(x_grid(range1,range2),y_grid(range1,range2),error_phy(range1,range2));
shading interp
colorbar
subplot(2,1,2)
error_temp = error_phy;
error_temp(range1,range2) = 100;
pcolor(x_grid(10:end-9,10:end-9),y_grid(10:end-9,10:end-9),error_temp(10:end-9,10:end-9))
colorbar
clim([0,100])
shading interp
%%
addpath('fitting\')
clc
for ih = 1:nh
    for ip = 1:np
        ph(ih*np+ip-np,1) =  20*ih;
        ph(ih*np+ip-np,2) = 0.9+0.2*ip;
    end
end
ph_new = [100,1.2;100,1.4;100,1.6;100,1.8];

%Iwant = sqrt(samp_sym.u.^2+samp_sym.v.^2);
Iwant = samp_sym.pressure;
%Itest = sqrt(samp_test.u.^2+samp_test.v.^2);
% Itest = squeeze(Iwant(:,:,5,2));

[An,Ds,phi,U0x] = POD_2d_hp(Iwant);
% Reconstruction
clc
nmodes = 12;%99.5% of energy=16,99.9% = 21
energyplot(Ds);

An_new_krig = krig_fit_sym(ph_new,ph,An,nmodes);
rng(0);
An_new_svm = svm_fit_sym(ph_new,ph,An,nmodes);
Rec_krig = recon_zh_2d(An_new_krig,phi,nmodes,U0x,nx,ny);
Rec_svm = recon_zh_2d(An_new_svm,phi,nmodes,U0x,nx,ny);
%% POD
test_idx = 1;
Itest = samp_test(test_idx).entropy;
Itest(isnan(Itest)) = 0;
% calculate the error
%error = rmse(Rec,samp_sym.entropy,"all")./mean(samp_sym.entropy,"all")
% Figure of flow field
load('model.mat');
clc
pc1 = Itest;
pc2 = squeeze(Rec_krig(:,:,test_idx));
pc1(~in_temp) = nan;
pc2(~in_temp) = nan;
pc = abs(pc1-pc2)./mean(abs(pc1(~isnan(pc1))),'all').*100;
pc(~in_temp) = nan;
error = max(pc(10:end-9,10:end-9),[],'all')
rerror = rmse(pc1(~isnan(pc1)),pc2(~isnan(pc2)),'all')./mean(pc1(~isnan(pc1)),'all')
% for i = 1:4
%     sim_temp(:,:,i) = samp_test(i).entropy;
% end
% rerror = rmse(sim_temp,Rec_krig,'all')./mean([samp_test.entropy],'all')
%pc = pc./1000000;
%pc1(900:1100,100:250) = 0;
figure(3)
subplot(3,1,1)
% pcolor(x_grid(10:end-9,10:end-9),y_grid(10:end-9,10:end-9),pc1(10:end-9,10:end-9));
pcolor(x_grid(10:end-9,10:end-9),y_grid(10:end-9,10:end-9),pc1(10:end-9,10:end-9));
xlabel('x/m');
ylabel('y/m');
% clim([700,1000])
% xlim([-0.23,0.23]);
% ylim([-0.05,0.15]);
% axis off;
cb = colorbar;
cb.Title.String = '             Entropy / J/Kg·K';
% cb.Title.String = '             Temperature / K';
shading interp;
set(gca, 'FontName', 'Times New Roman', 'FontSize',10);
subplot(3,1,2)
pcolor(x_grid(10:end-9,10:end-9),y_grid(10:end-9,10:end-9),pc2(10:end-9,10:end-9));
%pcolor(x_grid(200:800,350:500),y_grid(200:800,350:500),pc2(200:800,350:500)./1e6);
% clim([1.2,2])
xlabel('x/m');
ylabel('y/m');
cb = colorbar;
cb.Title.String = '             Entropy / J/Kg·K';
shading interp;
set(gca, 'FontName', 'Times New Roman', 'FontSize',10);
subplot(3,1,3)
pcolor(x_grid(10:end-9,10:end-9),y_grid(10:end-9,10:end-9),pc(10:end-9,10:end-9));

colormap(othercolor('BuDRd_18'));
cb = colorbar;
shading interp;
cb.Title.String = '             Relative error / %';
% xlim([-0.5,0.5]);
% ylim([-0.4,0.3]);
% clim([0,8]);
xlabel('x/m');
ylabel('y/m');
set(gca, 'FontName', 'Times New Roman', 'FontSize',10);
set(gcf,'Color',[1 1 1]);
set(gcf, 'Units', 'centimeters', 'Position', [30,10,7,15]);
%% Physical error
nmodes = 12;
An_new_svm = svm_fit_sym(ph,ph,An,nmodes);
Rec_svm = recon_zh_2d(An_new_svm,phi,nmodes,U0x,nx,ny);

clc
test_idx = 1;
pc = Rec_svm;
pc1 = [];u = []; v = [];rho = [];
for ih = 1:nh
    for ip = 1:np
        pc1(:,:,ih*np+ip-np) = squeeze(samp_sym.pressure(:,:,ih,ip));
        u(:,:,ih*np+ip-np) = samp_sym.u(:,:,ih,ip);
        v(:,:,ih*np+ip-np) = samp_sym.v(:,:,ih,ip);
        rho(:,:,ih*np+ip-np) = samp_sym.density(:,:,ih,ip);
    end
end
%%
xp = find(x_grid(:,1)>=0.35 & x_grid(:,1)<=0.45);
yp = find(y_grid(1,:)>=-0.22 & y_grid(1,:)<=-0.12);
u = u(xp,yp,:);
v = v(xp,yp,:);
Rec = Rec_svm(xp,yp,:);
rho = rho(xp,yp,:);

dudx = gradient(u,1)./dx;
dudy = gradient(u,2)./dy;
dvdx = gradient(v,1)./dx;
dvdy = gradient(v,2)./dy;
dpdx = gradient(Rec,1)./dx;
dpdy = gradient(Rec,2)./dy;
dudxx = gradient(dudx,1)./dx;
dudyy = gradient(dudy,2)./dy;
dvdyy = gradient(dvdy,2)./dy;
dvdxy = gradient(dvdx,2)./dy;
dvdxx = gradient(dvdx,1)./dx;
dudxy = gradient(dudx,2)./dy;
%先当近似无黏流动来算
pe = (dpdx+rho.*(u.*dudx+v.*dudy))./mean(abs(rho.*(u.*dudx+v.*dudy)),"all");
pe1 = (dpdy+rho.*(u.*dvdx+v.*dvdy))./mean(abs(rho.*(u.*dudx+v.*dudy)),"all");

% r = corr2(dpdy,rho.*(u.*dvdx+v.*dvdy))
pe_mean = mean(abs(pe),'all');
pe1_mean = mean(abs(pe1),'all');

er = 0.5*pe_mean+0.5*pe1_mean;
subplot(3,1,1)
pcolor(x_grid,y_grid,pc(:,:,1))
shading interp
colorbar
subplot(3,1,2)
pcolor(x_grid(xp,yp),y_grid(xp,yp),abs(pe(:,:,1)))
shading interp
colorbar
subplot(3,1,3)
pcolor(x_grid(xp,yp),y_grid(xp,yp),abs(pe1(:,:,1)))
shading interp
colorbar
%% PIGD
clc
load('svmmodel.mat');
learnrate = 0.001;
my_ks = [];
for i = 1:length(svmmdl)
    my_ks(i) = svmmdl(i).model.Epsilon;
end

rng(0)
tic
[er,~] = physics_e(my_ks,svmmdl,ph,An,nmodes,phi,U0x,nx,ny,u,v,rho,x_grid,y_grid);
toc
%
decay = 0.01;
[ks1,ps] = mapminmax(my_ks,1,2);
for i=1:10
    tic
    [er1,~] = physics_e(my_ks.*0.9999,svmmdl,ph,An,nmodes,phi,U0x,nx,ny,u,v,rho,x_grid,y_grid);
    [er2,~] = physics_e(my_ks.*1.0001,svmmdl,ph,An,nmodes,phi,U0x,nx,ny,u,v,rho,x_grid,y_grid);
    (er1+er2)./2 
    i
    % dg = (er2-er1)./(0.002.*ks1);
    dg = (er2-er1)./(0.0002.*my_ks);
    i_learnrate = learnrate./(1+decay*i);
    % ks1 = ks1-i_learnrate.*dg;
    my_ks = my_ks - i_learnrate.*dg;
    % my_ks = mapminmax('reverse',ks1,ps);
    % my_ks(my_ks<1e-3) = 1e-3;
    pause(0.0001);
    toc
end
[er_end,svmmdl_phy] = physics_e(my_ks,svmmdl,ph,An,nmodes,phi,U0x,nx,ny,u,v,rho,x_grid,y_grid);

%% Test for Physics-informed method
% [er_end,svmmdl_phy] = physics_e(my_ks,svmmdl,ph,An,nmodes,phi,U0x,nx,ny,u,v,rho,x_grid,y_grid);
for i = 1:nmodes
    An_new_phy(:,i) = predict(svmmdl(i).model,ph_new);
end
Rec_phy = recon_zh_2d(An_new_phy,phi,nmodes,U0x,nx,ny);
%
test_idx = 1;
Itest = samp_test(test_idx).pressure;
Itest(isnan(Itest)) = 0;
% calculate the error
%error = rmse(Rec,samp_sym.entropy,"all")./mean(samp_sym.entropy,"all")
% Figure of flow field
load('model.mat');
clc
pc1 = Itest;
pc2 = squeeze(Rec_phy(:,:,test_idx));
pc1(~in_temp) = nan;
pc2(~in_temp) = nan;
pc = abs(pc1-pc2)./mean(abs(pc1(~isnan(pc1))),'all').*100;
pc(~in_temp) = nan;
error = max(pc(10:end-9,10:end-9),[],'all')
rerror = rmse(pc1(~isnan(pc1)),pc2(~isnan(pc2)),'all')./mean(pc1(~isnan(pc1)),'all')
%pc = pc./1000000;
%pc1(900:1100,100:250) = 0;
figure(3)
subplot(3,1,1)
% pcolor(x_grid(10:end-9,10:end-9),y_grid(10:end-9,10:end-9),pc1(10:end-9,10:end-9));
pcolor(x_grid(10:end-9,10:end-9),y_grid(10:end-9,10:end-9),pc1(10:end-9,10:end-9)./1e6);
xlabel('x/m');
ylabel('y/m');
% clim([700,1000])
% xlim([-0.23,0.23]);
% ylim([-0.05,0.15]);
% axis off;
cb = colorbar;
cb.Title.String = '             Pressure / MPa';
% cb.Title.String = '             Temperature / K';
shading interp;
set(gca, 'FontName', 'Times New Roman', 'FontSize',10);
subplot(3,1,2)
pcolor(x_grid(10:end-9,10:end-9),y_grid(10:end-9,10:end-9),pc2(10:end-9,10:end-9));
%pcolor(x_grid(200:800,350:500),y_grid(200:800,350:500),pc2(200:800,350:500)./1e6);
% clim([1.2,2])
xlabel('x/m');
ylabel('y/m');
cb = colorbar;
cb.Title.String = '             Pressure / MPa';
shading interp;
set(gca, 'FontName', 'Times New Roman', 'FontSize',10);
subplot(3,1,3)
pcolor(x_grid(10:end-9,10:end-9),y_grid(10:end-9,10:end-9),pc(10:end-9,10:end-9));

colormap(othercolor('BuDRd_18'));
cb = colorbar;
shading interp;
cb.Title.String = '             Relative error / %';
% xlim([-0.5,0.5]);
% ylim([-0.4,0.3]);
clim([0,20]);
xlabel('x/m');
ylabel('y/m');
set(gca, 'FontName', 'Times New Roman', 'FontSize',10);
set(gcf,'Color',[1 1 1]);
set(gcf, 'Units', 'centimeters', 'Position', [30,10,7,15]);
%%
clc
addpath('figure\');
line_y = find(y_grid==yq(394));
data_fluent = readmatrix('E:\study materials\2024\2024autumn\PRV\2-data\data4sci1.xlsx','Sheet','3.2_Flow_parameter',...
    'Range','A2:H501');
figure(4)
throttle_plot_pressure(Rec_krig,Rec_svm,x_grid,y_grid,line_y,data_fluent)
