function [er,svmmdl_phy] = physics_e(my_ks,svmmdl,ph,An,nmodes,phi,U0x,nx,ny,u,v,rho,x_grid,y_grid)

% for i = 1:nmodes
%     if i == 1
%         kfunction = 'linear';
%     else
%         kfunction = 'polynomial';
%     end
%     svmmdl(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', 'none', ...
%     'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
%         'MaxTime', 300),'KernelFunction',kfunction,'Standardize',true,'Epsilon',svmmdl(i).model.Epsilon*0.1,...
%         'BoxConstraint',svmmdl(i).model.BoxConstraints(1),'KernelScale',my_ks(i));
%     An_new(:,i) = predict(svmmdl(i).model,ph);
% end
for i = 1:nmodes
    close all
    if i == 1
    svmmdl_phy(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', 'none', ...
    'KernelFunction','linear','Standardize',true,'Epsilon',my_ks(i),...
        'BoxConstraint',svmmdl(i).model.BoxConstraints(1),'KernelScale',svmmdl(i).model.KernelParameters.Scale);
    else
        params = hyperparameters('fitrsvm',ph,An(:,i));
        params(1).Range = [1e-20,1e20];
        params(1:5) = [];
        svmmdl_phy(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', {'Standardize'}, ...
    'KernelFunction','polynomial','Standardize',true,'Epsilon',my_ks(i),...
        'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
            'MaxTime', 1,'Verbose',0,'ShowPlots',false),'BoxConstraint',svmmdl(i).model.BoxConstraints(1),...
            'KernelScale',svmmdl(i).model.KernelParameters.Scale);
    end
    An_new(:,i) = predict(svmmdl_phy(i).model,ph);
end
Rec = recon_zh_2d(An_new,phi,nmodes,U0x,nx,ny);
dx = 0.001;
dy = 0.001;

xp = find(x_grid(:,1)>=-0.1 & x_grid(:,1)<=0.1);
yp = find(y_grid(1,:)>=-0.22 & y_grid(1,:)<=-0.12);
u = u(xp,yp);
v = v(xp,yp);
Rec = Rec(xp,yp);
rho = rho(xp,yp);

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

mu=3.17e-05;
pe = (dpdx+rho.*(u.*dudx+v.*dudy)-mu.*((dudxx+dudyy)+(dudxx+dvdxy)./3))./mean(abs(rho.*(u.*dudx+v.*dudy)),"all");
pe1 = (dpdy+rho.*(u.*dvdx+v.*dvdy)-mu.*((dvdxx+dvdyy)+(dudxy+dvdyy)./3))./mean(abs(rho.*(u.*dudx+v.*dudy)),"all");

pe_mean = mean(abs(pe),'all');
pe1_mean = mean(abs(pe1),'all');

er = 0.5*pe_mean+0.5*pe1_mean;


end