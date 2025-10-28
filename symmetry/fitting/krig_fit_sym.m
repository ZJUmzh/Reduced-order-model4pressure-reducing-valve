function An_new = krig_fit_sym(ph_new,ph,An,nmodes)

addpath('dace')
An(:,nmodes+1:end) = [];


% theta = ones(2,1)*100;
% lob = ones(2,1)*1e-1;
% upb = ones(2,1)*10000;

theta = 100;
lob = 1e-1;
upb = 1e6;

[krigmodel,krig_perf] = dacefit(ph,An,@regpoly2,@correxp,theta,lob,upb);
% [krigmodel,krig_perf] = dacefit(ph,An,@regpoly0,@correxp,theta,lob,upb);
[An_new,krig_mse] = predictor(ph_new,krigmodel);
end