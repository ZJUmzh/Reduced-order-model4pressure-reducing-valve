function An_new = svm_fit_PIGD(ph_new,ph,An,nmodes,phi,nx,ny,U0x)
% load('E:\study materials\2024\2024autumn\PRV\matlab_code\newdata\fitting\svmmdl_new.mat')

load('svmmodel.mat');
for i = 1:nmodes
    An_new(:,i) = predict(svmmdl(i).model,ph_new);
end
Rec_svm = recon_zh_2d(An_new,phi,nmodes,U0x,nx,ny);


end