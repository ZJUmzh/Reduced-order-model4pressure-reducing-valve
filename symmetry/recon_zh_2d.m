function Rec = recon_zh_2d(An,phi,n_modes,u0,nx,ny)

U_new_vec=zeros(size(phi,1),size(An,1));

for i=1:n_modes
    V{i}=An(:,i).*phi(:,i)';
    U_new_vec=U_new_vec+V{i}';
end
n = size(An,1);
Rec = zeros(nx,ny,n);
for i = 1:n
        Rec(:,:,i) = reshape(U_new_vec(:,i)+u0',nx,ny);
end
end