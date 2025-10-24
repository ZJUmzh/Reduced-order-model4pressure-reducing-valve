function Rec = recon_zh(An,phi,n_modes,u0,nx,ny,nz)

U_new_vec=zeros(size(phi,1),size(An,1));

for i=1:n_modes
    V{i}=An(:,i).*phi(:,i)';
    U_new_vec=U_new_vec+V{i}';
end
Rec = zeros(nx,ny,nz,size(An,1));
for i = 1:size(An,1)
    Rec(:,:,:,i) = reshape(U_new_vec(:,i)+u0',nx,ny,nz);
end
end